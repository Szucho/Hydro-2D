#include <cmath>
#include <iomanip>
#include <ios>
#include <iostream>
#include <fstream>
#include <filesystem>
#include <stdexcept>
#include "../headers/matrix.h"
#include "../headers/grid.h"
#include "../headers/HLLC.h"
#include "../headers/slopelim.h"
#include "../headers/io.h"

using namespace VecOps;

/*
 2D Hydrodynamics HLLC Euler Equation solver
 Bertalan Szuchovszky 03.03.2026

 Solves the 2D compressible Euler equations in conservative form
 on a uniform Cartesian grid using the HLLC approximate Riemann solver (see HLLC.cpp)
 with Superbee slope-limited interpolation at cell walls (see slope_limiters.cpp).

 Equation to be solved: dQ/dt + dF/dQ * div Q = 0 with dF/dQ = J Jacobi matrix 
 State vector: Q = [rho, rho*u, rho*v, rho*e_tot]
               e_tot = e_th + 0.5*(u^2+v^2),  e_th = 1/(gamma-1) * k_BT/(mu*m_p) - ideal gas law

 Spatial discretization: 2nd order finite volume (~MUSCL-Hancock, approximate Riemann solver)
 Time integration: Explicit Euler
 => 2nd order FV in space & 1st order in time
 Boundary conditions: Open, Closed, Periodic, Dirichlet (user specified)

 Usage: set init_cond() for your problem (DON'T FORGET THIS), then build and run.
 Build: clang++ -std=c++17 -O3 -march=native -I./headers sources/*.cpp -o builds/szuchydro2d.exe

 No CMake yet as I can't be bothered to write one
 Validated on Sod shock tube (Toro, Chapter 4)
*/


//FIRST: CFL condition -> need to check if dt > dt_{cfl}
double cfl_dt(const Grid& grid, double dx, double dy, double gamma, double CFL_max=0.5) {
    size_t nx = grid.rows() - 4;
    size_t ny = grid.cols() - 4;
    double smax_x = 0.0, smax_y = 0.0; //s characteristics
    for (size_t i = 2; i < nx+2; i++) {
        for (size_t j = 2; j < ny+2; j++) {
            Cell c = grid.getCell(i,j); //c contains [rho, rho u, rho v, rho e_tot]
            double rho = c[0];
            double u   = c[1]/rho;
            double v   = c[2]/rho;
            double p   = (gamma-1.0)*(c[3] - 0.5*rho*(u*u+v*v));
            double cs  = std::sqrt(gamma*p/rho); //adiabatic soundspeed
            smax_x = std::max(smax_x, std::abs(u)+cs);
            smax_y = std::max(smax_y, std::abs(v)+cs);
        }
    }
    return CFL_max / (smax_x/dx + smax_y/dy);
}

//HLLC step of a single grid point Q_{ij}
Vector HLLC_step(const Vector& Qi, 
                 const Vector& FluxXiph, 
                 const Vector& FluxYiph, 
                 const Vector& FluxXimh, 
                 const Vector& FluxYimh,
                 double dx, double dy, double dt){
  Vector Qi_new; //return the new state vector at ij gridpoint
  Qi_new = Qi - dt/dx*(FluxXiph - FluxXimh) - dt/dy*(FluxYiph-FluxYimh); //Euler timestep
  return Qi_new;
}




Grid timestep(const Grid& grid, double dx, double dy, double dt, double gamma){
  size_t nx = grid.rows();
  size_t ny = grid.cols();
  Grid grid_new(nx, ny);

  for (size_t i = 2; i<nx-2; i++){
    for(size_t j = 2; j<ny-2; j++){
      Vector Qi   = CellToVec(grid.getCell(i,j)); //bunch of Q_{ij} vals needed for the slope calculations 
      Vector Qim1 = CellToVec(grid.getCell(i-1,j));
      Vector Qim2 = CellToVec(grid.getCell(i-2,j));
      Vector Qip1 = CellToVec(grid.getCell(i+1,j));
      Vector Qip2 = CellToVec(grid.getCell(i+2,j));
      Vector Qjm1 = CellToVec(grid.getCell(i,j-1));
      Vector Qjm2 = CellToVec(grid.getCell(i,j-2));
      Vector Qjp1 = CellToVec(grid.getCell(i,j+1));
      Vector Qjp2 = CellToVec(grid.getCell(i,j+2));

      //slopes in x direction
      Vector sigma_im1_x, sigma_i_x, sigma_ip1_x;
      sigma_im1_x = sigma_minmod(Qim2, Qim1, Qi,   dx);
      sigma_i_x   = sigma_minmod(Qim1, Qi,   Qip1, dx);
      sigma_ip1_x = sigma_minmod(Qi,   Qip1, Qip2, dx);

      //slopes in the y direction
      Vector sigma_jm1_y, sigma_j_y, sigma_jp1_y;
      sigma_jm1_y = sigma_minmod(Qjm2, Qjm1, Qi,   dy);
      sigma_j_y   = sigma_minmod(Qjm1, Qi,   Qjp1, dy);
      sigma_jp1_y = sigma_minmod(Qi,   Qjp1, Qjp2, dy);
      
      //x interfaces
      Vector QL_imh = Qim1 + 0.5*dx * sigma_im1_x;  //left  state at i-1/2
      Vector QR_imh = Qi   - 0.5*dx * sigma_i_x;    //right state at i-1/2
      Vector QL_iph = Qi   + 0.5*dx * sigma_i_x;    //left  state at i+1/2
      Vector QR_iph = Qip1 - 0.5*dx * sigma_ip1_x;  //right state at i+1/2

      //y interfaces
      Vector QL_jmh = Qjm1 + 0.5*dy * sigma_jm1_y;  //left  state at j-1/2
      Vector QR_jmh = Qi   - 0.5*dy * sigma_j_y;    //right state at j-1/2
      Vector QL_jph = Qi   + 0.5*dy * sigma_j_y;    //left  state at j+1/2
      Vector QR_jph = Qjp1 - 0.5*dy * sigma_jp1_y;  //right state at j+1/2

      //fluxes -> HLLC method (see HLLC.cpp, Toro)
      Vector FXimh = FluxhllcX(QL_imh, QR_imh, gamma);
      Vector FXiph = FluxhllcX(QL_iph, QR_iph, gamma);
      Vector FYjmh = FluxhllcY(QL_jmh, QR_jmh, gamma);
      Vector FYjph = FluxhllcY(QL_jph, QR_jph, gamma); 
       
      //apply HLLC timestep at gridcell Q_{ij}
      Vector Q_new = HLLC_step(Qi, FXiph, FYjph, FXimh, FYjmh, dx, dy, dt);
      grid_new.setCell(i, j, VecToCell(Q_new)); //convert Q_new to Cell and then append it to the new grid
    }
  }
  return grid_new; //return new grid after dt timestep
}


int main(int argc, char*argv[]){
  using std::chrono::high_resolution_clock;
  using std::chrono::duration_cast;
  using std::chrono::duration;
  using std::chrono::milliseconds;
  auto t1 = high_resolution_clock::now();

  std::string paramfile = (argc > 1) ? argv[1] : "params.txt";
  SimParams par = readParams(paramfile);
  double cs2 = par.gas_p0/par.gas_rho0;

  double Nx = par.Nx, Ny = par.Ny;
  double Nt = par.Nt;
  double t0 = par.t0, tf = par.tf;
  double xmin = par.xmin, xmax = par.xmax;
  double ymin = par.ymin, ymax = par.ymax;
  double gamma = par.gamma;

  //step-size calculation
  double dt = (tf-t0)/Nt;
  double dx = (xmax-xmin)/Nx;
  double dy = (ymax-ymin)/Ny;
  
  GridBC bc;
  bc.left.type   = parseBC(par.bc_left);
  bc.right.type  = parseBC(par.bc_right);
  bc.top.type    = parseBC(par.bc_top);
  bc.bottom.type = parseBC(par.bc_bottom);
  validateBC(bc);

  //Nx,Ny were intentionally declared as doubles so that dt, dx, dy 
  //could be calculated, here I convert them to size_t unsigned type.
  //Conversion should go without issues as Nx, Ny are expected to be 
  //positive vals, if they are not, the safeguards throw error beforehand
  size_t nx = (size_t) (Nx + 4);
  size_t ny = (size_t) (Ny + 4);
  Grid grid(nx, ny); //ghost points included
  init_cond(grid, par);
  applyBC(grid, bc);

  //cfl condition check - either new dt or continue with CFL
  double cfl = cfl_dt(grid, dx, dy, gamma);
  if (dt > cfl){
    std::cout << "Warning: dt=" << dt << " exceeds CFL limit=" << cfl << "\n";
    std::cout << "Continue with CFL timestep? (y/n): ";
    std::string ans; 
    std::cin >> ans;
    if (ans == "y") {
        dt = cfl;
        Nt = (int)std::ceil((tf - t0) / dt);  //recalculate Nt based on CFL condition
        std::cout << "Using dt=" << dt << ", Nt=" << Nt << "\n";
    } else {
        std::cout << "Enter new Nt: "; //new Nt
        std::cin >> Nt;
        dt = (tf - t0) / Nt; //new dt based on new Nt
        if (dt > cfl)
            std::cout << "Still above CFL (" << cfl << "), negative density/pressure/energy is expected...\n";
    }
  }
  
  setupOutputDir(par.outdir);
  writeMetadata(par.outdir, par);
  auto gridfile = openGridFile(par.outdir, (size_t)Nx, (size_t)Ny);
  double t = t0;

  double dt_save = (tf - t0) / 100.0; // E.g., save 100 evenly spaced frames
  double t_next_save = t0;

  int print_interval = (int)Nt/10;
  if (print_interval==0) print_interval = 1;
  
  int n=0;
  try{
    while (t<tf){
      //Dynamic timestep selection to prevent CFL violation as flow accelerates
      double dt_cfl = cfl_dt(grid, dx, dy, gamma, 0.4); // Use safe CFL ~ 0.4
      dt = std::min(dt_cfl, tf - t);

      if (n%print_interval == 0 || n == (int)Nt-1){
        double progress = (double)n/((int)Nt-1)*100.0;
        std::cout << "Step: " << n << " / " << (int)Nt
                  << " [" <<std::fixed << std::setprecision(2) << progress << "%] "
                  << "t = " << t << std::endl;
      }

      if (t >= t_next_save) {
        writeFrame(gridfile, grid, t, (size_t)Nx, (size_t)Ny);
        gridfile.flush();
        t_next_save += dt_save;
        // std::cout << "Saved frame at t = " << t << std::endl;
      }

      grid = timestep(grid, dx, dy, dt, gamma);
      applyBC(grid, bc);
      n += 1;
      t += dt;
    }

    writeFrame(gridfile, grid, t, (size_t)Nx, (size_t)Ny);
    gridfile.flush();
    std::cout << "Simulation finished successfully at t = " << t << std::endl;

  } catch (const std::invalid_argument& e){
    std::cerr << "[CRASH] at step: "<< n << ", t = " << t << std::endl;
    throw;
  }
  
  
  auto t2 = high_resolution_clock::now();
  auto ms_int = duration_cast<milliseconds>(t2 - t1);
  duration<double, std::milli> ms_double = t2 - t1;
  std::cout << ms_int.count() << "ms\n";
  std::cout << ms_double.count() << "ms\n";
  return 0;
}
