#include <cmath>
#include <iostream>
#include <fstream>
#include <filesystem>
#include "../headers/matrix.h"
#include "../headers/grid.h"
#include "../headers/HLLC.h"
#include "../headers/slopelim.h"

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
 Build: clang++ -std=c++17 -O2 -I./headers sources/HLLC.cpp
        sources/slope_limiters.cpp sources/grid_setup.cpp sources/solver.cpp
        -o builds/solver
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
            Cell c = grid(i,j); //c contains [rho, rho u, rho v, rho e_tot]
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

//!!! CHANGE THIS FUNCTION ACCORDING TO YOUR PROBLEM !!!
void init_cond(Grid& grid, double xmin, double ymin, double dx, double dy, double gamma) {
  //modifies the initial grid elements == Cells (@ t=t0) as specified DO NOT FORGET TO CHANGE THIS!
  size_t nx = grid.rows() - 4;
  size_t ny = grid.cols() - 4;
  for (size_t i = 2; i < nx+2; i++) {
      for (size_t j = 2; j < ny+2; j++) {
          double x = xmin + (i-2)*dx;  //physical coordinates
          double y = ymin + (j-2)*dy;
          Cell c;
          //example
          if (x < 0.5) {
                c[0] = 1.0;              // rho
                c[1] = 0.0;              // rho*u
                c[2] = 0.0;              // rho*v
                c[3] = 1.0/(gamma-1.0);  // E = p/(gamma-1)
            } else {
                c[0] = 0.125;
                c[1] = 0.0;
                c[2] = 0.0;
                c[3] = 0.1/(gamma-1.0);
            }
          grid(i,j) = c;
      }
  }
}


BCType readBCType() { //BC type input reader + safety check if BC type is unknown
  std::string s;
  std::cin >> s;
  if      (s == "Open")      return BCType::Open;
  else if (s == "Closed")    return BCType::Closed;
  else if (s == "Periodic")  return BCType::Periodic;
  else if (s == "Dirichlet") return BCType::Dirichlet;
  else throw std::invalid_argument("Unknown BC type: " + s);
}

BoundaryCondition readBC() { //just save the BC type as a BoundaryCondition so that GridBC can be used
  BoundaryCondition bc;
  bc.type = readBCType();
  if (bc.type == BCType::Dirichlet) { //if Dirichlet -> the constant state vec at wall must be specified
      std::cout << "Enter fixed state (rho rho*u rho*v E): ";
      bc.Q_fixed[0] = 0; bc.Q_fixed[1] = 0; bc.Q_fixed[2] = 0; bc.Q_fixed[3] = 0; //if nothing -> 0 vector
      std::cin >> bc.Q_fixed[0] >> bc.Q_fixed[1] >> bc.Q_fixed[2] >> bc.Q_fixed[3];
  }
  return bc;
}


Grid timestep(const Grid& grid, double dx, double dy, double dt, double gamma){
  size_t nx = grid.rows();
  size_t ny = grid.cols();
  Grid grid_new(nx, ny);

  for (size_t i = 2; i<nx-2; i++){
    for(size_t j = 2; j<ny-2; j++){
      Vector Qi = CellToVec(grid(i,j)); //bunch of Q_{ij} vals needed for the slope calculations 
      Vector Qim1 = CellToVec(grid(i-1,j));
      Vector Qim2 = CellToVec(grid(i-2,j));
      Vector Qip1 = CellToVec(grid(i+1,j));
      Vector Qip2 = CellToVec(grid(i+2,j));
      Vector Qjm1 = CellToVec(grid(i,j-1));
      Vector Qjm2 = CellToVec(grid(i,j-2));
      Vector Qjp1 = CellToVec(grid(i,j+1));
      Vector Qjp2 = CellToVec(grid(i,j+2));

      //slopes in the x direction
      Vector sigma_im1_x = sigma_superbee(Qim2, Qim1, Qi,   dx);
      Vector sigma_i_x   = sigma_superbee(Qim1, Qi,   Qip1, dx);
      Vector sigma_ip1_x = sigma_superbee(Qi,   Qip1, Qip2, dx);
      //slopes in the y direction
      Vector sigma_jm1_y = sigma_superbee(Qjm2, Qjm1, Qi,   dy);
      Vector sigma_j_y   = sigma_superbee(Qjm1, Qi,   Qjp1, dy);
      Vector sigma_jp1_y = sigma_superbee(Qi,   Qjp1, Qjp2, dy);
      
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
      grid_new(i,j) = VecToCell(Q_new); //convert Q_new to Cell and then append it to the new grid
    }
  }
  return grid_new; //return new grid after dt timestep
}




//CLAUDE AI: "write grid into file every timestep with performance in mind"
//A: "For writeGrid, binary is the right choice for animation — fast to write,
//compact, and easy to read in Python with numpy.fromfile. 
//A simple format: write a header (nx, ny, t) then all cell data sequentially." 
//My idea: txt file but turns out binary is better...
//I'm not used to file writing so I just used AI here :/

void setupOutputDir(const std::string& dir) { //this was added later so that the user can specify
  if (!std::filesystem::exists(dir)) {      //a folder name, the simulation result will be saved there
      std::filesystem::create_directory(dir);
      std::cout << "Created output directory: " << dir << std::endl;
  }
}

void writeGrid(const Grid& grid, double t, int n, const std::string& outdir) {
  // filename: frame_0000.bin, frame_0001.bin etc.
  char filename[128];
  std::snprintf(filename, sizeof(filename), "%s/frame_%04d.bin", outdir.c_str(), n);
  
  std::ofstream f(filename, std::ios::binary);
  if (!f) throw std::runtime_error("Could not open file for writing");

  size_t nx = grid.rows() - 4;  // physical cells only
  size_t ny = grid.cols() - 4;

  // header
  f.write(reinterpret_cast<const char*>(&nx), sizeof(size_t));
  f.write(reinterpret_cast<const char*>(&ny), sizeof(size_t));
  f.write(reinterpret_cast<const char*>(&t),  sizeof(double));

  // data: physical cells only, skip ghost layers
  for (size_t i = 2; i < nx+2; i++) {
    for (size_t j = 2; j < ny+2; j++) {
        f.write(reinterpret_cast<const char*>(&grid(i,j).q), 4*sizeof(double));
    }
  }
}


int main(){

  std::cout << "!! Remember to change init_cond() in solver.cpp for your problem !!\n";
  double Nx, Ny, Nt, t0, tf, xmin, xmax, ymin, ymax, dx, dy, dt;
  
  //parameters of the grid as input
  std::cout << "Grid parameters\n";
  std::cout << "Nx, Ny, Nt, t0, tf, xmin, xmax, ymin, ymax:\n";
  std::cin >> Nx >> Ny >> Nt >> t0 >> tf >> xmin >> xmax >> ymin >> ymax;
  std::cout << std::endl;

  //some safeguards - no negative Nx, Ny, Nt and so on
  if (Nx <=0 || Ny <= 0 || Nt <= 0) throw std::invalid_argument("Nx, Ny, Nt cannot be negative!");
  if (xmax <= xmin || ymax <= ymin || tf <= t0){
    throw std::invalid_argument("xmax <= xmin or ymax <= ymin or tf <= t0!");
  }
  //TO DO? another idea would be chosing the max of the input xmax and xmin as xmax and the other is xmin

  //step-size calculation
  dt = (tf-t0)/Nt;
  dx = (xmax-xmin)/Nx;
  dy = (ymax-ymin)/Ny;
  
  double gamma;
  std::cout << "Gamma:";
  std::cin >> gamma;
  std::cout << std::endl;

  GridBC bc;
  std::cout << "Enter BC for left wall   (Open/Closed/Periodic/Dirichlet): "; bc.left   = readBC();
  std::cout << "Enter BC for right wall  (Open/Closed/Periodic/Dirichlet): "; bc.right  = readBC();
  std::cout << "Enter BC for top wall    (Open/Closed/Periodic/Dirichlet): "; bc.top    = readBC();
  std::cout << "Enter BC for bottom wall (Open/Closed/Periodic/Dirichlet): "; bc.bottom = readBC();
  validateBC(bc);

  //Nx,Ny were intentionally declared as doubles so that dt, dx, dy 
  //could be calculated, here I convert them to size_t unsigned type.
  //Conversion should go without issues as Nx, Ny are expected to be 
  //positive vals, if they are not, the safeguards throw error beforehand
  size_t nx = (size_t) (Nx + 4);
  size_t ny = (size_t) (Ny + 4);
  Grid grid(nx, ny); //ghost points included
  init_cond(grid, xmin, ymin, dx, dy, gamma); //initial condition applied to the grid

  //cfl condition check - either new dt or continue with CFL
  double cfl = cfl_dt(grid, dx, dy, gamma)
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
        if (dt > cfl) //
            std::cout << "Still above CFL (" << cfl << "), negative density/pressure/energy is expected...\n";
    }
  }
  
  std::string outdir;
  std::cout << "Output directory name: ";
  std::cin >> outdir;
  setupOutputDir(outdir);

  double t = t0;
  for (int n = 0; n<Nt; n++){
    // dt = std::min((tf-t0)/Nt, cfl_dt(grid, dx, dy, gamma));
    t = t0 + n*dt;
    writeGrid(grid, t, n, outdir);
    applyBC(grid, bc);
    grid = timestep(grid, dx, dy, dt, gamma);
    // t+=dt;
  }
  
  

  return 0;
}
