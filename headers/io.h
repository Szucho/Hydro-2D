#pragma once
#include <string>
#include <fstream>
#include <sstream>
#include <stdexcept>
#include <map>
#include <iomanip>
#include <iostream>
#include <cstdlib>
#include <cstdint>
#include "grid.h"
#include <filesystem>

constexpr double PI_IO = 3.141592653589793;


/*
 Input / Output header
 Bertalan Szuchovszky 20.04.2026.

 Only active for one instance, reads a params.txt file containing initial values of particle, gas, BC-s, ...
 Hanldes file writing into a specified folder (in params.txt), creates
  -> meta.txt containing metadata (grid params & physical constants i.e gamma)
  -> grid.bin containing all state vals on every gridpoint at every timestep
  -> particle.bin containing the position, velocity,... of particle
 There are some initial conditions available 
 !!! DO NOT FORGET TO CHANGE IT FOR A SPECIFIC INITIAL CONDITION !!!

*/

//everything the simulation needs to run
struct SimParams {
  //grid
  double Nx, Ny, Nt;
  double t0, tf;
  double xmin, xmax, ymin, ymax;
  //physics
  double gamma;
  //boundary conditions
  std::string bc_left, bc_right, bc_top, bc_bottom;
  //gas initial condition
  std::string gas_profile;
  double gas_rho0, gas_p0, gas_u0, gas_v0; 
  //output
  std::string outdir;
};


//parameter file reader
inline SimParams readParams(const std::string& filename) {
  std::ifstream f(filename);
  if (!f) throw std::runtime_error("Cannot open param file: " + filename);

  std::map<std::string, std::string> kv;
  std::string line;
  while (std::getline(f, line)) {
    auto comment = line.find('#');
    if (comment != std::string::npos) line = line.substr(0, comment);
    if (line.find_first_not_of(" \t\r\n") == std::string::npos) continue;
    std::istringstream ss(line);
    std::string key, val;
    if (ss >> key >> val) kv[key] = val;
  }

  auto get = [&](const std::string& k) -> std::string {
    auto it = kv.find(k);
    if (it == kv.end()) throw std::runtime_error("Missing parameter: " + k);
    return it->second;
  };
  auto getd = [&](const std::string& k) {
    try { return std::stod(get(k)); }
    catch (const std::invalid_argument&) {
      throw std::runtime_error("Parameter '" + k + "' is not a valid number");
    }
  };

  SimParams p;

  //grid
  p.Nx   = getd("Nx");    p.Ny  = getd("Ny");   p.Nt  = getd("Nt");
  p.t0   = getd("t0");    p.tf  = getd("tf");
  p.xmin = getd("xmin");  p.xmax = getd("xmax");
  p.ymin = getd("ymin");  p.ymax = getd("ymax");

  //physics
  p.gamma = getd("gamma");

  //BC
  p.bc_left   = get("bc_left");
  p.bc_right  = get("bc_right");
  p.bc_top    = get("bc_top");
  p.bc_bottom = get("bc_bottom");

  //gas IC
  p.gas_profile = get("gas_profile");
  p.gas_rho0    = getd("gas_rho0");
  p.gas_p0      = getd("gas_p0");
  p.gas_u0      = getd("gas_u0");
  p.gas_v0      = getd("gas_v0");


  //output
  p.outdir = get("outdir");

  //sanity checks
  if (p.Nx <= 0 || p.Ny <= 0 || p.Nt <= 0)
    throw std::invalid_argument("Nx, Ny, Nt must be positive");
  if (p.xmax <= p.xmin || p.ymax <= p.ymin)
    throw std::invalid_argument("xmax <= xmin or ymax <= ymin");
  if (p.tf <= p.t0)
    throw std::invalid_argument("tf must be greater than t0");
  if (p.gas_rho0 <= 0 || p.gas_p0 <= 0)
    throw std::invalid_argument("gas_rho0 and gas_p0 must be positive");
  return p;
}


//BC string -> BC type
inline BCType parseBC(const std::string& s) {
  if      (s == "Open")      return BCType::Open;
  else if (s == "Closed")    return BCType::Closed;
  else if (s == "Periodic")  return BCType::Periodic;
  else if (s == "Dirichlet") return BCType::Dirichlet;
  else throw std::invalid_argument("Unknown BC type: " + s);
}


//gas initial conditions
inline void init_cond(Grid& grid, const SimParams& par) {
  size_t nx = grid.rows() - 4;
  size_t ny = grid.cols() - 4;
  double dx = (par.xmax - par.xmin) / par.Nx;

  if (par.gas_profile == "uniform") {
    for (size_t i = 2; i < nx+2; i++) {
      for (size_t j = 2; j < ny+2; j++) {
        Cell c = grid.getCell(i,j);
        c[0] = par.gas_rho0;
        c[1] = par.gas_rho0 * par.gas_u0;
        c[2] = par.gas_rho0 * par.gas_v0;
        c[3] = par.gas_p0 / (par.gamma - 1.0)
             + 0.5 * par.gas_rho0 * (par.gas_u0*par.gas_u0 + par.gas_v0*par.gas_v0);
        grid.setCell(i,j,c);
      }
    }
  } else if (par.gas_profile == "Sod_x"){
    //Sod shock tube vals:
    const double rho_L = par.gas_rho0,     p_L = par.gas_p0,     u_L = 0.0, v_L = 0.0;
    const double rho_R = par.gas_rho0/8.0, p_R = 0.1*par.gas_p0, u_R = 0.0, v_R = 0.0;
    const double x_mid = 0.5 * (par.xmin + par.xmax);

    for (size_t i = 2; i < nx + 2; i++) {
      //cell-center x coordinate (accounting for 2 ghost layers)
      double x = par.xmin + (static_cast<double>(i - 2) + 0.5) * dx;

      //select Left or Right state based on x location
      double rho = (x < x_mid) ? rho_L : rho_R;
      double p   = (x < x_mid) ? p_L   : p_R;
      double u   = (x < x_mid) ? u_L   : u_R;
      double v   = (x < x_mid) ? v_L   : v_R;

      //compute total energy per unit volume: E = p/(gamma-1) + 0.5*rho*(u^2 + v^2)
      double E = p / (par.gamma - 1.0) + 0.5 * rho * (u * u + v * v);

      for (size_t j = 2; j < ny + 2; j++) {
        Cell c = grid.getCell(i, j);
        c[0] = rho;
        c[1] = rho * u;
        c[2] = rho * v;
        c[3] = E;
        grid.setCell(i, j, c);
      }
    } 
  } else {
    throw std::runtime_error("Unknown gas_profile: '" + par.gas_profile + "'. Available: uniform, gaussian_z");
  }
}



//output directory
inline void setupOutputDir(const std::string& dir) {
  std::error_code ec;
  // create_directories creates the full path and doesn't error if it exists
  std::filesystem::create_directories(dir, ec); 
  
  if (ec) {
    throw std::runtime_error("Could not create output directory: " + dir + " - " + ec.message());
  }
  std::cout << "Output directory: " << dir << "\n";
}


//metadata file txt containing some important values of the simulation
inline void writeMetadata(const std::string& outdir, const SimParams& par) {
  std::string path = outdir + "/meta.txt";
  std::ofstream f(path);
  if (!f) throw std::runtime_error("Cannot open meta.txt: " + path);
  f << std::setprecision(17) << std::scientific;
  f << "xmin "  << par.xmin  << "\n"
    << "xmax "  << par.xmax  << "\n"
    << "ymin "  << par.ymin  << "\n"
    << "ymax "  << par.ymax  << "\n"
    << "gamma " << par.gamma << "\n"
    << "Nx "    << (int)par.Nx    << "\n"
    << "Ny "    << (int)par.Ny    << "\n"
    << "cs_iso " << par.gas_p0/par.gas_rho0 <<"\n";
}


//grid output, single binary file, all frames sequential suggested by Claude
inline std::ofstream openGridFile(const std::string& outdir, size_t nx, size_t ny) {
  std::string path = outdir + "/grid.bin";
  std::ofstream f(path, std::ios::binary);
  if (!f) throw std::runtime_error("Cannot open grid output file: " + path);
  uint64_t nx64 = nx; uint64_t ny64 = ny;
  f.write(reinterpret_cast<const char*>(&nx64), sizeof(uint64_t));
  f.write(reinterpret_cast<const char*>(&ny64), sizeof(uint64_t));
  return f;
}

inline void writeFrame(std::ofstream& f, const Grid& grid, double t, size_t nx, size_t ny) {
  f.write(reinterpret_cast<const char*>(&t), sizeof(double));
  for (size_t i = 2; i < nx+2; i++)
    for (size_t j = 2; j < ny+2; j++)
      f.write(reinterpret_cast<const char*>(grid.cell(i,j)), 4 * sizeof(double));
}
