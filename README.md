# Hydro-2D
2D HLLC hydrodynamics code. Solves the 2D compressible Euler equations in conservative form on a uniform Cartesian grid using the HLLC Riemann solver with Superbee slope-limited interpolation at cell walls.

## Builds
Includes the solver.exe file. To compile solve.cpp I recommend using the following command in the terminal window: clang++ -std=c++17 -O2 -I./headers sources/HLLC.cpp sources/slope_limiters.cpp sources/grid_setup.cpp sources/solver.cpp -o builds/solver 

Feel free to compile it however You want to.

## Headers
1. **grid.h:** the header file of grid_setup.cpp found in sources, also includes the Cell struct and Grid class. The idea behind the grid setup was using AoS method, a Cell holds the 4 vals of the state vector - So in memory it looks like $\{$rho_00, rhou_00, rhov_00, E_00,  rho_01, rhou_01, rhov_01, E_01, ...$\}$. The other functions enforce the boundary conditions provided by the user on the walls + some safety checks
2. **HLLC.h:** header file of HLLC.cpp found in sources.
3. **matrix.h:** contains numpy array like logic for matrices and vectors, the matrix class was NEVER used as I wanted to experiment with the Cell & Grid AoS method. VecOps namespace is used almost everywhere and it just means some operator overrides for vector<double>-s
4. **slopelim.h** header file of slope_limiters.cpp

## Scripts
Contains animate.py, which creates an animation of the solution using the grid values at every timestep saved in binary files for more efficiency. Also produces a 1D snippet at y = ymax/2 of the solution at the last timestep - this was important for the Sod shock tube test, disable it if not needed.

## Sources
1. **grid_setup.cpp:** contains functions used for the enforcing of the boundary conditions and some safety checks especially for the periodic BC.
2. **HLLC.cpp:** is the main logic of the solver, contains the methods used for the HLLC flux calculation at cell walls.
3. **slope_limiters:** only contains the superbee slope limiter, logic behind the extrapolation of the state vector till the cell walls.
4. **solver.cpp:** is everything put together, sets grid and initial conditions + BC, calculates CFL condition then the Euler time-steps. To calculate the time-step it first uses the superbee slope limiter to interpolate the state vector values at the cell walls, then uses the HLLC method to calculate the HLLC flux. After this, we get a new grid with $Q^{n+1} = Q^{n} + dt/dx (F^{n}_{i+1/2,j}-F^{n}_{i-1/2,j) + dt/dy(F^{n}_{i,j+1/2}-F^{n}_{i,j-1/2})$ The new grid values are then saved in a binary file.
