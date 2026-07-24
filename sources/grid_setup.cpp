#include "../headers/grid.h"
#include <cmath>
#include <stdexcept>

/*
 Grid BC application for 2D Grids
 Bertalan Szuchovszky 02.03.2026

 Grid, BCType, Cell, GridBC, BoundaryCondition are defined in the grid.h header file

 validateBC: checks if the Periodic BC was set on a pair of opposing walls - if not, error message
  Input:
    GridBC
  Output:
    nothing, throws error message if periodic BC was set incorrectly

 appyBC: applies chosen BCType on a specific wall allowing for different BC at different walls
  -> Open:      what goes out the wall vanishes (Neumann BC: derivative = 0 at walls)
  -> Closed:    what collides with the wall bounces back (normal momentum flipping)
  -> Periodic:  what goes out one wall comes back on the opposing side
  -> Dirichlet: given constant state vector at walls
  Input:
    Grid & GridBC - the whole grid and the boundary conditions
  Output:
    nothing, modifies the state vector at grid walls depending on BC type
*/


void validateBC(const GridBC& bc) { //If periodic, the other wall needs to be set to periodic BC aswell
    if ((bc.left.type == BCType::Periodic) != (bc.right.type == BCType::Periodic))
        throw std::invalid_argument("Periodic BC must be applied to both left and right walls");
    if ((bc.bottom.type == BCType::Periodic) != (bc.top.type == BCType::Periodic))
        throw std::invalid_argument("Periodic BC must be applied to both bottom and top walls");
}





void applyBC(Grid& grid, const GridBC& bc) {
  validateBC(bc); //check if Periodic BC was used correctly
  size_t Nx = grid.rows() - 4; //ghost points not included
  size_t Ny = grid.cols() - 4;
  //every wall has a BC -> switch case for every wall depending on BCType
  //then just apply the BC at the walls

  //left wall
  for (size_t j = 0; j < grid.cols(); ++j) {
    switch (bc.left.type) {
      case BCType::Open:{
        grid.copyCell(0,j, 2,j);
        grid.copyCell(1,j, 2,j);
        break;
      }
      case BCType::Closed:
        grid.copyCell(1,j, 2,j); grid(1,j,1) *=-1.0; //flip normal momentum rho*u
        grid.copyCell(0,j, 3,j); grid(0,j,1) *=-1.0;
        break;
      case BCType::Periodic:
        grid.copyCell(0,j, Nx,  j); //goes out this side enters at the other side
        grid.copyCell(1,j, Nx+1,j);
        break;
      case BCType::Dirichlet:
        grid.setCell(0,j, bc.left.Q_fixed); //const BC value
        grid.setCell(1,j, bc.left.Q_fixed);
        break;
    }
  }
  //right wall
  for (size_t j = 0; j < grid.cols(); ++j) {
    switch (bc.right.type) {
      case BCType::Open:{
        grid.copyCell(Nx+2, j, Nx+1, j);
        grid.copyCell(Nx+3, j, Nx+1, j);
        break;
      }
      case BCType::Closed:
        grid.copyCell(Nx+2,j, Nx+1,j); grid(Nx+2,j,1) *= -1.0;
        grid.copyCell(Nx+3,j, Nx,  j); grid(Nx+3,j,1) *= -1.0;
        break;
      case BCType::Periodic:
        grid.copyCell(Nx+2,j, 2,j);
        grid.copyCell(Nx+3,j, 3,j);
        break;
      case BCType::Dirichlet:
        grid.setCell(Nx+2,j, bc.right.Q_fixed);
        grid.setCell(Nx+3,j, bc.right.Q_fixed);
        break;
    }
  }
  //bottom wall
  for (size_t i = 0; i < grid.rows(); ++i) {
    switch (bc.bottom.type) {
      case BCType::Open:
        grid.copyCell(i,0, i,2);
        grid.copyCell(i,1, i,2);
        break;
      case BCType::Closed:
        grid.copyCell(i,1, i,2); grid(i,1,2) *= -1.0;
        grid.copyCell(i,0, i,3); grid(i,0,2) *= -1.0;
        break;
      case BCType::Periodic:
        grid.copyCell(i,0, i,Ny);
        grid.copyCell(i,1, i,Ny+1);
        break;
      case BCType::Dirichlet:
        grid.setCell(i,0, bc.bottom.Q_fixed);
        grid.setCell(i,1, bc.bottom.Q_fixed);
        break;
    }
  }
  //top wall
  for (size_t i = 0; i < grid.rows(); ++i) {
    switch (bc.top.type) {
      case BCType::Open:
        grid.copyCell(i,Ny+2, i,Ny+1);
        grid.copyCell(i,Ny+3, i,Ny+1);
        break;
      case BCType::Closed:
        grid.copyCell(i,Ny+2, i,Ny+1); grid(i,Ny+2,2) *= -1.0;
        grid.copyCell(i,Ny+3, i,Ny  ); grid(i,Ny+3,2) *= -1.0;
        break;
      case BCType::Periodic:
        grid.copyCell(i,Ny+2, i,2);
        grid.copyCell(i,Ny+3, i,3);
        break;
      case BCType::Dirichlet:
        grid.setCell(i,Ny+2, bc.top.Q_fixed);
        grid.setCell(i,Ny+3, bc.top.Q_fixed);
        break;
    }
  }
}
