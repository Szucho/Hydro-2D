#ifndef GRID_H
#define GRID_H

#include <vector>
#include <stdexcept>
#include "matrix.h"

/*
 Grid BC application for 2D Grids
 Bertalan Szuchovszky 02.03.2026

 Cell: cell point struct, at each cell point we have a Q state vector, Cell[0] = rho, Cell[1] = rho u ...
 Grid: the whole grid, every element of the grid is a cell Grid(i,j)[0] = rho, Grid(i,j)[1] = rho u ...
 BCType: Boundary condition type - just names, the solver just needs the names
 BoundaryCondition: sets the BCType, for Dirichlet we need Cell Q_fixed constant state vector
 GridBC: BC type at each grid wall

 The other 2 functions are defined in grid_setup.cpp
*/



//at each cell we have a 4D state vector Q
struct Cell {
    double q[4];
    double& operator[](int c)       { return q[c]; } //access the state vector elements
    double  operator[](int c) const { return q[c]; }
};

//Grid consists of Cells
class Grid {
private:
    std::vector<Cell> data_;
    size_t rows_, cols_;
public:
    Grid(size_t rows, size_t cols) //give it a size
        : data_(rows * cols), rows_(rows), cols_(cols) {}
    Cell& operator()(size_t i, size_t j)             { return data_[i*cols_ + j]; } //access elements
    const Cell& operator()(size_t i, size_t j) const { return data_[i*cols_ + j]; }
    size_t rows() const { return rows_; } //access row/col size
    size_t cols() const { return cols_; }
};

enum class BCType {Open, Closed, Periodic, Dirichlet}; //BC types available

struct BoundaryCondition {
    BCType type;
    Cell Q_fixed;  //Dirichlet
};

struct GridBC {
    BoundaryCondition left, right, top, bottom; //BC on all walls
};



inline Vector CellToVec(const Cell& c) {
    return {c[0], c[1], c[2], c[3]};
}
inline Cell VecToCell(const Vector& v) {
    Cell c; c[0]=v[0]; c[1]=v[1]; c[2]=v[2]; c[3]=v[3];
    return c;
}


//defined in grid_setup.cpp
void validateBC(const GridBC& bc); //validateBC checks if Periodic BC was set correctly
void applyBC(Grid& grid, const GridBC& bc); //applies BC type at chosen wall

#endif // GRID_H
