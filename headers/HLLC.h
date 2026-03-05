#ifndef HLLC_H
#define HLLC_H

#include "matrix.h"

/*
 HLLC method for 2D Euler HD equation header file
 Bertalan Szuchovszky 26.02.2026
 
 state vector Q = [rho, rho*u, rho*v, rho*e_tot] 
 with e_tot = e_th + 0.5*|v|^2 & e_th = 1/(gamma-1) k_BT/mu m_p thermal energy

 We only need the hllc flux from the HLLC.cpp file
 FluxhllcX : numerical flux at an interface in x direction
 FluxhllcY : numerical flux at an interface in y direction

 Input:
  QL, QR: left and right state vectors at the interface!!! USE SLOPE LIMITERS TO GET THEM!
  gamma: (f+2)/f adiabatic index
 Returns:
  HLLC flux vectors in x direction or y direction
*/

Vector xFlux(const Vector& Q, double gamma);
Vector yFlux(const Vector& Q, double gamma);
Vector FluxhllcX(const Vector& QL, const Vector& QR, double gamma);
Vector FluxhllcY(const Vector& QL, const Vector& QR, double gamma);

#endif // HLLC_H
