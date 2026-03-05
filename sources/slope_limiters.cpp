#include <cmath>
#include "../headers/matrix.h"
#include <stdexcept>

/*
 Superbee slope limiter script
 Bertalan Szuchovszky 26.02.2026

 Based on FSM lecture notes, surely this can be found online aswell
 The main function is sigma superbee:
    -> calculates superbee slope based on Q_{i-1}, Q_{i}, Q_{i+1} and dx or dy
    -> will be used to estimate Q_L and Q_R values at cell interfaces
  Every other function is a helper function
*/

using namespace VecOps;

// minmod(a,b) = b if a*b>0, |a| > |b|
//             = a if a*b>0, |a| <= |b|
//             = 0 if a*b<0
double minmod(double a, double b){
  if (a*b<0.0) return 0.0;
  else if (std::abs(a)>std::abs(b)) return b;
  else return a;
  }
// minmod(a,b) = a if a*b>0, |a| > |b|
//             = b if a*b>0, |a| <= |b|
//             = 0 if a*b<0
double maxmod(double a, double b){
  if (a*b<0.0) return 0.0;
  else if (std::abs(a)>std::abs(b)) return a;
  else return b;
}

//superbee: maxmod(s1,s2) with s1,s2 being minmod slopes with the proper arguments
double superbee_one(double qim1, double qi, double qip1, double dx){ //dx or dy depending on direction
  double s1 = minmod((qi-qim1)/dx, 2*(qip1-qi)/dx);
  double s2 = minmod(2*(qi-qim1)/dx, (qip1-qi)/dx);
  return maxmod(s1,s2);
}

//Q = [rho, rho u, rho v, rho e_tot], i index means i-th gridpoint
//vectorial version, calculate superbee slope for every element of state vector
Vector sigma_superbee(const Vector& Qim1, const Vector& Qi, const Vector& Qip1, double dk){
  if (Qim1.size()!=Qi.size() || Qim1.size()!=Qip1.size() || Qi.size()!=Qip1.size()){
    throw std::invalid_argument("Vector size mismatch");
  }
  double s0 = superbee_one(Qim1[0], Qi[0], Qip1[0], dk); //k = x or y depending on direction
  double s1 = superbee_one(Qim1[1], Qi[1], Qip1[1], dk); //in practice:
  double s2 = superbee_one(Qim1[2], Qi[2], Qip1[2], dk); //sigma_x = sigma_superbee(Qim1,Qi,Qip1,dx)
  double s3 = superbee_one(Qim1[3], Qi[3], Qip1[3], dk); //sigma_y = sigma_superbee(Qim1,Qi,Qip1,dy)
  return {s0, s1, s2, s3};
}
