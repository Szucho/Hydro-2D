#include <cmath>
#include "../headers/matrix.h"
#include <stdexcept>


/*
  HLLC method for 2D Euler HD equation
  Bertalan Szuchovszky 26.02.2026

  For details read Toro - Riemann Solvers and Numerical Methods for Fluid Dynamics (3rd ed)
  I implemented chapter 10 using my matrix.h header for numpy like arrays.
  This code serves to calculate the hllc flux with given Q state vector and gamma adiabatic index.
  The Toro book uses SK instead of lambdas, but to stay consistent with the FSM lecture notation
  I will use lambda_{m}, lambda_{0} and lambda_{+} instead.
*/

using namespace VecOps; //vector operation namespace for std::vector<double> in matrix.h

struct Primitive_Vals{
  double rho, u, v, p, cs;
};

static Primitive_Vals QtoPrim(const Vector& Q, double gamma){
  if (Q.size() != 4) throw std::invalid_argument("State vector must be 4D");

  Primitive_Vals q;
  q.rho = Q[0];
  q.u = Q[1]/q.rho; //rho u / rho
  q.v = Q[2]/q.rho; //rho v / rho
  //Q[3] = rho e_tot = rho(e_th + 0.5(u^2 + v^2)), e_th = 1/(gamma-1)k_BT/mu m_p
  //p = rho k_BT/mu m_p ideal gas law -> p = (gamma-1)*e_th = (gamma-1)rho(e_tot-0.5(u^2+v^2)) 
  q.p = (gamma-1.0)*(Q[3]-q.rho * 0.5*(q.u*q.u + q.v*q.v));
  q.cs = std::sqrt(gamma * q.p/q.rho);

  if (q.rho <0.0||q.p<0.0) throw std::invalid_argument("Density or pressure is negative");
  return q;
}


//flux in the x direction - Toro chapter 10
static Vector xFlux(const Vector&Q, double gamma){
  Primitive_Vals q = QtoPrim(Q, gamma);
  double e_tot = Q[3]/q.rho;
  double h_tot = e_tot + q.p/q.rho;
  //Fx = (rho u, rho u^2 + p, rho u v, rho u h_tot)
  return{
    q.rho*q.u,
    q.rho*q.u*q.u + q.p,
    q.rho*q.v*q.u,
    q.rho*q.u*h_tot
  };
}

//flux in the y direction
static Vector yFlux(const Vector&Q, double gamma){
  Primitive_Vals q = QtoPrim(Q, gamma);
  double e_tot = Q[3]/q.rho;
  double h_tot = e_tot + q.p/q.rho;
  //Fy = (rho u, rho u v, rho v^2 + p, rho v h_tot)
  return{
    q.rho*q.u,
    q.rho*q.u*q.v,
    q.rho*q.v*q.v + q.p,
    q.rho*q.v*h_tot
  };
}


struct lambdas{
  double lm, l0, lp; //lambda_{-}, lambda_{0}, lambda_{+} eigenvals of the Jacobi matrix
};


static lambdas WavecX(const Primitive_Vals& L, const Primitive_Vals& R){ //L:left, R:right
  //Toro 10.48
  double SL = std::min(L.u-L.cs, R.u-R.cs); //lambda_{-} = u-cs but let it be the min of these
  double SR = std::max(L.u+L.cs, R.u+R.cs); //lambda_{+} = u+cs and the max of these 

  double Sstar = (R.p-L.p + L.rho*L.u*(SL-L.u) - R.rho*R.u*(SR-R.u))/(L.rho*(SL-L.u)-R.rho*(SR-R.u)); //Toro eq. 10.36
  return {SL, Sstar, SR};
}


static Vector Qtilde(const Vector&Q, double gamma, double SK, double Sstar){
  //Q is Q_K with K=L || K=R => q is also q_K
  Primitive_Vals q = QtoPrim(Q, gamma);
  double mul = q.rho*(SK-q.u)/(SK-Sstar);
  double E = Q[3];
  //Toro 10.39
  return{
    mul,
    mul*Sstar,
    mul*q.v,
    mul*(E/q.rho + (Sstar-q.u)*(Sstar + q.p/(q.rho*(SK-q.u))))
  };
}


Vector FluxhllcX(const Vector& QL, const Vector& QR, double gamma){
  //left and right primitive vals
  Primitive_Vals L = QtoPrim(QL, gamma);
  Primitive_Vals R = QtoPrim(QR, gamma);

  lambdas w = WavecX(L, R);

  //Toro 10.26
  if (0.0<=w.lm){
    return xFlux(QL, gamma);
  }
  
  if (w.lp<=0.0){
    return xFlux(QR, gamma);
  }

  if (w.lm<=0.0 && 0.0<=w.l0){
    Vector QtL = Qtilde(QL, gamma, w.lm, w.l0); //Toro 10.27
    return xFlux(QL, gamma) + w.lm*(QtL-QL);
  }

  Vector QtR = Qtilde(QR, gamma, w.lp, w.l0);
  return xFlux(QR, gamma) + w.lp*(QtR-QR); //Toro 10.29
}


//trick to get the flux in y direction instead of writing the whole thing again
static Vector rotateQxy(const Vector& Q) {
    //[rho, rho*u, rho*v, E]  ->  [rho, rho*v, rho*u, E]
    return {Q[0], Q[2], Q[1], Q[3]};
}

Vector FluxhllcY(const Vector& QL, const Vector& QR, double gamma) {
    //rotate so that y becomes the normal direction
    Vector QL_rot = rotateQxy(QL);
    Vector QR_rot = rotateQxy(QR);

    //solve as an x direction flux but now with rotated Q
    Vector Frot = FluxhllcX(QL_rot, QR_rot, gamma);

    //rotate flux back:  G = [rho*v, rho*u*v, rho*v^2+p, (E+p)*v]
    return {Frot[0], Frot[2], Frot[1], Frot[3]};
}
