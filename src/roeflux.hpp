
#ifndef ROE_FLUX_HPP_
#define ROE_FLUX_HPP_

#include <math.h>
#include <limits>
#include <cstdio>
#include <cstdlib>
#include <cstring>
#include <sys/time.h>
#include <iostream>
#include <fstream>

template<typename flux_t,typename state_t, typename normal_t>
void roeflux_kernel( flux_t  & F, const state_t & qL,
		     const state_t & qR, const normal_t n, double g)
{
// Computes the flux for the shallow water equations  
// INPUTS:
//    qL: conservative state vector in left cell
//    qR: conservative state vector in right cell
//    n: normal pointing from the left cell to the right cell

// OUTPUTS:
//  F   : the flux out of the left cell (into the right cell)
  double es = 1.e-30;
  double hL = qL[0];
  double uL = qL[1]/(hL + es);
  double vL = qL[2]/(hL + es);
  double unL = uL*n[0] + vL*n[1];

  double pL = 0.5*g*std::pow(hL,2.);
  double FL[ 3 ];
  FL[0] = hL*unL;
  FL[1] = qL[1]*unL + pL*n[0];
  FL[2] = qL[2]*unL + pL*n[1];

  double hR = qR[0];
  double uR = qR[1]/(hR + es);
  double vR = qR[2]/(hR + es);
  double unR = uR*n[0] + vR*n[1];
  double pR = 0.5*g*std::pow(hR,2.);
  double FR[3];
  FR[0] = hR*unR;
  FR[1] = qR[1]*unR + pR*n[0];
  FR[2] = qR[2]*unR + pR*n[1];

  // rho average
  double hm = 0.5*(hL + hR);
  double um = (unL*std::pow(hL,0.5) + unR*std::pow(hR,0.5) )/( std::pow(hL,0.5) + std::pow(hR,0.5)  + es);
  // eigenvalues
  double smax = abs(um) + abs(std::pow(g*hm, 0.5));
  // flux assembly
  F[0]    = 0.5*(FL[0]+FR[0])-0.5*smax*(qR[0] - qL[0]);
  F[1]    = 0.5*(FL[1]+FR[1])-0.5*smax*(qR[1] - qL[1]);
  F[2]    = 0.5*(FL[2]+FR[2])-0.5*smax*(qR[2] - qL[2]);

}

template<typename jac_t,typename state_t, typename normal_t>
void roeflux_jacobian( jac_t  & JL , jac_t & JR, const state_t & qL,
		     const state_t & qR, const normal_t n, double g){

// Computes the flux Jacobian for the shallow water equations  
// INPUTS:
//    qL: conservative state vector in left cell
//    qR: conservative state vector in right cell
//    n: normal pointing from the left cell to the right cell

// OUTPUTS:
//  F   : the flux out of the left cell (into the right cell)

  double es = 1.e-30;

  double hL = qL[0];
  double uL = qL[1]/(hL + es);
  double vL = qL[2]/(hL + es);
  double unL = uL*n[0] + vL*n[1];

  double hR = qR[0];
  double uR = qR[1]/(hR + es);
  double vR = qR[2]/(hR + es);
  double unR = uR*n[0] + vR*n[1];

  // rho average
  double hm = 0.5*(hL + hR);
  double um = (unL*std::pow(hL,0.5) + unR*std::pow(hR,0.5) )/( std::pow(hL,0.5) + std::pow(hR,0.5)  + es);
  // eigenvalues
  double smax = abs(um) + abs(std::pow(g*hm, 0.5));


  double termL = (n[0]*qL[1] + n[1]*qL[2]) / std::pow(qL[0],2.);

  double termR = (n[0]*qR[1] + n[1]*qR[2]) / std::pow(qR[0],2.);

  double hL_sqrt = std::pow(hL_sqrt,0.5);
  double hR_sqrt = std::pow(hR_sqrt,0.5);

  double hsqrt_un = hL_sqrt*unL + hR_sqrt*unR;
  double dsmaxL[3];
  double dsmaxR[3];

  dsmaxL[0] = - std::abs( hsqrt_un ) / (2.*hL_sqrt*std::pow(hL_sqrt + hR_sqrt,2.) ) + 
              (0.5*unL / hL_sqrt - hL_sqrt*termL )*hsqrt_un  / ( (hL_sqrt + hR_sqrt)*std::abs( hsqrt_un ) ) + 
              g/(std::pow(2.,3./2.)*std::pow(g*(hL + hR),0.5) );
  dsmaxL[1] = n[0]*hsqrt_un / ( hL_sqrt*(hL_sqrt + hR_sqrt) * std::abs( hsqrt_un) );
  dsmaxL[2] = n[1]*hsqrt_un / ( hL_sqrt*(hL_sqrt + hR_sqrt) * std::abs( hsqrt_un) );


  dsmaxR[0] = - std::abs( hsqrt_un ) / (2.*hR_sqrt*std::pow(hL_sqrt + hR_sqrt,2.) ) + 
              (0.5*unR / hR_sqrt - hR_sqrt*termR )*hsqrt_un  / ( (hL_sqrt + hR_sqrt)*std::abs( hsqrt_un ) ) + 
              g/(std::pow(2.,3./2.)*std::pow(g*(hL + hR),0.5) );
  dsmaxR[1] = n[0]*hsqrt_un / ( hR_sqrt*(hL_sqrt + hR_sqrt) * std::abs( hsqrt_un) );
  dsmaxR[2] = n[1]*hsqrt_un / ( hR_sqrt*(hL_sqrt + hR_sqrt) * std::abs( hsqrt_un) );

  // jacobian w.r.p to the left state
  JL[0,0] = -0.5*dsmaxL[0]*(qR[0] - qL[0]) + 0.5*(n[0]*uL + n[1]*vL - qL[0]*termL)  + 0.5*smax;
  JL[0,1] = 0.5*n[0] - 0.5*dsmaxL[1]*(qR[0] - qL[0]);
  JL[0,2] = 0.5*n[1] - 0.5*dsmaxL[2]*(qR[0] - qL[0]);

  JL[1,0] = 0.5*(g*n[0]*qL[0] - qL[1]*termL) - 0.5*dsmaxL[0]*(qR[1] - qL[1]);
  JL[1,1] = n[0]*uL  + 0.5*n[1]*vL + 0.5*smax -0.5*dsmaxL[1]*(qR[1] - qL[1]);
  JL[1,2] = 0.5*n[1]*uL  -0.5*dsmaxL[2]*(qR[1] - qL[1]);

  JL[2,0] = 0.5*(g*n[1]*qL[0] - qL[2]*termL) - 0.5*dsmaxL[0]*(qR[2] - qL[2]);
  JL[2,2] = 0.5*n[0]*vL  -0.5*dsmaxL[1]*(qR[2] - qL[2]);
  JL[2,1] = n[1]*vL + 0.5*n[0]*uL + 0.5*smax -0.5*dsmaxL[2]*(qR[2] - qL[2]);

   // jacobian w.r.p to the right state


  JR[0,0] = -0.5*dsmaxR[0]*(qR[0] - qL[0]) + 0.5*(n[0]*uR + n[1]*vR - qR[0]*termR)  - 0.5*smax;
  JR[0,1] = 0.5*n[0] - 0.5*dsmaxR[1]*(qR[0] - qL[0]);
  JR[0,2] = 0.5*n[1] - 0.5*dsmaxR[2]*(qR[0] - qL[0]);

  JR[1,0] = 0.5*(g*n[0]*qR[0] - qR[1]*termR) - 0.5*dsmaxR[0]*(qR[1] - qL[1]);
  JR[1,1] = n[0]*uR + 0.5*n[1]*vR - 0.5*smax - 0.5*dsmaxR[1]*(qR[1] - qL[1]);
  JR[1,2] = 0.5*n[1]*uR -0.5*dsmaxR[2]*(qR[1] - qL[1]);

  JR[2,0] = 0.5*(g*n[1]*qR[0]  - qR[2]*termR) - 0.5*dsmaxR[0]*(qR[2] - qL[2]);
  JR[2,2] = 0.5*n[0]*vR - 0.5*dsmaxR[1]*(qR[2] - qL[2]);
  JR[2,1] = n[1]*vR + 0.5*n[0]*uR - 0.5*smax - 0.5*dsmaxR[2]*(qR[2] - qL[2]);

}

#endif
