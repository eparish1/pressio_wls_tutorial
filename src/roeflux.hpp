
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
void roeflux_kernel( flux_t  & F, const state_t & UL,
		     const state_t & UR, const normal_t n, double g)
{
// PURPOSE: This function calculates the flux for the Euler equations
// using the Roe flux function
// INPUTS:
//    UL: conservative state vector in left cell
//   UR: conservative state vector in right cell
//    n: normal pointing from the left cell to the right cell

// OUTPUTS:
//  F   : the flux out of the left cell (into the right cell)
  double es = 1.e-30;
  double hL = UL[0];
  double uL = UL[1]/(hL + es);
  double vL = UL[2]/(hL + es);
  double unL = uL*n[0] + vL*n[1];

  double pL = 0.5*g*std::pow(hL,2.);
  double FL[ 3 ];
  FL[0] = hL*unL;
  FL[1] = UL[1]*unL + pL*n[0];
  FL[2] = UL[2]*unL + pL*n[1];

  double hR = UR[0];
  double uR = UR[1]/(hR + es);
  double vR = UR[2]/(hR + es);
  double unR = uR*n[0] + vR*n[1];
  double pR = 0.5*g*std::pow(hR,2.);
  double FR[3];
  FR[0] = hR*unR;
  FR[1] = UR[1]*unR + pR*n[0];
  FR[2] = UR[2]*unR + pR*n[1];

  // rho average
  double hm = 0.5*(hL + hR);
  double um = (unL*std::pow(hL,0.5) + unR*std::pow(hR,0.5) )/( std::pow(hL,0.5) + std::pow(hR,0.5)  + es);
  // eigenvalues
  double smax = abs(um) + abs(std::pow(g*hm, 0.5));
  // flux assembly
  F[0]    = 0.5*(FL[0]+FR[0])-0.5*smax*(UR[0] - UL[0]);
  F[1]    = 0.5*(FL[1]+FR[1])-0.5*smax*(UR[1] - UL[1]);
  F[2]    = 0.5*(FL[2]+FR[2])-0.5*smax*(UR[2] - UL[2]);

}


void roeflux_jacobian(){

  JL[0,0] = -0.5*dsmaxL[0]*(qR[0] - qL[0]) + 0.5*(n[0]*qL[1] + n[1]*qL[2])  + 0.5*smax;
  JL[0,1] = 0.5*n[0]*qL[0] - 0.5*dsmaxL[1]*(qR[0] - qL[0]);
  JL[0,2] = 0.5*n[1]*qL[0] - 0.5*dsmaxL[2]*(qR[0] - qL[0]);


}

#endif
