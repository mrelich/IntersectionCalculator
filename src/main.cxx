
#include "Calculator.h"
#include <iostream>

int main()
{

  // Some basic information
  float n_ice = 1.78;
  float n_air = 1.0;
  TVector3 p_track = TVector3(-1,0,-1);
  TVector3 p_antenna = TVector3(1,0,1);

  // Ice properties
  float tilt = 30 * TMath::Pi() / 180.;
  //float tilt = 0 * TMath::Pi() / 180.;
  TVector3 normal = TVector3(sin(tilt),0,cos(tilt));
  TVector3 center = TVector3(0,0,0);
  float x_dim     = 2.0; // [m]
  float y_dim     = 1.0; // [m]
  float z_dim     = 0.3; // [m]
  
  // Create calculator with normal in z direction
  Calculator* calc = new Calculator(normal,
				    center,
				    x_dim,
				    y_dim,
				    z_dim,
				    n_ice,
				    n_air,
				    true);


  // Now have some testing  
  TVector3 p_int = calc->getIntPoint(p_track,
				     p_antenna);

  std::cout<<"Final: "<<p_int.X()<<" "<<p_int.Y()<<" "<<p_int.Z()<<std::endl;

  delete calc;
  
  return 0;
}
