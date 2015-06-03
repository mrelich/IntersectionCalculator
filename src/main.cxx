
#include "Calculator.h"
#include <iostream>

int main()
{

  // Some basic information
  double n_ice = 1.78;
  double n_air = 1.0;
  TVector3 p_track = TVector3(-1,0,-1);
  TVector3 p_antenna = TVector3(1,0,1);

  // Ice properties
  double tilt = 30 * TMath::Pi() / 180.;
  //double tilt = 0 * TMath::Pi() / 180.;
  TVector3 normal = TVector3(sin(tilt),0,cos(tilt));
  TVector3 center = TVector3(0,0,0);
  double x_dim     = 1.0; // [m]
  double y_dim     = 0.3; // [m]
  double z_dim     = 0.3; // [m]
  
  // Create calculator with normal in z direction
  Calculator* calc = new Calculator(normal,
				    center,
				    x_dim,
				    y_dim,
				    z_dim,
				    n_ice,
				    n_air,
				    false); //true);


  // Now have some testing  
  TVector3 p_int = calc->getIntPoint(p_track,
				     p_antenna);

  std::cout<<"Final: "<<p_int.X()<<" "<<p_int.Y()<<" "<<p_int.Z()<<std::endl;

  delete calc;
  
  return 0;
}
