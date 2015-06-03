
#include "Calculator.h"

using namespace std;

//---------------------------------------------------//
// Constructor
//---------------------------------------------------//
Calculator::Calculator(TVector3 v_norm, TVector3 planeCenter,
		       double xwidth, double ywidth, double zwidth,
		       double ind0, double ind1, bool dbg):
  zRotation(NULL),
  identity(NULL),
  m_tolerance(1),
  m_dbg(false)
{

  // Set variables
  norm   = v_norm;
  planeC = planeCenter;
  x_ice  = xwidth;
  y_ice  = ywidth;
  z_ice  = zwidth;
  n0 = ind0;
  n1 = ind1;
  m_dbg = dbg;

  // Set identidy matrix
  TMatrixD id = TMatrixT<double>(3,3);
  id(0,0) = 1;
  id(1,1) = 1;
  id(2,2) = 1;
  identity = new TMatrixT<double>(id);


}

//---------------------------------------------------//
// Destructor
//---------------------------------------------------//
Calculator::~Calculator()
{

  if( zRotation ) zRotation->Delete();

}

//---------------------------------------------------//
// Get interaction point
//---------------------------------------------------//
TVector3 Calculator::getIntPoint(TVector3 pt, TVector3 pa)
{

  // Set the rotation matrix
  setZRotation(TVector3(0,0,1), norm);

  // Rotate the points pt and pa
  TVector3 ptp = *zRotation * pt;
  TVector3 pap = *zRotation * pa;

  // Shift by the z position of the plane.
  TVector3 zshift = TVector3(0,0,z_ice/2.);
  TVector3 ptpp = ptp - zshift;
  TVector3 papp = pap - zshift;


  // Now minimize to find the interaction point
  TVector3 pipp = scan(ptpp,papp);

  // Snell's law check
  float ang0 = atan(sqrt(pow(ptpp.X()-pipp.X(),2)+pow(ptpp.Y()-pipp.Y(),2))/fabs(ptpp.Z()));
  float ang1 = atan(sqrt(pow(papp.X()-pipp.X(),2)+pow(papp.Y()-pipp.Y(),2))/fabs(papp.Z()));
  float expectedAngle = (ang0 < asin(n1/n0)) ? asin(n0*sin(ang0)/n1) : 99999;
  
  if(m_dbg){
    cout<<"------------------------- Snell's Law Check --------------------------"<<endl;
    pipp.Print();
    cout<<"Angles: "<<ang0*180/TMath::Pi()<<" "<<ang1*180/TMath::Pi()<<endl;
    cout<<"Expected refracted angle: "<<expectedAngle*180/TMath::Pi()<<endl;
    cout<<n0 * sin(ang0)<<endl;
    cout<<n1 * sin(ang1)<<endl;

  }

  if( ang0 >= asin(n1/n0) || fabs(1 - expectedAngle/ang1) > m_tolerance ){
    cout<<"Fail Snell's law test: "<<ang0*180/TMath::Pi()
	<<" max allowd: "<<asin(n1/n0)*180/TMath::Pi()<<endl;
    cout<<"\tDifference with expected: "<<fabs(1-expectedAngle/ang1)<<endl;
    pipp.SetXYZ(-9999,-9999,-9999);
    return pipp;
  }
    
  // Now need to shift and rotate back
  setZRotation(norm, TVector3(0,0,1));
  pipp += zshift;
  if(*zRotation == *identity) return pipp;
  return (*zRotation) * pipp;
  
}

//---------------------------------------------------//
// Scan x-y plane to find interaction point
//---------------------------------------------------//
TVector3 Calculator::scan(TVector3 pt, TVector3 pa)
{

  // Brute force scan
  int nxsteps = 1000;
  double xmin  = planeC.X() - x_ice/2.;
  double xmax  = planeC.X() + x_ice/2.;
  double xstep = (xmax-xmin)/nxsteps;

  int nysteps = 1000;
  double ymin  = planeC.Y() - y_ice/2.;
  double ymax  = planeC.Y() + y_ice/2.;
  double ystep = (ymax-ymin)/nysteps;

  // Loop over and minimize the gradient in each 
  // coordinate.
  TVector3 pi = TVector3(0,0,0); // int point
  
  double gradX = 9999;
  //double gradY = 9999;
  TVector3 temp = TVector3(0,0,0);
  for(int ix = 0; ix<nxsteps; ++ix){
    double xi = xmin + ix * xstep;
    
    for(int iy = 0; iy<nysteps; ++iy){
      double yi = ymin + iy * ystep;
    
      temp.SetXYZ(xi,yi,0);
      double gx = fabs(getGradient(pt,pa,temp,0));
      double gy = fabs(getGradient(pt,pa,temp,1));

      //if( gx <= gradX && gy <= gradY ){
      if( fabs(gx+gy) <= gradX ){
	//gradX = gx;
	gradX = fabs(gx+gy);
	//gradY = gy;
	pi.SetXYZ(xi,yi,0);
      }
	

    }// end loop over ysteps

  }// end loop over xsteps



  /*
  temp.SetXYZ(-0.94,0,0);
  double gx = fabs(getGradient(pt,pa,temp,0));		  
  cout<<"Gradient: "<<gx<<endl;
  pi.SetXYZ(temp.X(),0,0);

  temp.SetXYZ(-0.10,0,0);
  gx = fabs(getGradient(pt,pa,temp,0));		  
  cout<<"Gradient: "<<gx<<endl;
  pi.SetXYZ(temp.X(),0,0);
  */
  return pi;

}

//---------------------------------------------------//
// Method to calculate the gradient
//---------------------------------------------------//
double Calculator::getGradient(TVector3 pt, TVector3 pa, TVector3 pi, 
	   int option)
{

  double grad = 9999;

  // Decide what coordinate we are looking at
  double ct = 0, ca = 0, ci = 0;
  if(option == 0){
    ct = pt.X();
    ca = pa.X();
    ci = pi.X();
  }
  else if( option == 1 ){
    ct = pt.Y();
    ca = pa.Y();
    ci = pi.Y();
  }
  else return grad;

  // The point should lie somewhere between 
  // the original points
  if( !((ct <= ci && ci <= ca) || (ca <= ci && ci <= ct)) ) return grad;

  // Get the constant factors
  if( (pt-pi).Mag() == 0 ) return grad;
  if( (pa-pi).Mag() == 0 ) return grad;

  double t_const = n0 / ( (pt-pi).Mag());
  double a_const = n1 / ( (pa-pi).Mag());

  // Get gradient for this term
  grad = t_const * fabs(ct-ci) - a_const * fabs(ci-ca);


  if(false && m_dbg){
    cout<<endl;
    cout<<"\tct: "<<ct<<" ca: "<<ca<<" ci: "<<ci<<endl;
    cout<<"\tpt: "<<pt.X()<<" "<<pt.Y()<<" "<<pt.Z()<<endl;
    cout<<"\tpa: "<<pa.X()<<" "<<pa.Y()<<" "<<pa.Z()<<endl;
    cout<<"\tpi: "<<pi.X()<<" "<<pi.Y()<<" "<<pi.Z()<<endl;
    cout<<"\tMagt: "<<(pt-pi).Mag()<<endl;
    cout<<"\tMaga: "<<(pa-pi).Mag()<<endl;
    cout<<endl;
    //cout<<"\ttrack sign:  "<<t_sign<<endl;
    cout<<"\ttrack const: "<<t_const<<endl;
    cout<<"\ttrack diff:  "<<(ct-ci)<<endl;
    //cout<<"\tant sign:    "<<a_sign<<endl;
    cout<<"\tant const:   "<<a_const<<endl;
    cout<<"\tant diff:    "<<(ci-ca)<<endl;
  }  

  return grad;


}


//---------------------------------------------------//
// Set first rotation matrix
// This matrix will rotate the normal vector into
// the z-direction
//---------------------------------------------------//
void Calculator::setZRotation(TVector3 k, TVector3 normal)
{

  // To do this we first calculate some pieces
  // R = I + [v]_x + [v]_x^2 * (1-C)/S^2
  // v = n x k
  // s = ||v|| 
  // c = n * k
  // [v]_x = skew symmetrix matrix


  // Set k
  //TVector3 k = TVector3(0,0,1);

  // Get V
  TVector3 v = normal.Cross( k );

  // Get S
  double s = v.Mag();

  // Get C
  double c = normal.Dot( k );

  // Setup the skew symmetrix matrix
  TMatrixD v_x = TMatrixT<double>(3,3);
  v_x(0,1) = -v(2);
  v_x(0,2) = v(1);
  v_x(1,0) = v(2);
  v_x(1,2) = -v(0);
  v_x(2,0) = -v(1);
  v_x(2,1) = v(0);

  // Set the first rotation
  // Special case: n || k
  if( normal*k == 1 ) zRotation = new TMatrixT<double>(*identity);  
  else                zRotation = new TMatrixT<double>(*identity + v_x + v_x*v_x*((1-c)/(s*s)));    

  // Debug info
  if( m_dbg ){
    cout<<"-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-="<<endl;
    zRotation->Print();
    cout<<"-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-="<<endl;
  }

}


