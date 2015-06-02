
#include "Calculator.h"

using namespace std;

//---------------------------------------------------//
// Constructor
//---------------------------------------------------//
Calculator::Calculator(TVector3 v_norm, TVector3 planeCenter,
		       float xwidth, float ywidth, float zwidth,
		       float ind0, float ind1, bool dbg):
  zRotation(NULL),
  identity(NULL),
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
  TMatrixF id = TMatrixF(3,3);
  id(0,0) = 1;
  id(1,1) = 1;
  id(2,2) = 1;
  identity = new TMatrixF(id);


  // Set the rotation matrix
  //setZRotation();

}

//---------------------------------------------------//
// Destructor
//---------------------------------------------------//
Calculator::~Calculator()
{

  //if( firstRotation ) firstRotation->Delete();
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

  if(m_dbg){
    ptp.Print();
    pap.Print();
    TVector3 temp = (*zRotation * norm);
    norm.Print();
    temp.Print();    
    temp = *zRotation * temp;
    temp.Print();
  }

  // Shift by the z position of the plane.
  // Here this isn't really setup to handle yet,
  // so assume plane is at 0 already
  // TODO: apply z-shift.

  // Now minimize to find the interaction point
  TVector3 pip = scan(ptp,pap);

  cout<<"---------------- Snell Check ----------------------"<<endl;
  cout<<n0 * sin(atan(sqrt(pow(ptp.X()-pip.X(),2)+pow(ptp.Y()-pip.Y(),2))/fabs(ptp.Z())))<<endl;
  cout<<n1 * sin(atan(sqrt(pow(pap.X()-pip.X(),2)+pow(pap.Y()-pip.Y(),2))/fabs(pap.Z())))<<endl;
  cout<<pip.X()<<" "<<pip.Y()<<" "<<pip.Z()<<endl;
  cout<<"Normal: "<<norm.X()<<" "<<norm.Y()<<" "<<norm.Z()<<endl;

  // Now need to rotate back
  setZRotation(norm, TVector3(0,0,1));
  (*zRotation * ptp).Print();
  (*zRotation * pap).Print();

  if(*zRotation == *identity) return pip;
  return (*zRotation) * pip;
  
}

//---------------------------------------------------//
// Scan x-y plane to find interaction point
//---------------------------------------------------//
TVector3 Calculator::scan(TVector3 pt, TVector3 pa)
{

  // Brute force scan
  int nxsteps = 1000;
  float xmin  = planeC.X() - x_ice/2.;
  float xmax  = planeC.X() + x_ice/2.;
  float xstep = (xmax-xmin)/nxsteps;

  int nysteps = 1000;
  float ymin  = planeC.Y() - y_ice/2.;
  float ymax  = planeC.Y() + y_ice/2.;
  float ystep = (ymax-ymin)/nysteps;

  // Loop over and minimize the gradient in each 
  // coordinate.
  TVector3 pi = TVector3(0,0,0); // int point
  float gradX = 9999;
  float gradY = 9999;
  TVector3 temp = TVector3(0,0,0);
  for(int ix = 0; ix<nxsteps; ++ix){
    float xi = xmin + ix * xstep;
    
    for(int iy = 0; iy<nysteps; ++iy){
      float yi = ymin + iy * ystep;
    
      temp.SetXYZ(xi,yi,0);
      float gx = fabs(getGradient(pt,pa,temp,0));
      float gy = fabs(getGradient(pt,pa,temp,1));

      //cout<<"xi: "<<xi<<" gx: "<<gx<<" current grad: "<<gradX<<endl;
      if( gx <= gradX && gy <= gradY ){
	gradX = gx;
	gradY = gy;
	pi.SetXYZ(xi,0,0);
      }
	

    }// end loop over ysteps

  }// end loop over xsteps

  return pi;

}

//---------------------------------------------------//
// Method to calculate the gradient
//---------------------------------------------------//
float Calculator::getGradient(TVector3 pt, TVector3 pa, TVector3 pi, 
	   int option)
{

  float grad = 9999;

  // Decide what coordinate we are looking at
  float ct = 0, ca = 0, ci = 0;
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

  // Get the constant factors
  if( (pt-pi).Mag() == 0 ) return grad;
  if( (pa-pi).Mag() == 0 ) return grad;

  float t_const = n0 / ( (pt-pi).Mag());
  float a_const = n1 / ( (pa-pi).Mag());

  // Determine the sign
  int t_sign = ci > ct ? -1 : 1;
  int a_sign = ci > ca ? -1 : 1;
  
  // Get gradient for this term
  grad = -t_const * t_sign * (ct - ci) -
    a_const * a_sign * (ca - ci);
  
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
  float s = v.Mag();

  // Get C
  float c = normal.Dot( k );

  // Setup the skew symmetrix matrix
  TMatrixF v_x = TMatrixF(3,3);
  v_x(0,1) = -v(2);
  v_x(0,2) = v(1);
  v_x(1,0) = v(2);
  v_x(1,2) = -v(0);
  v_x(2,0) = -v(1);
  v_x(2,1) = v(0);

  // Set the first rotation
  // Special case: n || k
  if( normal*k == 1 ) zRotation = new TMatrixF(*identity);  
  else                zRotation = new TMatrixF(*identity + v_x + v_x*v_x*((1-c)/(s*s)));    

  // Debug info
  if( m_dbg ){
    cout<<"-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-="<<endl;
    zRotation->Print();
    cout<<"-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-="<<endl;
  }

}


