
#include "Calculator.h"

using namespace std;

//---------------------------------------------------//
// Constructor
//---------------------------------------------------//
Calculator::Calculator(float nx, float ny, float nz,
		       float ind0, float ind1, bool dbg):
  firstRotation(NULL),
  secondRotation(NULL),
  m_dbg(false)
{

  // Set variables
  norm = TVector3(nx,ny,nz);
  n0 = ind0;
  n1 = ind1;
  m_dbg = dbg;

  // Set the rotation matrix
  setFirstRotation();

}

//---------------------------------------------------//
// Destructor
//---------------------------------------------------//
Calculator::~Calculator()
{

  if( firstRotation ) firstRotation->Delete();
  if( secondRotation ) secondRotation->Delete();

}

//---------------------------------------------------//
// Get interaction point
//---------------------------------------------------//
void Calculator::getIntPoint(float &xi, float &yi, float &zi,
			     float x0, float y0, float z0,
			     float x1, float y1, float z1)
{

  // Setup points
  TVector3 p0 = TVector3(x0,y0,z0);
  TVector3 p1 = TVector3(x1,y1,z1);
  
  // Apply first rotation (p for prime)
  TVector3 p0p = *firstRotation * p0;
  TVector3 p1p = *firstRotation * p1;
  
  // Set the second rotation and apply 
  // (pp for double prime)
  setSecondRotation(p0p, p1p);
  TVector3 p0pp = *secondRotation * p0p;
  TVector3 p1pp = *secondRotation * p1p;

  // Now calculate beta
  float beta = (p0pp.Mag2() + p1pp.Mag2() - (p0pp-p1pp).Mag())/(2*p0pp.Mag()*p1pp.Mag());
  
  // Solve for theta_i and theta_r
  float theta_r = -(TMath::Sin(TMath::ACos(beta)))/(beta + n1/n0);
  float theta_i = TMath::Pi() + theta_r + TMath::ACos(beta);
  
  // Now get zipp
  float zipp = (1/(1-TMath::Tan(theta_i)/TMath::Tan(theta_r))) *
    (p0pp.z()+p1pp.x()-p1pp.x()*TMath::Tan(theta_i)-p1pp.z()*
     TMath::Tan(theta_i)/TMath::Tan(theta_r));
  float xipp = p1pp.x() - (p1pp.z() - zipp)/TMath::Tan(theta_r);

  // Construct pipp (pi'')
  TVector3 pipp = TVector3(xipp,0,zipp);
  
  // Now rotate back
  TVector3 pi = rotateBack(pipp);
  
  // Set variables
  xi = pi.x();
  yi = pi.y();
  zi = pi.z();

  // Print some debug info
  if( m_dbg ){
    cout<<"---------------------------------------"<<endl;
    p0.Print();
    p1.Print();
    cout<<"---------------------------------------"<<endl;
    p0p.Print();
    p1p.Print();
    cout<<"---------------------------------------"<<endl;
    p0pp.Print();
    p1pp.Print();
    cout<<"---------------------------------------"<<endl;
    cout<<"Beta: "<<beta<<endl;
    cout<<"Ref:  "<<theta_r<<endl;
    cout<<"Inc:  "<<theta_i<<endl;
    cout<<"---------------------------------------"<<endl;
    pi.Print();
  }

}

//---------------------------------------------------//
// Set first rotation matrix
// This matrix will rotate the normal vector into
// just the z-direction
//---------------------------------------------------//
void Calculator::setFirstRotation()
{

  // To do this we first calculate some pieces
  // R = I + [v]_x + [v]_x^2 * (1-C)/S^2
  // v = n x k
  // s = ||v|| 
  // c = n * k
  // [v]_x = skew symmetrix matrix


  // Set k
  TVector3 k = TVector3(0,0,1);

  // Get V
  TVector3 v = norm.Cross( k );

  // Get S
  float s = v.Mag();

  // Get C
  float c = norm.Dot( k );

  // Set identidy matrix
  TMatrixF identity = TMatrixF(3,3);
  identity(0,0) = 1;
  identity(1,1) = 1;
  identity(2,2) = 1;
  
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
  if( norm*k == 1 ) firstRotation = new TMatrixF(identity);  
  else              firstRotation = new TMatrixF(identity - v_x + v_x*v_x*((1-c)/(s*s)));    

  // Debug info
  cout<<"Printing first"<<endl;
  if( m_dbg ) firstRotation->Print();
  cout<<"After first"<<endl;
}

//---------------------------------------------------//
// Set Second Rotation
// This is meant to rotate the x and y positions
// such that the two points lie in the x-y plane
//---------------------------------------------------//
void Calculator::setSecondRotation(TVector3 p0,
				   TVector3 p1)
{

  // Condition is that phi = atan(-y0/x0)
  // And then the rotation matrix is normal
  // 2D matrix with angle phi.
  
  // Get angle
  float phi = TMath::ATan(-p0.y()/p0.x());
  
  // Set matrix
  TMatrixF temp = TMatrixF(3,3);
  temp(0,0) = TMath::Cos(phi);
  temp(0,1) = -TMath::Sin(phi);
  temp(1,0) = TMath::Sin(phi);
  temp(1,1) = TMath::Cos(phi);
  temp(2,2) = 1;

  // set second matrix
  if( secondRotation ) secondRotation->Delete();
  secondRotation = new TMatrixF(temp);

  // Debug info
  cout<<"Printing second"<<endl;
  if( m_dbg ) secondRotation->Print();
  cout<<"After second"<<endl;
}

//---------------------------------------------------//
// Rotate vector back
//---------------------------------------------------//
TVector3 Calculator::rotateBack(TVector3 intPoint)
{

  // Rotate by negative second rotation
  TVector3 origPoint = (*secondRotation) * intPoint * (-1);

  // Rotate by negative of first
  origPoint = (*firstRotation) * origPoint*(-1);

  // Now return
  return origPoint;

}
