
#include "Calculator.h"

using namespace std;

//---------------------------------------------------//
// Constructor
//---------------------------------------------------//
Calculator::Calculator(float nx, float ny, float nz,
		       float x0, float y0, float z0,
		       float ind0, float ind1, bool dbg):
  firstRotation(NULL),
  secondRotation(NULL),
  m_dbg(false)
{

  // Set variables
  norm  = TVector3(nx,ny,nz);
  plane = TVector3(x0,y0,z0);
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
  // Also apply to plane point
  TVector3 p0p    = *firstRotation * p0;
  TVector3 p1p    = *firstRotation * p1;
  TVector3 planep = *firstRotation * plane; 
  
  // Set the second rotation and apply 
  // (pp for double prime)
  setSecondRotation(p0p, p1p);
  TVector3 p0pp    = *secondRotation * p0p;
  TVector3 p1pp    = *secondRotation * p1p;
  TVector3 planepp = *secondRotation * planep;

  // Now translate the two points such that we
  // remove the z dependence.  By this point, 
  // the plane already lies in the x-z plane only
  TVector3 p0ppp = p0pp - planepp;
  TVector3 p1ppp = p1pp - planepp;

  // Now calculate beta
  float beta = (p0ppp.Mag2() + p1ppp.Mag2() - (p0ppp-p1ppp).Mag2())/(2*p0ppp.Mag()*p1ppp.Mag());


  // Solve for theta_i and theta_r
  float theta_r = TMath::ATan(-(TMath::Sin(TMath::ACos(beta)))/(beta + n1/n0));
  float theta_i = theta_r + TMath::ACos(beta) - TMath::Pi();
  
  // Now get xippp that satisfies the constraint
  float xippp = p1ppp.z() * tan(theta_r) + p1ppp.x();
  if( p0ppp.x() < p1ppp.x() && !(xippp+p0ppp.x() < xippp && xippp < xippp+p0ppp.x()))
    xippp = p1ppp.x() - p1ppp.z() * tan(theta_r);
  else if(p1ppp.x() < p0ppp.x() && !(xippp+p1ppp.x() < xippp && xippp < xippp+p0ppp.x()))
    xippp = p1ppp.x() - p1ppp.z() * tan(theta_r);
  if( p0ppp.x() < p1ppp.x() && !(xippp+p0ppp.x() < xippp && xippp < xippp+p1ppp.x())){
    cout<<"There appears to be no solution"<<endl;
    xippp = 0;
  }
  else if( p1ppp.x() < p0ppp.x() && !(xippp+p1ppp.x() < xippp && xippp < xippp+p0ppp.x())){
    cout<<"There appears to be no solution"<<endl;
    xippp = 0;
  }

  // Construct pippp (pi''')
  TVector3 pippp = TVector3(xippp,0,0);
  
  // Now translate and rotate back
  TVector3 pi = rotateBack(pippp+planepp);
  
  // Set variables
  xi = pi.x();
  yi = pi.y();
  zi = pi.z();

  // Print some debug info
  if( m_dbg ){
    cout<<"---------------------------------------"<<endl;
    p0.Print();
    p1.Print();
    plane.Print();
    cout<<"---------------------------------------"<<endl;
    p0p.Print();
    p1p.Print();
    planep.Print();
    cout<<"---------------------------------------"<<endl;
    p0pp.Print();
    p1pp.Print();
    planepp.Print();
    cout<<"---------------------------------------"<<endl;
    p0ppp.Print();
    p1ppp.Print();
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
  if( m_dbg ){
    cout<<"-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-="<<endl;
    firstRotation->Print();
    cout<<"-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-="<<endl;
  }

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
  if( m_dbg ){
    cout<<"-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-="<<endl;
    secondRotation->Print();
    cout<<"-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-="<<endl;
  }

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
