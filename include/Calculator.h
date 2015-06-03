#ifndef Calculator_h
#define Calculator_h

//-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-//
// This tool can be used to calculate the interaction point on a //
// plane when refraction is involved. The user must input:       //
//   * Indices of refraction of materials                        //
//   * The two points                                            //
//   * Normal vector to plane                                    //
//-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-//

#include "TVector3.h"
#include "TMatrix.h"
#include "TMatrixT.h"
#include "TMath.h"
#include <iostream>

typedef TMatrixT<double> TMatrixD; // should already be included...

class Calculator
{

 public:
  
  // Constructor
  Calculator(TVector3 v_norm,      // Normal to the plane
	     TVector3 planeCenter, // Center of the plane position
	     double xwidth,         // distance of cube in x
	     double ywidth,         // distance of cube in y
	     double zwidth,         // distance of cube in z
	     double ind0,           // index of refraction 0
	     double ind1,           // index of refraction 1
	     bool dbg = false);        

  // Destructor
  virtual ~Calculator();


  // Method to retrieve intersection point
  TVector3 getIntPoint(TVector3 pt,    // Track position
		       TVector3 pa);   // Antenna position
  

  // Set the debug flag
  void setDebug(){ m_dbg = true; };

 private:

  // Calculate the rotation matrix for z direction
  void setZRotation(TVector3 k, TVector3 normal);

  // Scan values for int point
  TVector3 scan(TVector3 pt, TVector3 pa);

  // Get the gradient for a point
  double getGradient(TVector3 pt,    // track point
		    TVector3 pa,    // antenna point
		    TVector3 pi,    // potential interaction point
		    int option);

  // Variables
  TVector3 norm;            // Normal vector to plane
  TVector3 planeC;          // Center of the plane
  double x_ice;              // x width of iceblock
  double y_ice;              // y width of iceblock
  double z_ice;              // z width of iceblock

  double n0;                 // Index of refraction for material 0
  double n1;                 // Index of refraction for material 1

  TMatrixD* zRotation;      // Rotate into z direction
  TMatrixD* identity;       // Identity matix

  bool m_dbg;               // Debug flag

};

#endif
