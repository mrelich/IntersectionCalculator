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
#include "TMath.h"
#include <iostream>

class Calculator
{

 public:
  
  // Constructor
  Calculator(float nx, float ny, float nz, // Normal vector
	     float ind0,                   // Index of refraction 0
	     float ind1,                   // Index of refraction 1
	     bool dbg = false);

  // Destructor
  virtual ~Calculator();

  // Method to retrieve intersection point
  void getIntPoint(float &xi, float &yi, float &zi, // Intersection Point
		   float x0, float y0, float z0,    // Point inside material 0
		   float x1, float y1, float z1);   // Point inside material 1
  

  // Set the debug flag
  void setDebug(){ m_dbg = true; };
 private:

  // Calculate the rotation matrix for z direction
  void setFirstRotation();

  // Calculate the rotation matrix for y direction
  void setSecondRotation(TVector3 p0, TVector3 p1);

  // Rotate the intersection point back to normal coords
  TVector3 rotateBack(TVector3 intPoint);

  // Variables
  TVector3 norm;          // Normal vector to plane

  float n0;               // Index of refraction for material 0
  float n1;               // Index of refraction for material 1

  TMatrixF* firstRotation;  // Rotate into z direction
  TMatrixF* secondRotation; // Rotate into x-y plane

  bool m_dbg;               // Debug flag

};

#endif
