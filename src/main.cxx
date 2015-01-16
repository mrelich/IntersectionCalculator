
#include "Calculator.h"
#include <iostream>

int main()
{
  
  // Create calculator with normal in z direction
  Calculator* calc = new Calculator(0,0,1,1.75,1,true);

  // some vars
  float xi = 0, yi = 0, zi = 0;

  // Now have some testing  
  calc->getIntPoint(xi,yi,zi,
		    -1, 0, -1,
		    1, 0, 1);

  std::cout<<"Final: "<<xi<<","<<yi<<","<<zi<<std::endl;

  delete calc;
  
  return 0;
}
