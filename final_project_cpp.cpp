// C++ final project main file - PHYS30762

#include <iostream>
#include <string>
#include <vector>
#include "final_project_header.h"

int main()
{
  four_momentum up_quark_4_vec(4, 2, 2, 2);

  quark up_quark("up");
  up_quark.set_colour("red");
  up_quark.set_four_momenta(up_quark_4_vec);
  up_quark.set_mass();

  std::cout << "mass: " << up_quark.get_mass() << std::endl;

  return 0;
}