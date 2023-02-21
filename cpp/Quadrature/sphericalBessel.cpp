#include <cmath>
#include <iostream>
#include <boost/math/special_functions.hpp>

#include "sphericalBessel.h"

double SphericalBessel(unsigned order, double x) { return boost::math::sph_bessel(order, x); }

// Using boost because apple clang doesn't have these yet even though part of c++17
