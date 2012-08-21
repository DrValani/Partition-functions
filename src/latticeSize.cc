#ifndef LATTICESIZE_CC
#define LATTICESIZE_CC
/**
 * Declare the size of the lattice here
 */
#include <math.h>
#include "latticeSize.cc"

const int Q = 2;
const int xSize = 3;
const int ySize = 3;
const int zSize = 3;
const int ARRSIZE = 8;// This value is Q^xSize fill in by hand
const bool BConds1 = 1;
const bool BConds2 = 1;
const bool BConds3 = 0;
const int PRINT_PHASE = 5;


#endif
