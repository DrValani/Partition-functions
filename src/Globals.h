#ifndef GLOBALS_H
#define GLOBALS_H
#include <vector>
#include <time.h>
#include <string>
#include <iostream>
#include <sstream>
#include <fstream>
#include "latticeSize.cc"
using namespace std;

extern time_t startTime;
extern const int MAXPOWER;
extern int shiftxARR[];
extern int reflectxARR[];
extern int transposeARR[][xSize];
extern int subARR[][ySize]; // For the subclasses
extern std::vector<std::vector<int> > transformARR;


typedef vector<int> G_TABLE1;
typedef vector<G_TABLE1> G_TABLE2;
typedef vector<G_TABLE2> G_TABLE3;
extern G_TABLE2 GROUPTABLE;

//#ifndef MAXPOLYSIZE
//#define MAXPOLYSIZE (3 * xSize * ySize * zSize) + 1// Max power of lattice
//#endif

extern const int MAX_VPOL_SIZE = 80000;
extern string VPOL_FILENAME;

extern ostream & operator<<(ostream & out, const vector<pair<int, int> > &s);

extern ostream & operator<<(ostream & out, const vector<int> &allCombos);

#endif
