#ifndef Globals_CC
#define Globals_CC

/*
 * Globals.cc
 *
 *  Created on: 19 Aug 2012
 *      Author: yogi
 */

#include "latticeSize.cc"
#include "Globals.h"

time_t startTime;
const int MAXPOWER = xSize + 1;
int shiftxARR[ARRSIZE];
int reflectxARR[ARRSIZE];
int transposeARR[ARRSIZE][xSize];
int subARR[xSize][ySize]; // For the subclasses
std::vector<std::vector<int> > transformARR;


typedef vector<int> G_TABLE1;
typedef vector<G_TABLE1> G_TABLE2;
typedef vector<G_TABLE2> G_TABLE3;
G_TABLE2 GROUPTABLE;

#ifndef MAXPOLYSIZE
#define MAXPOLYSIZE (3 * xSize * ySize * zSize) + 1// Max power of lattice
#endif
string VPOL_FILENAME;

ostream & operator<<(ostream & out, const vector<pair<int, int> > &s) {
    for (vector<pair<int, int> >::size_type i = 0; i < s.size(); i++)
        out << "(" << s[i].first << ", " << s[i].second << ")\t";
    return out;
}

ostream & operator<<(ostream & out, const vector<int> &allCombos) {
    for (unsigned int i = 0; i < allCombos.size(); i++) {
        out << allCombos[i] << endl;
    }
    out << endl;
    return out;
}

#endif

