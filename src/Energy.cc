//============================================================================
// Name        : 1.cpp
// Author      : Yogendra Valani
// Version     :
// Copyright   : Your copyright notice
// Description : Hello World in C++, Ansi-style
//============================================================================

#include <vector>
#include <time.h>
#include <cstdlib>
#include "latticeSize.cc"
#include "Globals.cc"
#include "myTimer.cc"

template <class T, class S>
int search_list(T &elist, S &s);

// }}}
#include "Grid.cc"
#include "Polynomial.h"
#include "refPolynomial.cc"
#include "Storage.cc"
#include "Storage_y.cc"


#include "searchMethods.cc"
#include "fileMethods.cc"
#include "commonMethods.cc"
#include "buildGroups.cc"
#include "build2d_9.cc"
#include "build2d_10.cc"
#include "AddLayers2.cc"
#include "AddLayers.cc"

// }}}

// {{{ Main method

int main(int argc, char** argv) {
    time(&startTime);
    VPOL_FILENAME = create_filename(zSize);
    string vpol_filename = "/tmp";
    VPOL_FILENAME = concat(vpol_filename, VPOL_FILENAME);
    vector<vector<Storage> >  vlist1(ySize + 1);
    vector<vector<Storage> >  vpol_vlist(ySize + 1);


    makeGroups();
    G_TABLE1 transformGroups;
    for (int i = 0; i < TRANSFORMGROUP.size(); i++) {
        transformGroups.push_back(i);
    }

    G_TABLE2 ylistTable;
    vector<vector<Polynomial> > vpol(ySize + 1);
    build_2dLattice(vlist1, transformGroups, ylistTable, vpol[1], vpol_vlist[0]); //Make the 2d lattice. (8)
    cout << "built 2d lattice\n";

    cout << "create total ref pol\n";
    time(&startTime);
    refPolynomial total;
    if (true) {
        vector<refPolynomial> vrefpol;
        int vlist_index = vlist1[0].size() - 1;
        if (ySize == 1)
            vlist_index = 0; //vpol_vlist[0].size()-1;
        make_refpol(0,
                transformGroups,
                vrefpol,
                vlist1[0],
                vpol_vlist[0],
                vlist_index
                );
        total = vrefpol[0];
    }

    time(&startTime);
    if(BConds3){
        total.ylist_multiplicity(ylistTable);
        buildAdd(ylistTable, transformGroups, vlist1, vpol, vpol_vlist);
    }
    else{
        buildAdd(ylistTable, transformGroups, vlist1, vpol, vpol_vlist, total);
    }
    return 0;
}

// }}}
