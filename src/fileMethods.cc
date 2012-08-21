#ifndef FILEMETHODS_CC
#define FILEMETHODS_CC
#include <string>
#include <iostream>
#include <sstream>
#include <fstream>
#include "latticeSize.cc"
#include "Globals.cc"
#include "Polynomial.cc"
using namespace std;



string concat(string string1, string string2) {
    ostringstream s1, s2;
    s1 << string1;
    s2 << string2;
    return (s1.str() + s2.str());
}

string create_filename(const int zs) {
    ostringstream x, y, z, b1, b2, b3, q;
    q << Q;
    x << xSize;
    y << ySize;
    z << zs;

    b1 << BConds2;
    b2 << BConds1;
    b3 << BConds3;

    string fname = "";
    fname = fname + q.str() + "-";
    if (xSize < 10)
        fname = fname + "0" + x.str() + "x";
    else
        fname = fname + x.str() + "x";

    if (ySize < 10)
        fname = fname + "0" + y.str() + "x";
    else
        fname = fname + y.str() + "x";
    if (zs < 10)
        fname = fname + "0" + z.str() + "-";
    else
        fname = fname + z.str() + "-";

    fname = fname + b1.str() + b2.str() + b3.str();
    return fname;

}

void logError(string err) {

    string fname = "polynomials/";
    string fname1 = create_filename(zSize);
    fname = concat(fname, fname1);
    ofstream file_op;
    file_op.open("Error.log", ios::out | ios::app);
    if (!file_op.is_open())
        cout << " Could not write to file Error.log\n";
    else {
        file_op << endl << fname << endl;
        file_op << "Error:\n";
        file_op << err << endl;
        cout << "\n Error with " << fname << endl;
        cout << err;
        cout << "Appended to Error.log" << endl;
    }
    file_op.close();
}

void print2File(Polynomial poly) {
    string fname = "polynomials/";
    string fname1 = create_filename(zSize);
    fname = concat(fname, fname1);
    ofstream file_op;
    file_op.open(fname.c_str(), ios::out | ios::trunc);
    if (!file_op.is_open())
        cout << " Could not write to file: " << fname;
    else {
        file_op << poly;
        cout << "\nSaved pol to " << fname << endl;
    }
    file_op.close();
}

void print2File(Polynomial poly, int zs) {
    string fname = "polynomials/";
    string fname1 = create_filename(zs);
    fname = concat(fname, fname1);

    ofstream file_op;
    file_op.open(fname.c_str(), ios::out | ios::trunc);
    if (!file_op.is_open())
        cout << " Could not write to file: " << fname;
    else {
        file_op << poly;
        cout << "\nSaved pol to " << fname << endl;
    }
    file_op.close();
}

void truncateFile(string filename) {
    ofstream file_op;
    file_op.open(filename.c_str(), ios::out | ios::trunc);
    if (!file_op.is_open())
        cout << " Could not truncate file: " << filename;
    file_op.close();
}

void vpol_to_file(string filename, const vector<Polynomial> &vpol) {
    ofstream file_op;
    file_op.open(filename.c_str(), ios::out | ios::trunc);
    if (!file_op.is_open())
        cout << " Could not write to file: " << filename;
    else {
        for (int i = 0; i < vpol.size(); i++) {
            file_op << -1 << endl;
            file_op << vpol[i];
        }
    }
    file_op.close();
}

void pol_to_vpol_file(string filename, const Polynomial &p) {
    ofstream file_op;
    file_op.open(filename.c_str(), ios::out | ios::app);
    if (!file_op.is_open())
        cout << " Could not write to file: " << filename;
    else {
        file_op << -1 << endl;
        file_op << p;
    }
}

void file_to_vpol(vector<Polynomial> &vpol, string filename) {
    ifstream file_op;
    file_op.open(filename.c_str());
    if (!file_op.is_open())
        cout << " Could not read file: " << filename;
    else {
        mpz_class coeff;
        int vpol_index = -1;
        int expo = 0;
        while (file_op >> coeff) {
            if (coeff == -1) {
                vpol_index++;
                vpol[vpol_index].clear();
                expo = 0;
            } else {
                vpol[vpol_index].push_back(coeff, expo);
                expo++;
            }
        }
    }
    file_op.close();
}

void file_to_vpol(vector<Polynomial> &vpol, const G_TABLE2 &ylist, string filename) {
    if (vpol.size() < ylist.size()) {
        Polynomial p;
        while (vpol.size() < ylist.size()) {
            vpol.push_back(p);
        }
    }
    ifstream file_op;
    file_op.open(filename.c_str());
    if (!file_op.is_open())
        cout << " Could not read file: " << filename;
    else {
        mpz_class coeff;
        int vpol_index = -1;
        int expo = 0;
        int ylist_index = -1;
        while (file_op >> coeff) {
            if (coeff == -1) {
                ylist_index++;
                vpol_index = ylist[ylist_index][ySize];
                vpol[vpol_index].clear();
                expo = 0;
            } else {
                vpol[vpol_index].push_back(coeff, expo);
                expo++;
            }
        }
    }
    file_op.close();
}
#endif
