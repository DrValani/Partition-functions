#ifndef Storage_CC
#define Storage_CC
#include <vector>
#include <set>
#include <iostream>
#include <sstream>
#include <fstream>
#include "latticeSize.cc"
#include "Globals.cc"
#include "Polynomial.cc"
#include "Source.h"
using namespace std;
// {{{ Storage class

//class Storage;
//ostream & operator<<(ostream & out, const Storage &s);

/**
 *   For storing vlists the size is 2^xsize
 */

void Storage::createGroup(const Storage &s1, const Storage &s2) {
	for (int i = 0; i < ARRSIZE; i++) {
		int j = ARRSIZE - 1;
		while (j >= 0) {
			if (s1[j] == s2[i]) {
				list[i] = j;
				break;
			}
			j--;
		}
	}
}


template<class T>
T Storage::applyGroup(const T &s1) const {
	T temp;
	for (int i = 0; i < ARRSIZE; i++) {
		int pow = s1[this->list[i]];
		temp.push_back(i, pow);
	}
	return temp;
}


bool Storage::equals(const Storage &s) const {
	for (int i = 0; i < ARRSIZE; i++)
		if (list[i] != s[i]) return false;
	return true;
}

void Storage::clear() {
	for (int i = 0; i < ARRSIZE; i++)
		list[i] = -1;
}

void Storage::push_back(const int index, int pow) {
	this->list[index] = pow;
}

int Storage::searchGroup(vector<Storage> &allg1) const {

	unsigned int i = 0;
	while (i < allg1.size()) {
		if (allg1[i].equals(*this)) {
			return i;
		}
		i++;
	}
	allg1.push_back(*this);
	return allg1.size() - 1;
}

// }}}

ostream & operator<<(ostream & out, const Storage &s) {
    for (int i = 0; i < ARRSIZE; i++) {
        out << s[i] << " ";
    }
    out << " -- groupid =" << s.getId();
    return out;
}

vector<Storage> TRANSFORMGROUP;


#endif

