#ifndef Storage_y_CC
#define Storage_y_CC

#include <set>
#include <vector>
#include <iostream>
#include <sstream>
#include <fstream>
#include <algorithm>
#include "latticeSize.cc"
#include "Storage.cc"
#include "Grid.cc"
using namespace std;

// {{{ Storage_y class

/**
 *   For storing vlists the size is 2^xsize
 */
class Storage_y {
public:

	Storage_y() :
			pol_no(-1), delta1(0) {
	}
	void push_back(const int index, int pow);
	void print() const;
	void clear();
	void partialClear();
	bool operator==(const Storage_y &s) const;

	int operator[](const int i) const {
		return list[i];
	}
	bool operator<(const Storage_y &s) const;
	bool equals(const Storage_y &s) const;  // compare everything
	bool y_equals(const Storage_y &s) const;  // compare everything
	bool x_equals(const Storage_y &s) const;  // compare everything
	Storage_y order_it() const;
	Storage_y order_it_asending() const;

	int getPol_no() const {
		return pol_no;
	}

	void setPol_no(int id) {
		pol_no = id;
	}
	Storage_y transform(const int &k) const;
	Storage_y transpose() const;
	Storage_y reflect() const;
	Storage_y shift() const;
	Storage_y reflect_x() const;
	Storage_y shift_x() const;
	void calc_subgraphs();

	int getDelta() const {
		return delta1;
	}
	int getDelta(const vector<int> &ylist_bc3,
	    const vector<Storage> &GrowthMatrix) const;

	vector<pair<int, int> > get_vpr() const {
		return vpr;
	}

	vector<pair<int, int> > get_inversevpr() const {
		return inverse_vpr;
	}
	void output() const;
	void setSmallest_x();
	void setSmallest();
	void y_less(const Storage_y &s);

private:
	int list[ySize];
	int pol_no;
	vector<pair<int, int> > vpr;
	vector<pair<int, int> > inverse_vpr;
	int delta1;
};
Grid gd_y(xSize, ySize);

int Storage_y::getDelta(const vector<int> &ylist_bc3,
    const vector<Storage> &GrowthMatrix) const {
	int delta = 0;
	for (int i = 0; i < ySize; i++) {
		delta += GrowthMatrix[ARRSIZE][list[i]];
		if (BConds3) delta += GrowthMatrix[ylist_bc3[i]][list[i]];
	}
	for (int i = 1; i < ySize; i++) {
		delta += GrowthMatrix[list[i - 1]][list[i]];
	}
	if (BConds1) delta += GrowthMatrix[list[0]][list[ySize - 1]];
	return delta;
}

void Storage_y::y_less(const Storage_y &s) {
	Storage_y temp = s;
	Storage_y tempR = s.reflect();
	for (int i = 0; i < ySize; i++) {
		if (temp < *this) {
			*this = temp;
		}
		if (tempR < *this) *this = tempR;
		temp = temp.shift();
		tempR = temp.reflect();
	}
}

void Storage_y::setSmallest_x() {
	Storage_y temp = *this;
	Storage_y tempR = this->reflect_x();
	vector<Storage_y> tempT(transformARR.size()), tempRT(transformARR.size());
	for (int i = 0; i < transformARR.size(); i++) {
		tempT[i] = this->transform(i);
		tempRT[i] = tempT[i].reflect_x();
	}
	//  Storage_y tempT = this->transform();
	//Storage_y tempRT = tempT.reflect_x();
	int i = 0;
	if (!BConds2) i = xSize - 1;
	for (; i < xSize; i++) {
		this->y_less(temp);
		this->y_less(tempR);
		for (int i = 0; i < transformARR.size(); i++) {
			this->y_less(tempT[i]);
			this->y_less(tempRT[i]);
		}
		temp = temp.shift_x();
		tempR = temp.reflect_x();
		for (int i = 0; i < transformARR.size(); i++) {
			tempT[i] = temp.transform(i);
			tempRT[i] = tempT[i].reflect_x();
		}
	}
}

void Storage_y::setSmallest() {
	Storage_y temp1 = *this;
	temp1.setSmallest_x();
	if (temp1 < *this) *this = temp1;
	if (xSize == ySize && BConds1 == BConds2) {
		Storage_y temp = this->transpose();
		temp.setSmallest_x();
		if (temp < *this) *this = temp;
	}
}

bool Storage_y::operator<(const Storage_y &s) const {
	for (int i = 0; i < ySize; i++) {
		if (list[i] < s[i]) return true;
		else if (list[i] > s[i]) return false;
	}
	return false;
}

void Storage_y::output() const {
	for (int i = 0; i < ySize; i++)
		cout << list[i] << " ";
	cout << "\t" << pol_no;
	cout << endl;
}

void Storage_y::calc_subgraphs() {
	for (int i = 0; i < ySize; i++) {
		for (int j = 0; j < xSize; j++) {
			gd_y.setGrid(j, i, transposeARR[list[i]][j]);
		}
	}
	this->delta1 = gd_y.calcEnergy_brute2d_1();
	//gd_y.calculate_edge_subgraphs(vpr, false);
	//this->delta1 = 0;
	//for(vector<pair<int, int> >::size_type i = 0; i < vpr.size(); i++){
	//  this->delta1 = this->delta1 + vpr[i].second;
	//}
	//gd_y.calculate_edge_subgraphs(inverse_vpr, true);
}
int tempARR[xSize][ySize];

Storage_y Storage_y::transpose() const {
	Storage_y temp;
	for (int ind = 0; ind < xSize; ind++) {
		for (int j = 0; j < ySize; j++) {
			int i = list[j];
			tempARR[ind][j] = transposeARR[i][ySize - 1 - ind];
		}
	}
	for (int k = 0; k < ySize; k++) {
		// search for it;
		bool found = false;
		for (int j = 0; j < ARRSIZE && !found; j++) {
			bool silu = true;
			for (int l = 0; l < ySize && silu; l++) {
				if (transposeARR[j][l] != tempARR[k][l]) {
					silu = false;
				}
			}
			if (silu) {
				found = true;
				temp.push_back(k, j);
			}
		}
	}
	return temp;
}

Storage_y Storage_y::transform(const int &k) const {
	Storage_y temp;
	for (int i = 0; i < ySize; i++) {
		int j = transformARR[k][this->list[i]];  //ARRSIZE - 1 - this->list[i];
		//int j = ARRSIZE - 1 - this->list[i];
		temp.push_back(i, j);
	}
	return temp;
}

Storage_y Storage_y::reflect() const {
	Storage_y temp = *this;
	for (int i = 0; i < ySize; i++)
		temp.push_back(i, list[ySize - 1 - i]);
	return temp;
}

Storage_y Storage_y::reflect_x() const {
	Storage_y temp = *this;
	for (int i = 0; i < ySize; i++)
		temp.push_back(i, reflectxARR[list[i]]);
	return temp;
}

Storage_y Storage_y::shift_x() const {
	Storage_y temp = *this;
	for (int i = 0; i < ySize; i++) {
		temp.push_back(i, shiftxARR[list[i]]);
	}
	return temp;
}

Storage_y Storage_y::shift() const {
	if (BConds1) {
		Storage_y temp = *this;
		for (int i = 0; i < ySize - 1; i++) {
			temp.push_back(i, list[i + 1]);
		}
		temp.push_back(ySize - 1, list[0]);
		return temp;
	}
	return *this;
}

Storage_y Storage_y::order_it_asending() const {
	Storage_y temp;
	multiset<int> s(list, list + ySize);
	multiset<int>::iterator sittr, eittr;
	sittr = s.begin();
	eittr = s.end();
	int i = 0;
	for (; sittr != eittr; sittr++) {
		temp.push_back(i, *sittr);
		i++;
	}
	return temp;
}

Storage_y Storage_y::order_it() const {
	Storage_y temp;
	multiset<int> s(list, list + ySize);
	multiset<int>::iterator sittr, eittr;
	sittr = s.begin();
	eittr = s.end();
	int i = ySize;
	for (; sittr != eittr; sittr++) {
		i--;
		temp.push_back(i, *sittr);
	}
	return temp;
}

bool Storage_y::y_equals(const Storage_y &s) const {
	Storage_y temp = s;
	Storage_y tempR = s.reflect();
	for (int i = 0; i < ySize; i++) {
		if (this->equals(temp)) return true;
		if (this->equals(tempR)) return true;
		temp = temp.shift();
		tempR = temp.reflect();
	}
	return false;
}

bool Storage_y::x_equals(const Storage_y &s) const {
	Storage_y thisOrdered = this->order_it();
	Storage_y temp = s;
	Storage_y tempR = s.reflect_x();
	vector<Storage_y> tempT(transformARR.size()), tempRT(transformARR.size());
	for (int i = 0; i < transformARR.size(); i++) {
		tempT[i] = s.transform(i);
		tempRT[i] = tempT[i].reflect_x();
	}
	//  Storage_y tempT = s.transform();
	//  Storage_y tempRT = tempT.reflect_x();
	int j = 0;
	if (!BConds2) j = xSize - 1;
	for (; j < xSize; j++) {
		if (thisOrdered.equals(temp.order_it())) if (this->y_equals(temp)) return true;
		if (!tempR.equals(temp)) if (thisOrdered.equals(tempR.order_it())) if (this->y_equals(
		    tempR)) return true;

		for (int i = 0; i < transformARR.size(); i++) {
			if (thisOrdered.equals(tempT[i].order_it())) if (this->y_equals(tempT[i])) return true;
			if (thisOrdered.equals(tempRT[i].order_it())) if (this->y_equals(
			    tempRT[i])) return true;
		}
		temp = temp.shift_x();
		tempR = temp.reflect_x();
		for (int i = 0; i < transformARR.size(); i++) {
			tempT[i] = temp.transform(i);
			tempRT[i] = tempT[i].reflect_x();
		}
	}
	return false;
}

bool Storage_y::operator==(const Storage_y &s) const {
	if (delta1 == s.delta1) if (this->vpr == s.vpr
	    && this->inverse_vpr == s.inverse_vpr) {
		if (this->x_equals(s)) return true;
		else if (xSize == ySize && BConds1 == BConds2) {
			Storage_y temp = s.transpose();
			if (!s.equals(temp)) {
				return this->x_equals(temp);
			}
		}
	}
	return false;
}

bool Storage_y::equals(const Storage_y &s) const {
	for (int i = 0; i < ySize; i++)
		if (list[i] != s[i]) return false;
	return true;
}

void Storage_y::clear() {
	for (int i = 0; i < ySize; i++)
		list[i] = -1;
}

/**
 *  Don't clear the group_id
 */
void Storage_y::partialClear() {
	for (int i = 0; i < ySize; i++)
		list[i] = -1;
}

void Storage_y::push_back(const int index, int pow) {
	this->list[index] = pow;
}

void Storage_y::print() const {
	for (int i = 0; i < ySize; i++)
		cout << list[i] << " ";
	cout << endl;
}
ostream & operator<<(ostream & out, const Storage_y &s) {
    for (int i = 0; i < ySize; i++)
        out << s[i] << " ";
    out << " delta: " << s.getDelta() << endl;
    Storage_y temp = s;
    temp.calc_subgraphs();
    out << "vpr: " << temp.get_vpr() << endl;
    out << "ivpr: " << temp.get_inversevpr() << endl;
    return out;
}

// }}}

#endif

