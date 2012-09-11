#ifndef refPolynomial_CC
#define refPolynomial_CC
#include <map>
#include <vector>
#include <algorithm>
#include "Polynomial.cc"
#include "Storage.cc"
using namespace std;
// {{{ refPolynomial Class

/**
 *  A reference to the polynomial class.
 *  max power is the highest power that we increase the polynomial by
 *  expo = ((polynomial position in vpol) * maxpower) + power to increase by
 *  eg. if raising vpol[15] to the power of 5. with maxpower at 10. then
 *  expo = 15 * 10 + 5 = 155;
 *  coeff is the number of times this has to be done.
 */
class refPolynomial {
public:

	refPolynomial() {
	}

	void clear() {
		s.clear();
	}
	void push_back(int expo, int coeff);
	ostream &output(ostream &out) const;
	bool operator==(const refPolynomial &p) const;
	void sum(Polynomial &p, const vector<Polynomial> &vPol, const int &maxpower,
	    Polynomial &p1) const;
	void sum(Polynomial &p, const vector<Polynomial> &vPol, const int &maxpower,
	    const int &Hamiltonian, Polynomial &p1) const;
	void sum(Polynomial &p, const vector<Polynomial> &vPol) const;
	void sum(refPolynomial &p, vector<refPolynomial> &vrefPol);
	void sum(const refPolynomial &p, vector<refPolynomial> &vrefPol,
	    const int &maxpower1, const int &maxpower2);
	void raise(const refPolynomial &p, const int &power, const int &maxpower,
	    const int &maxpower1, const int &coeff);

	map<int, int>::iterator begin() {
		return s.begin();
	}

	map<int, int>::iterator end() {
		return s.end();
	}

	map<int, int>::const_iterator begin() const {
		return s.begin();
	}

	map<int, int>::const_iterator end() const {
		return s.end();
	}

	void next(map<int, int>::iterator &st) {
		++st;
	}

	void next(map<int, int>::const_iterator &st) const {
		++st;
	}

	void ylist_multiplicity(vector<vector<int> > &ylistTable) const;

private:
	map<int, int> s;

	void buildvectorRefs(vector<int>& counter,
			vector<vector<int> >& vref) const;
};


void refPolynomial::sum(refPolynomial &p, vector<refPolynomial> &s2list) {
	map<int, int>::iterator st = p.begin(), en = p.end(), s2st, s2en;
	while (st != en) {
		s2st = s2list[(*st).first].begin();
		s2en = s2list[(*st).first].end();
		while (s2st != s2en) {
			this->push_back((*s2st).first, (*s2st).second * (*st).second);
			s2list[(*st).first].next(s2st);
		}
		//*this += s2list[(*st).first]  * (*st).second;
		p.next(st);
	}
}

void refPolynomial::sum(Polynomial &p, const vector<Polynomial> &vPol,
    const int &maxpower, Polynomial &p1) const {

	  int Hamiltonian = 0;
	  sum(p, vPol, maxpower, Hamiltonian, p1);
}

void refPolynomial::buildvectorRefs(vector<int>& counter,
		vector<vector<int> >& vref) const {

	vector<int> temp(1, 0);
	map<int, int>::const_iterator st = s.begin(), en = s.end();
	while (st != en) {
		bool found = false;
		for (vector<int>::size_type j = 0; j < counter.size() && !found; j++) {
			if ((*st).second == counter[j]) {
				found = true;
				vref[j].push_back((*st).first);
			}
		}
		if (!found) {
			counter.push_back((*st).second);
			temp[0] = (*st).first;
			vref.push_back(temp);
		}
		++st;
	}
}

void refPolynomial::sum(Polynomial &p, const vector<Polynomial> &vPol,
    const int &maxpower, const int &Hamiltonian, Polynomial &p1) const {
	vector<int> counter;
	vector < vector<int> > vref;
	buildvectorRefs(counter, vref);

	for (vector<int>::size_type i = 0; i < counter.size(); i++) {
		p1.clear();
		for (vector<int>::size_type k = 0; k < vref[i].size(); k++) {
			int place = vref[i][k] / maxpower;
			int power = (vref[i][k] % maxpower) + Hamiltonian;
			p1.plusPower(vPol[place], power);
		}
		p.multiPlus(p1, counter[i]);
	}
}

void refPolynomial::sum(Polynomial &p, const vector<Polynomial> &vPol) const {
	vector<int> counter;
	vector < vector<int> > vref;
	buildvectorRefs(counter, vref);

	Polynomial p1;
	for (vector<int>::size_type i = 0; i < counter.size(); i++) {
		p1.clear();
		for (vector<int>::size_type k = 0; k < vref[i].size(); k++)
			p1 += vPol[vref[i][k]];
		p.multiPlus(p1, counter[i]);
	}
}

void refPolynomial::raise(const refPolynomial &refpol2, const int &power,
    const int &maxpower, const int &maxpower1, const int &coeff) {
	map<int, int>::const_iterator st = refpol2.begin(), en = refpol2.end(), s2st,
	    s2en;
	while (st != en) {
		int place = (*st).first / maxpower1;
		int pow = ((*st).first % maxpower1) + power;
		int expo = (place * maxpower) + pow;
		this->push_back(expo, (*st).second * coeff);
		++st;
	}
}

void refPolynomial::sum(const refPolynomial &refpol2, vector<refPolynomial> &s2list,
    const int &maxpower1, const int &maxpower2) {
	const int maxpower = maxpower1 + maxpower2;
	map<int, int>::const_iterator st = refpol2.begin(), en = refpol2.end();
	while (st != en) {
		int place = (*st).first / maxpower2;
		int power = (*st).first % maxpower2;
		this->raise(s2list[place], power, maxpower, maxpower1, (*st).second);
		++st;
	}

}

void refPolynomial::ylist_multiplicity(
    vector<vector<int> > & ylistTable) const {
	for (int i = 0; i < ylistTable.size(); i++) {
		ylistTable[i][ySize] = 0;
	}
	map<int, int>::const_iterator st = s.begin(), en = s.end();
	while (st != en) {
		ylistTable[(*st).first][ySize] = ylistTable[(*st).first][ySize]
		    + (*st).second;
		++st;
	}
}

bool refPolynomial::operator==(const refPolynomial &p) const {
	return s == p.s;
}

void refPolynomial::push_back(int expo, int coeff) {
	s[expo] += coeff;
}

ostream &refPolynomial::output(ostream &out) const {
	map<int, int>::const_iterator st = s.begin(), en = s.end();
	out << "z(x) = ";
	while (st != en) {
		out << (*st).second << "x^" << (*st).first << " + ";
		++st;
	}
	out << "0";
	return out;
}

ostream & operator<<(ostream & out, const refPolynomial &p) {
    return p.output(out);
}

// }}}
#endif

