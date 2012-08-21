/*
 * Polynomial.h
 *
 *  Created on: 20 Aug 2012
 *      Author: yogi
 */

#ifndef POLYNOMIAL_H_
#define POLYNOMIAL_H_
#include <gmp.h>
#include <gmpxx.h>
#include <iostream>
#include <vector>
#include <sstream>
#include <fstream>
#include "latticeSize.cc"
using namespace std;

class Polynomial {
public:
	Polynomial() :
			degree(0) {
			summands[0] = 0;
	}
	Polynomial operator=(const Polynomial &p);
	Polynomial operator+(const Polynomial &p) const;
	Polynomial operator-(const Polynomial &p);
	Polynomial operator*(const Polynomial &p);  //Mulitply two polynomials
	Polynomial operator*(const int &x) const;  // Mutiply polynomial with integer
	Polynomial operator^(const int &x) const;
	bool operator==(const Polynomial &p) const;
	void operator+=(const Polynomial &p);
	bool equals(const Polynomial &p) const;
	void clear();
	void push_back(const int coeff, const int expo);
	template<class T>
	void push_back(const T &coeff, const int &expo);
	template<class T>
	void push_back(const T *coeff, const int &expo);
	mpz_class operator[](int i) const {
		return summands[i];
	}
	int getDegree() const {
		return degree;
	}
	void setDegree(int x) {
		degree = x;
	}
	void plusPower(const Polynomial &p, const int pow);
	bool isZero() {
		return ((degree == 0) && (summands[0] == 0));
	}
	void multiPlus(const Polynomial &p, const int &ct);
private:
	mpz_class summands[MAXPOLYSIZE];
	int degree;
};



#endif /* POLYNOMIAL_H_ */
