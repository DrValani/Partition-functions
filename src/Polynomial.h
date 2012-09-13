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

using namespace std;


class Polynomial {
public:
	Polynomial(const int polynomialSize);
	Polynomial(const Polynomial &p);
	~Polynomial(){delete[] summands;}
	Polynomial operator=(const Polynomial &p);
	Polynomial operator+(const Polynomial &p) const;
	Polynomial operator-(const Polynomial &p);
	Polynomial operator*(const Polynomial &p);  //Mulitply two polynomials
	Polynomial operator*(const int &x) const;  // Mutiply polynomial with integer
	void operator+=(const Polynomial &p);
	void clear();
	void push_back(const int coeff, const int expo);
	template<class T>
	void push_back(const T &coeff, const int &expo);
	template<class T>
	void push_back(const T *coeff, const int &expo);
	mpz_class operator[](int i) const {return summands[i];}
	int getDegree() const {return degree;}
	void plusPower(const Polynomial &p, const int pow);
	bool isZero() {return ((degree == 0) && (summands[0] == 0));}
	void multiPlus(const Polynomial &p, const int &ct);
private:
	const int MAXPOLYSIZE;
	const bool isInBounds(const int &expo){return expo < MAXPOLYSIZE;}
	mpz_class * summands;
	int degree;
};

#endif /* POLYNOMIAL_H_ */
