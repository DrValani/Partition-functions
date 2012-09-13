// {{{ Polynomial Class

#ifndef Polynomial_CC
#define Polynomial_CC

#include "Polynomial.h"
using namespace std;

/** 
 * This class stores the polynomials using mpz ints. The size of each 
 * polynomial is the max power of the lattice.
 */

Polynomial::Polynomial(const int polynomialSize) : 	MAXPOLYSIZE(polynomialSize),degree(0){
	 summands = new mpz_class[polynomialSize]; summands[0] = 0;
}
Polynomial::Polynomial(const Polynomial &p):
	MAXPOLYSIZE (p.MAXPOLYSIZE), degree(p.degree){
	summands = new mpz_class[MAXPOLYSIZE];
	for(int i = 0; i <= degree; i++){
		summands[i] = p.summands[i];
	}
}
/**
 * This is the same as this = this + (p * ct)
 * Multiply a given polynomial by an int then add it.
 */
void Polynomial::multiPlus(const Polynomial &p, const int &ct) {
	if (ct == 1) {
		for (int i = 0; i <= p.getDegree(); i++) {
			this->push_back(p[i], i);
		}
	} else
		for (int i = 0; i <= p.getDegree(); i++) {
			this->push_back(p[i] * ct, i);
		}
}

void Polynomial::clear() {
	this->degree = 0;
	//  for(int i = 0; i < MAXPOLYSIZE; i++)
	summands[0] = 0;
}

template<class T>
void Polynomial::push_back(const T &coeff, const int &expo) {
	if (isInBounds(expo)) {

	if (expo > this->degree) {
		for (int i = this->degree + 1; i <= expo; i++)
			this->summands[i] = 0;
		this->degree = expo;
	}
	this->summands[expo] = (this->summands[expo] + coeff);
	}else
		throw new runtime_error("Increase MAXPOLYSIZE");

	}

template<class T>
void Polynomial::push_back(const T *coeff, const int &expo) {

	if (isInBounds(expo)) {
		if (expo > this->degree) {
			for (int i = this->degree + 1; i <= expo; i++)
				this->summands[i] = 0;
			this->degree = expo;
		}
		this->summands[expo] = (this->summands[expo] + *coeff);
	}
	else
		throw new runtime_error("Increase MAXPOLYSIZE");

}

void Polynomial::push_back(const int coeff, const int expo) {
	if (isInBounds(expo)) {
		if (expo > this->degree) {
			for (int i = this->degree + 1; i <= expo; i++)
				this->summands[i] = 0;
			this->degree = expo;
		}
		this->summands[expo] = (this->summands[expo] + coeff);

	} else
		throw new runtime_error("Increase MAXPOLYSIZE");
}

/**
 * This function raises polynomial p to the power of pow, and then adds this 
 * new polynoimial to *this.
 */
void Polynomial::plusPower(const Polynomial &p, const int pow) {
	// const mpz_class *ip;
	//ip = p.summands;
	for (int i = 0; i <= p.getDegree(); i++) {
		this->push_back(p[i], i + pow);
		//++ip;
	}
	//  ip = 0;
	//  delete ip;
}
void Polynomial::operator+=(const Polynomial &p) {
	const mpz_class *endp, *ip;
	ip = p.summands;
	endp = p.summands + 1 + p.degree;
	int i = 0;
	for (; ip != endp; ++ip) {
		this->push_back(ip, i);
		i++;
	}
	endp = 0;
	ip = 0;
	delete endp;
	delete ip;
}

Polynomial Polynomial::operator=(const Polynomial &p) {
	this->degree = p.getDegree();
	for (int i = p.getDegree(); i >= 0; i--)
		this->summands[i] = p[i];
	return *this;
}

ostream & operator<<(ostream & out, const Polynomial &p) {
	for (int i = 0; i <= p.getDegree(); i++) {
		out << p[i];
		out << endl;
	}
	return out;
}

ostream & operator,(ostream & out, const Polynomial &p) {
	int i = p.getDegree();
	out << "z(x) = ";
	out << p[i] << "x^{" << i << "}";

	for (i--; i >= 0; i--)
		if (p[i] != 0) {
			out << " + " << p[i];
			out << "x^{" << i << "}";
		}
	return out;
}

#endif

// }}}

