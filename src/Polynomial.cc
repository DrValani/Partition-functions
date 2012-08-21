// {{{ Polynomial Class

#ifndef Polynomial_CC
#define Polynomial_CC

#ifndef MAXPOLYSIZE
#define MAXPOLYSIZE (3 * xSize * ySize * zSize) + 1// Max power of lattice
#endif


#include "Polynomial.h"
using namespace std;





/** 
 * This class stores the polynomials using mpz ints. The size of each 
 * polynomial is the max power of the lattice.
 */
/**
 * This is the same as this = this + (p * ct)
 * Multiply a given polynomial by an int then add it.
 */
void Polynomial::multiPlus(const Polynomial &p, const int &ct) {
	if (ct == 1) {
		for (int i = 0; i <= p.getDegree(); i++) {
			this->push_back(p[i], i);
		}
	} else for (int i = 0; i <= p.getDegree(); i++) {
		this->push_back(p[i] * ct, i);
	}
}

void Polynomial::clear() {
	this->degree = 0;
	//  for(int i = 0; i < MAXPOLYSIZE; i++)
	summands[0] = 0;
}

Polynomial Polynomial::operator=(const Polynomial &p) {
	this->degree = p.getDegree();
	for (int i = p.getDegree(); i >= 0; i--)
		this->summands[i] = p[i];
	return *this;
}
/**
 * Most polynomial's if they differ, are most likely to have similar coeffients
 * on the ends but different coeffs in the middle. So search from the middle 
 * outwards.
 *
 */
bool Polynomial::operator==(const Polynomial &p) const {
	if (this->degree == p.getDegree()) {
		int halfway = this->degree / 2;
		for (int i = halfway; i >= 0; i--) {
			if (this->summands[i] != p[i]) return false;
		}
		for (int i = halfway + 1; i <= this->degree; i++) {
			if (this->summands[i] != p[i]) return false;
		}
		return true;
	}
	return false;
}

bool Polynomial::equals(const Polynomial &p) const {
	if (this->degree == p.getDegree()) {
		for (int i = this->degree; i >= 0; i--) {
			if (this->summands[i] != p[i]) return false;
		}
		return true;
	}
	return false;
}

template<class T>
void Polynomial::push_back(const T &coeff, const int &expo) {
	if (expo >= MAXPOLYSIZE) {
		cout << "error 1: increase MAXPOLYSIZE" << endl;
		cout << MAXPOLYSIZE << endl;
		cout << expo << endl;
		cout << this->degree << endl;
		cout, *this;
		throw new runtime_error("error 1: increase MAXPOLYSIZE");
	}
	if (expo > this->degree) {
		for (int i = this->degree + 1; i <= expo; i++)
			this->summands[i] = 0;
		this->degree = expo;
	}
	this->summands[expo] = (this->summands[expo] + coeff);
}

template<class T>
void Polynomial::push_back(const T *coeff, const int &expo) {
	if (expo >= MAXPOLYSIZE) {
		cout << "error 2: increase MAXPOLYSIZE" << endl;
		cout << MAXPOLYSIZE << endl;
		cout << expo + this->degree << endl;
		throw new runtime_error("error 2: increase MAXPOLYSIZE");
	}

	if (expo > this->degree) {
		for (int i = this->degree + 1; i <= expo; i++)
			this->summands[i] = 0;
		this->degree = expo;
	}
	this->summands[expo] = (this->summands[expo] + *coeff);
}
void Polynomial::push_back(const int coeff, const int expo) {
	if (expo >= MAXPOLYSIZE) {
		cout << "error 3: increase MAXPOLYSIZE" << endl;
		throw new runtime_error("error 3: increase MAXPOLYSIZE" );
	}
	if (expo > this->degree) {
		for (int i = this->degree + 1; i <= expo; i++)
			this->summands[i] = 0;
		this->degree = expo;
	}
	this->summands[expo] = (this->summands[expo] + coeff);
}
/*
 // Multiply polynomial with an integer
 Polynomial Polynomial::operator*(const int &x)const{
 Polynomial temp = *this;
 for(int i = 0; i <= this->degree; i++)
 temp.summands[i] = (this->summands[i] * x);
 return temp;
 }
 */
Polynomial Polynomial::operator^(const int &x) const {
	Polynomial temp;
	temp.setDegree(this->getDegree() + x);
	for (int i = 0; i <= x; i++)
		temp.summands[i] = 0;
	if (temp.getDegree() >= MAXPOLYSIZE) {
		cout << "error 4: increase MAXPOLYSIZE" << endl;
		throw new runtime_error("error 4: increase MAXPOLYSIZE");
	}
	for (int i = 0; i <= this->degree; i++) {
		temp.summands[i + x] = this->summands[i];
	}
	return temp;
}
Polynomial Polynomial::operator+(const Polynomial &p) const {
	Polynomial temp;
	if (this->degree >= p.getDegree()) {
		temp.setDegree(this->degree);
		int i = 0;
		for (; i <= p.getDegree(); i++)
			temp.summands[i] = this->summands[i] + p.summands[i];
		for (; i <= this->getDegree(); i++)
			temp.summands[i] = this->summands[i];
	} else {
		temp.setDegree(p.getDegree());
		int i = 0;
		for (; i <= this->getDegree(); i++)
			temp.summands[i] = this->summands[i] + p.summands[i];
		for (; i <= p.getDegree(); i++)
			temp.summands[i] = p.summands[i];
	}
	return temp;
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

ostream & operator<<(ostream & out, const Polynomial &p) {
    for (int i = 0; i <= p.getDegree(); i++) {
        out << p[i];
        out << endl;
    }
    return out;
}

ostream & operator, (ostream & out, const Polynomial &p) {
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

