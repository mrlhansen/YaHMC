#ifndef SPINOR_H
#define SPINOR_H

#include <suN.h>
#include <complex.h>

struct Spinor {
	double re[REPR_DIM][4];
	double im[REPR_DIM][4];

	Spinor();

	Spinor& operator*=(const double&);
	Spinor& operator*=(const complex&);
	friend Spinor operator*(const double&, const Spinor&);
	friend Spinor operator*(const complex&, const Spinor&);
	friend Spinor operator*(const suNf&, const Spinor&);

	Spinor& operator-=(const Spinor&);
	friend Spinor operator-(const Spinor&, const Spinor&);

	Spinor& operator+=(const Spinor&);
	friend Spinor operator+(const Spinor&, const Spinor&);

	void gamma0();
	void gamma1();
	void gamma2();
	void gamma3();
	void gamma5();

	void g5_sigma(int, int);
	void g5_one_minus_gamma(int);
	void p_plus();
	void p_minus();

	double sqnorm();
	void randomize();

	complex get(int s, int c)
	{
		return complex(re[c][s], im[c][s]);
	}
};

complex inner_product(Spinor&, Spinor&);
double inner_product_re(Spinor&, Spinor&);

#endif
