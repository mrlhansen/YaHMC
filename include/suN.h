#ifndef SUN_H
#define SUN_H

#include <cstdio>
#include <complex.h>
#include <repr.h>

struct suNa {
	double av[NG];

	suNa();
	suNa(int);
	double sqnorm();
	void random();

	suNa& operator+=(const suNa&);
	friend suNa operator+(const suNa&, const suNa&);

	suNa& operator-=(const suNa&);
	friend suNa operator-(const suNa&, const suNa&);

	suNa& operator*=(const double&);
	friend suNa operator*(const double&, const suNa&);
};

struct suNg {
	double re[NC*NC];
	double im[NC*NC];

	suNg();
	suNg(suNa&);

	void reunitarize();
	void unit();
	void random();
	double sqnorm();

	int read(FILE*);
	int write(FILE*);

	suNg& operator+=(const suNg&);
	friend suNg operator+(const suNg&, const suNg&);

	suNg& operator-=(const suNg&);
	friend suNg operator-(const suNg&, const suNg&);

	suNg& operator*=(const suNg&);
	friend suNg operator*(const suNg&, const suNg&);

	suNg& operator*=(const complex&);
	friend suNg operator*(const complex&, const suNg&);

	suNg& operator*=(const double&);
	friend suNg operator*(const double&, const suNg&);
};

struct suNf {
	double re[REPR_DIM*REPR_DIM];
	double im[REPR_DIM*REPR_DIM];

	suNf();
	double sqnorm();

	suNf& operator+=(const suNf&);
	friend suNf operator+(const suNf&, const suNf&);

	suNf& operator-=(const suNf&);
	friend suNf operator-(const suNf&, const suNf&);

	suNf& operator*=(const suNf&);
	friend suNf operator*(const suNf&, const suNf&);

	suNf& operator*=(const complex&);
	friend suNf operator*(const complex&, const suNf&);

	suNf& operator*=(const double&);
	friend suNf operator*(const double&, const suNf&);
};

suNa algebra_project(suNg&, double = 1.f);
suNa algebra_project(suNf&, double = 1.f);
suNg exp(suNa&, double = 1.f);

suNg dagger(suNg);
suNg transpose(suNg);
double trace(suNg);
double trace(suNg&, suNg&);
complex ctrace(suNg);
suNg exp(suNg&);

suNf dagger(suNf);
double trace(suNf);
double trace(suNf&, suNf&);

#endif
