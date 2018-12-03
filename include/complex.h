#ifndef COMPLEX_H
#define COMPLEX_H

#include <cmath>

struct complex {
	double re;
	double im;

	complex()
	{
		re = 0;
		im = 0;
	}

	complex(double real, double imag)
	{
		re = real;
		im = imag;
	}

	friend complex operator+(const complex &lhs, const complex &rhs)
	{
		complex u(lhs);
		u += rhs;
		return u;
	}

	complex& operator+=(const complex &rhs)
	{
		re += rhs.re;
		im += rhs.im;
		return *this;
	}

	friend complex operator-(const complex &lhs, const complex &rhs)
	{
		complex u(lhs);
		u -= rhs;
		return u;
	}

	complex& operator-=(const complex &rhs)
	{
		re -= rhs.re;
		im -= rhs.im;
		return *this;
	}

	friend complex operator*(const int &lhs, const complex &rhs)
	{
		complex u(rhs);
		u *= lhs;
		return u;
	}

	complex& operator*=(const int &rhs)
	{
		re *= rhs;
		im *= rhs;
		return *this;
	}

	friend complex operator*(const double &lhs, const complex &rhs)
	{
		complex u(rhs);
		u *= lhs;
		return u;
	}

	friend complex operator*(const complex &lhs, const double &rhs)
	{
		complex u(lhs);
		u *= rhs;
		return u;
	}

	complex& operator*=(const double &rhs)
	{
		re *= rhs;
		im *= rhs;
		return *this;
	}

	friend complex operator*(const complex &lhs, const complex &rhs)
	{
		complex u(lhs);
		u *= rhs;
		return u;
	}

	complex& operator*=(const complex &rhs)
	{
		double rtmp = re;
		re = rtmp*rhs.re - im*rhs.im;
		im = rtmp*rhs.im + im*rhs.re;
		return *this;
	}

	friend complex operator/(const complex &lhs, const double &rhs)
	{
		complex u(lhs);
		u /= rhs;
		return u;
	}

	complex& operator/=(const double &rhs)
	{
		re /= rhs;
		im /= rhs;
		return *this;
	}

	friend complex operator/(const complex &lhs, const complex &rhs)
	{
		complex u(lhs);
		u /= rhs;
		return u;
	}

	complex& operator/=(const complex &rhs)
	{
		double rtmp = re;
		double sq = rhs.re*rhs.re + rhs.im*rhs.im;
		re = (rtmp*rhs.re + im*rhs.im) / sq;
		im = (im*rhs.re - rtmp*rhs.im) / sq;
		return *this;
	}
};

#define cabs(v) \
	sqrt(v.re*v.re + v.im*v.im)

#define cabs2(v) \
	(v.re*v.re + v.im*v.im)

#define conj(v) \
	complex(v.re, -v.im)

#endif
