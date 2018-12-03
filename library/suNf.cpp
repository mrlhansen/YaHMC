#include <suN.h>
#include <random.h>
#include <cstdio>
#include <cmath>

suNf::suNf()
{
	for(int i = 0; i < REPR_DIM*REPR_DIM; i++)
	{
		re[i] = 0;
		im[i] = 0;
	}
}

suNf operator+(const suNf &lhs, const suNf &rhs)
{
	suNf ret(lhs);
	ret += rhs;
	return ret;
}

suNf& suNf::operator+=(const suNf &rhs)
{
	for(int i = 0; i < REPR_DIM*REPR_DIM; i++)
	{
		re[i] += rhs.re[i];
		im[i] += rhs.im[i];
	}
	return *this;
}

suNf operator-(const suNf &lhs, const suNf &rhs)
{
	suNf ret(lhs);
	ret -= rhs;
	return ret;
}

suNf& suNf::operator-=(const suNf &rhs)
{
	for(int i = 0; i < REPR_DIM*REPR_DIM; i++)
	{
		re[i] -= rhs.re[i];
		im[i] -= rhs.im[i];
	}
	return *this;
}

suNf operator*(const suNf &lhs, const suNf &rhs)
{
	suNf ret;

	for(int i = 0; i < REPR_DIM; i++)
	{
		for(int j = 0; j < REPR_DIM; j++)
		{
			for(int k = 0; k < REPR_DIM; k++)
			{
				ret.re[i*REPR_DIM+j] += lhs.re[i*REPR_DIM+k]*rhs.re[k*REPR_DIM+j] - lhs.im[i*REPR_DIM+k]*rhs.im[k*REPR_DIM+j];
				ret.im[i*REPR_DIM+j] += lhs.re[i*REPR_DIM+k]*rhs.im[k*REPR_DIM+j] + lhs.im[i*REPR_DIM+k]*rhs.re[k*REPR_DIM+j];
			}
		}
	}

	return ret;
}

suNf& suNf::operator*=(const suNf &rhs)
{
	suNf ret;

	for(int i = 0; i < REPR_DIM; i++)
	{
		for(int j = 0; j < REPR_DIM; j++)
		{
			for(int k = 0; k < REPR_DIM; k++)
			{
				ret.re[i*REPR_DIM+j] += re[i*REPR_DIM+k]*rhs.re[k*REPR_DIM+j] - im[i*REPR_DIM+k]*rhs.im[k*REPR_DIM+j];
				ret.im[i*REPR_DIM+j] += re[i*REPR_DIM+k]*rhs.im[k*REPR_DIM+j] + im[i*REPR_DIM+k]*rhs.re[k*REPR_DIM+j];
			}
		}
	}

	*this = ret;
	return *this;
}

suNf operator*(const complex &lhs, const suNf &rhs)
{
	suNf ret(rhs);
	ret *= lhs;
	return ret;
}

suNf& suNf::operator*=(const complex &rhs)
{
	double rtmp;

	for(int i = 0; i < REPR_DIM*REPR_DIM; i++)
	{
		rtmp = re[i];
		re[i] = rtmp*rhs.re - im[i]*rhs.im;
		im[i] = rtmp*rhs.im + im[i]*rhs.re;
	}

	return *this;
}

suNf operator*(const double &lhs, const suNf &rhs)
{
	suNf ret(rhs);
	ret *= lhs;
	return ret;
}

suNf& suNf::operator*=(const double &rhs)
{
	for(int i = 0; i < REPR_DIM*REPR_DIM; i++)
	{
		re[i] *= rhs;
		im[i] *= rhs;
	}
	return *this;
}

double suNf::sqnorm()
{
	double ret = 0;
	for(int i = 0; i < REPR_DIM*REPR_DIM; i++)
	{
		ret += re[i]*re[i] + im[i]*im[i];
	}
	return ret;
}

suNf dagger(suNf g)
{
	suNf ret;
	for(int i = 0; i < REPR_DIM; i++)
	{
		for(int j = 0; j < REPR_DIM; j++)
		{
			ret.re[i*REPR_DIM+j] =  g.re[j*REPR_DIM+i];
			ret.im[i*REPR_DIM+j] = -g.im[j*REPR_DIM+i];
		}
	}
	return ret;
}

double trace(suNf g)
{
	double ret = 0;
	for(int i = 0; i < REPR_DIM; i++)
	{
		ret += g.re[i*REPR_DIM+i];
	}
	return ret;
}

double trace(suNf &lhs, suNf &rhs)
{
	double re = 0;
	for(int i = 0; i < REPR_DIM; i++)
	{
		for(int k = 0; k < REPR_DIM; k++)
		{
			re += lhs.re[i*REPR_DIM+k]*rhs.re[k*REPR_DIM+i] - lhs.im[i*REPR_DIM+k]*rhs.im[k*REPR_DIM+i];
		}
	}
	return re;
}
