#include <spinor.h>
#include <random.h>
#include <cmath>
#include <simd.h>

Spinor::Spinor()
{
	for(int i = 0; i < REPR_DIM; i++)
	{
		re[i][0] = 0;
		re[i][1] = 0;
		re[i][2] = 0;
		re[i][3] = 0;
		im[i][0] = 0;
		im[i][1] = 0;
		im[i][2] = 0;
		im[i][3] = 0;
	}
}

Spinor operator*(const suNf &lhs, const Spinor &rhs)
{
	Spinor u;

#ifdef ENABLE_AVX
	avx_vector res_re;
	avx_vector res_im;
	avx_vector rhs_re[REPR_DIM];
	avx_vector rhs_im[REPR_DIM];
	avx_vector lhs_re;
	avx_vector lhs_im;
	avx_vector tmp;

	for(int i = 0; i < REPR_DIM; i++)
	{
		avx_load(rhs_re[i], rhs.re[i]);
		avx_load(rhs_im[i], rhs.im[i]);
	}

	for(int i = 0; i < REPR_DIM; i++)
	{
		avx_set_zero(res_re);
		avx_set_zero(res_im);

		for(int j = 0; j < REPR_DIM; j++)
		{
			avx_set_dbl(lhs_re, lhs.re[i*REPR_DIM+j]);
			avx_set_dbl(lhs_im, lhs.im[i*REPR_DIM+j]);

			avx_mul(tmp, lhs_re, rhs_re[j]);
			avx_add_assign(res_re, tmp);
			avx_mul(tmp, lhs_im, rhs_im[j]);
			avx_sub_assign(res_re, tmp);

			avx_mul(tmp, lhs_re, rhs_im[j]);
			avx_add_assign(res_im, tmp);
			avx_mul(tmp, lhs_im, rhs_re[j]);
			avx_add_assign(res_im, tmp);
		}

		avx_store(u.re[i], res_re);
		avx_store(u.im[i], res_im);
	}
#else
	for(int i = 0; i < REPR_DIM; i++)
	{
		for(int j = 0; j < REPR_DIM; j++)
		{
			for(int k = 0; k < 4; k++)
			{
				u.re[i][k] += lhs.re[i*REPR_DIM+j]*rhs.re[j][k] - lhs.im[i*REPR_DIM+j]*rhs.im[j][k];
				u.im[i][k] += lhs.re[i*REPR_DIM+j]*rhs.im[j][k] + lhs.im[i*REPR_DIM+j]*rhs.re[j][k];
			}
		}
	}
#endif

	return u;
}

Spinor operator*(const double &lhs, const Spinor &rhs)
{
	Spinor u(rhs);
	u *= lhs;
	return u;
}

Spinor& Spinor::operator*=(const double &rhs)
{
	for(int i = 0; i < REPR_DIM; i++)
	{
		re[i][0] *= rhs;
		re[i][1] *= rhs;
		re[i][2] *= rhs;
		re[i][3] *= rhs;
		im[i][0] *= rhs;
		im[i][1] *= rhs;
		im[i][2] *= rhs;
		im[i][3] *= rhs;
	}
	return *this;
}

Spinor operator*(const complex &lhs, const Spinor &rhs)
{
	Spinor u(rhs);
	u *= lhs;
	return u;
}

Spinor& Spinor::operator*=(const complex &rhs)
{
	double rtmp;

	for(int i = 0; i < REPR_DIM; i++)
	{
		for(int k = 0; k < 4; k++)
		{
			rtmp = re[i][k];
			re[i][k] = rtmp*rhs.re - im[i][k]*rhs.im;
			im[i][k] = rtmp*rhs.im + im[i][k]*rhs.re;
		}
	}

	return *this;
}

Spinor operator+(const Spinor &lhs, const Spinor &rhs)
{
	Spinor u(lhs);
	u += rhs;
	return u;
}

Spinor& Spinor::operator+=(const Spinor &rhs)
{
	for(int i = 0; i < REPR_DIM; i++)
	{
		re[i][0] += rhs.re[i][0];
		re[i][1] += rhs.re[i][1];
		re[i][2] += rhs.re[i][2];
		re[i][3] += rhs.re[i][3];
		im[i][0] += rhs.im[i][0];
		im[i][1] += rhs.im[i][1];
		im[i][2] += rhs.im[i][2];
		im[i][3] += rhs.im[i][3];
	}

	return *this;
}

Spinor operator-(const Spinor &lhs, const Spinor &rhs)
{
	Spinor u(lhs);
	u -= rhs;
	return u;
}

Spinor& Spinor::operator-=(const Spinor &rhs)
{
	for(int i = 0; i < REPR_DIM; i++)
	{
		re[i][0] -= rhs.re[i][0];
		re[i][1] -= rhs.re[i][1];
		re[i][2] -= rhs.re[i][2];
		re[i][3] -= rhs.re[i][3];
		im[i][0] -= rhs.im[i][0];
		im[i][1] -= rhs.im[i][1];
		im[i][2] -= rhs.im[i][2];
		im[i][3] -= rhs.im[i][3];
	}

	return *this;
}

double Spinor::sqnorm()
{
	double res = 0;

	for(int i = 0; i < REPR_DIM; i++)
	{
		res += re[i][0]*re[i][0] + im[i][0]*im[i][0];
		res += re[i][1]*re[i][1] + im[i][1]*im[i][1];
		res += re[i][2]*re[i][2] + im[i][2]*im[i][2];
		res += re[i][3]*re[i][3] + im[i][3]*im[i][3];
	}

	return res;
}

void Spinor::gamma0()
{
	Spinor u(*this);
	for(int i = 0; i < REPR_DIM; i++)
	{
		re[i][0] = -u.re[i][2];
		re[i][1] = -u.re[i][3];
		re[i][2] = -u.re[i][0];
		re[i][3] = -u.re[i][1];
		im[i][0] = -u.im[i][2];
		im[i][1] = -u.im[i][3];
		im[i][2] = -u.im[i][0];
		im[i][3] = -u.im[i][1];
	}
}

void Spinor::gamma1()
{
	Spinor u(*this);
	for(int i = 0; i < REPR_DIM; i++)
	{
		re[i][0] =  u.im[i][3];
		re[i][1] =  u.im[i][2];
		re[i][2] = -u.im[i][1];
		re[i][3] = -u.im[i][0];
		im[i][0] = -u.re[i][3];
		im[i][1] = -u.re[i][2];
		im[i][2] =  u.re[i][1];
		im[i][3] =  u.re[i][0];
	}
}

void Spinor::gamma2()
{
	Spinor u(*this);
	for(int i = 0; i < REPR_DIM; i++)
	{
		re[i][0] = -u.re[i][3];
		re[i][1] =  u.re[i][2];
		re[i][2] =  u.re[i][1];
		re[i][3] = -u.re[i][0];
		im[i][0] = -u.im[i][3];
		im[i][1] =  u.im[i][2];
		im[i][2] =  u.im[i][1];
		im[i][3] = -u.im[i][0];
	}
}

void Spinor::gamma3()
{
	Spinor u(*this);
	for(int i = 0; i < REPR_DIM; i++)
	{
		re[i][0] =  u.im[i][2];
		re[i][1] = -u.im[i][3];
		re[i][2] = -u.im[i][0];
		re[i][3] =  u.im[i][1];
		im[i][0] = -u.re[i][2];
		im[i][1] =  u.re[i][3];
		im[i][2] =  u.re[i][0];
		im[i][3] = -u.re[i][1];
	}
}

void Spinor::gamma5()
{
	for(int i = 0; i < REPR_DIM; i++)
	{
		re[i][2] *= -1;
		im[i][2] *= -1;
		re[i][3] *= -1;
		im[i][3] *= -1;
	}
}

void Spinor::g5_sigma(int mu, int nu)
{
	Spinor u(*this);

	if(mu == 0 && nu == 1)
	{
		for(int i = 0; i < REPR_DIM; i++)
		{
			re[i][0] =  u.im[i][1];
			re[i][1] =  u.im[i][0];
			re[i][2] =  u.im[i][3];
			re[i][3] =  u.im[i][2];
			im[i][0] = -u.re[i][1];
			im[i][1] = -u.re[i][0];
			im[i][2] = -u.re[i][3];
			im[i][3] = -u.re[i][2];
		}
	}
	else if(nu == 0 && mu == 1)
	{
		for(int i = 0; i < REPR_DIM; i++)
		{
			re[i][0] = -u.im[i][1];
			re[i][1] = -u.im[i][0];
			re[i][2] = -u.im[i][3];
			re[i][3] = -u.im[i][2];
			im[i][0] =  u.re[i][1];
			im[i][1] =  u.re[i][0];
			im[i][2] =  u.re[i][3];
			im[i][3] =  u.re[i][2];
		}
	}
	else if(mu == 0 && nu == 2)
	{
		for(int i = 0; i < REPR_DIM; i++)
		{
			re[i][0] = -u.re[i][1];
			re[i][1] =  u.re[i][0];
			re[i][2] = -u.re[i][3];
			re[i][3] =  u.re[i][2];
			im[i][0] = -u.im[i][1];
			im[i][1] =  u.im[i][0];
			im[i][2] = -u.im[i][3];
			im[i][3] =  u.im[i][2];
		}
	}
	else if(nu == 0 && mu == 2)
	{
		for(int i = 0; i < REPR_DIM; i++)
		{
			re[i][0] =  u.re[i][1];
			re[i][1] = -u.re[i][0];
			re[i][2] =  u.re[i][3];
			re[i][3] = -u.re[i][2];
			im[i][0] =  u.im[i][1];
			im[i][1] = -u.im[i][0];
			im[i][2] =  u.im[i][3];
			im[i][3] = -u.im[i][2];
		}
	}
	else if(mu == 0 && nu == 3)
	{
		for(int i = 0; i < REPR_DIM; i++)
		{
			re[i][0] =  u.im[i][0];
			re[i][1] = -u.im[i][1];
			re[i][2] =  u.im[i][2];
			re[i][3] = -u.im[i][3];
			im[i][0] = -u.re[i][0];
			im[i][1] =  u.re[i][1];
			im[i][2] = -u.re[i][2];
			im[i][3] =  u.re[i][3];
		}
	}
	else if(nu == 0 && mu == 3)
	{
		for(int i = 0; i < REPR_DIM; i++)
		{
			re[i][0] = -u.im[i][0];
			re[i][1] =  u.im[i][1];
			re[i][2] = -u.im[i][2];
			re[i][3] =  u.im[i][3];
			im[i][0] =  u.re[i][0];
			im[i][1] = -u.re[i][1];
			im[i][2] =  u.re[i][2];
			im[i][3] = -u.re[i][3];
		}
	}
	else if(mu == 1 && nu == 2)
	{
		for(int i = 0; i < REPR_DIM; i++)
		{
			re[i][0] = -u.im[i][0];
			re[i][1] =  u.im[i][1];
			re[i][2] =  u.im[i][2];
			re[i][3] = -u.im[i][3];
			im[i][0] =  u.re[i][0];
			im[i][1] = -u.re[i][1];
			im[i][2] = -u.re[i][2];
			im[i][3] =  u.re[i][3];
		}
	}
	else if(nu == 1 && mu == 2)
	{
		for(int i = 0; i < REPR_DIM; i++)
		{
			re[i][0] =  u.im[i][0];
			re[i][1] = -u.im[i][1];
			re[i][2] = -u.im[i][2];
			re[i][3] =  u.im[i][3];
			im[i][0] = -u.re[i][0];
			im[i][1] =  u.re[i][1];
			im[i][2] =  u.re[i][2];
			im[i][3] = -u.re[i][3];
		}
	}
	else if(mu == 1 && nu == 3)
	{
		for(int i = 0; i < REPR_DIM; i++)
		{
			re[i][0] = -u.re[i][1];
			re[i][1] =  u.re[i][0];
			re[i][2] =  u.re[i][3];
			re[i][3] = -u.re[i][2];
			im[i][0] = -u.im[i][1];
			im[i][1] =  u.im[i][0];
			im[i][2] =  u.im[i][3];
			im[i][3] = -u.im[i][2];
		}
	}
	else if(nu == 1 && mu == 3)
	{
		for(int i = 0; i < REPR_DIM; i++)
		{
			re[i][0] =  u.re[i][1];
			re[i][1] = -u.re[i][0];
			re[i][2] = -u.re[i][3];
			re[i][3] =  u.re[i][2];
			im[i][0] =  u.im[i][1];
			im[i][1] = -u.im[i][0];
			im[i][2] = -u.im[i][3];
			im[i][3] =  u.im[i][2];
		}
	}
	else if(mu == 2 && nu == 3)
	{
		for(int i = 0; i < REPR_DIM; i++)
		{
			re[i][0] = -u.im[i][1];
			re[i][1] = -u.im[i][0];
			re[i][2] =  u.im[i][3];
			re[i][3] =  u.im[i][2];
			im[i][0] =  u.re[i][1];
			im[i][1] =  u.re[i][0];
			im[i][2] = -u.re[i][3];
			im[i][3] = -u.re[i][2];
		}
	}
	else if(nu == 2 && mu == 3)
	{
		for(int i = 0; i < REPR_DIM; i++)
		{
			re[i][0] =  u.im[i][1];
			re[i][1] =  u.im[i][0];
			re[i][2] = -u.im[i][3];
			re[i][3] = -u.im[i][2];
			im[i][0] = -u.re[i][1];
			im[i][1] = -u.re[i][0];
			im[i][2] =  u.re[i][3];
			im[i][3] =  u.re[i][2];
		}
	}
}

void Spinor::g5_one_minus_gamma(int mu)
{
	Spinor u(*this);

	if(mu == 0)
	{
		for(int i = 0; i < REPR_DIM; i++)
		{
			re[i][0] =  u.re[i][0] + u.re[i][2];
			re[i][1] =  u.re[i][1] + u.re[i][3];
			re[i][2] = -u.re[i][2] - u.re[i][0];
			re[i][3] = -u.re[i][3] - u.re[i][1];
			im[i][0] =  u.im[i][0] + u.im[i][2];
			im[i][1] =  u.im[i][1] + u.im[i][3];
			im[i][2] = -u.im[i][2] - u.im[i][0];
			im[i][3] = -u.im[i][3] - u.im[i][1];
		}
	}
	else if(mu == 1)
	{
		for(int i = 0; i < REPR_DIM; i++)
		{
			re[i][0] =  u.re[i][0] - u.im[i][3];
			re[i][1] =  u.re[i][1] - u.im[i][2];
			re[i][2] = -u.re[i][2] - u.im[i][1];
			re[i][3] = -u.re[i][3] - u.im[i][0];
			im[i][0] =  u.im[i][0] + u.re[i][3];
			im[i][1] =  u.im[i][1] + u.re[i][2];
			im[i][2] = -u.im[i][2] + u.re[i][1];
			im[i][3] = -u.im[i][3] + u.re[i][0];
		}
	}
	else if(mu == 2)
	{
		for(int i = 0; i < REPR_DIM; i++)
		{
			re[i][0] =  u.re[i][0] + u.re[i][3];
			re[i][1] =  u.re[i][1] - u.re[i][2];
			re[i][2] = -u.re[i][2] + u.re[i][1];
			re[i][3] = -u.re[i][3] - u.re[i][0];
			im[i][0] =  u.im[i][0] + u.im[i][3];
			im[i][1] =  u.im[i][1] - u.im[i][2];
			im[i][2] = -u.im[i][2] + u.im[i][1];
			im[i][3] = -u.im[i][3] - u.im[i][0];
		}
	}
	else if(mu == 3)
	{
		for(int i = 0; i < REPR_DIM; i++)
		{
			re[i][0] =  u.re[i][0] - u.im[i][2];
			re[i][1] =  u.re[i][1] + u.im[i][3];
			re[i][2] = -u.re[i][2] - u.im[i][0];
			re[i][3] = -u.re[i][3] + u.im[i][1];
			im[i][0] =  u.im[i][0] + u.re[i][2];
			im[i][1] =  u.im[i][1] - u.re[i][3];
			im[i][2] = -u.im[i][2] + u.re[i][0];
			im[i][3] = -u.im[i][3] - u.re[i][1];
		}
	}
}

void Spinor::p_minus()
{
	Spinor u(*this);
	u *= 0.5;

	for(int i = 0; i < REPR_DIM; i++)
	{
		re[i][0] = u.re[i][0] + u.re[i][2];
		re[i][1] = u.re[i][1] + u.re[i][3];
		re[i][2] = u.re[i][2] + u.re[i][0];
		re[i][3] = u.re[i][3] + u.re[i][1];
		im[i][0] = u.im[i][0] + u.im[i][2];
		im[i][1] = u.im[i][1] + u.im[i][3];
		im[i][2] = u.im[i][2] + u.im[i][0];
		im[i][3] = u.im[i][3] + u.im[i][1];
	}
}

void Spinor::p_plus()
{
	Spinor u(*this);
	u *= 0.5;

	for(int i = 0; i < REPR_DIM; i++)
	{
		re[i][0] = u.re[i][0] - u.re[i][2];
		re[i][1] = u.re[i][1] - u.re[i][3];
		re[i][2] = u.re[i][2] - u.re[i][0];
		re[i][3] = u.re[i][3] - u.re[i][1];
		im[i][0] = u.im[i][0] - u.im[i][2];
		im[i][1] = u.im[i][1] - u.im[i][3];
		im[i][2] = u.im[i][2] - u.im[i][0];
		im[i][3] = u.im[i][3] - u.im[i][1];
	}
}

void Spinor::randomize()
{
	for(int i = 0; i < REPR_DIM; i++)
	{
		re[i][0] = rand_gaussian();
		im[i][0] = rand_gaussian();
		re[i][1] = rand_gaussian();
		im[i][1] = rand_gaussian();
		re[i][2] = rand_gaussian();
		im[i][2] = rand_gaussian();
		re[i][3] = rand_gaussian();
		im[i][3] = rand_gaussian();
	}
}

complex inner_product(Spinor &a, Spinor &b)
{
	complex c;

	for(int i = 0; i < REPR_DIM; i++)
	{
		c.re += a.re[i][0]*b.re[i][0] + a.im[i][0]*b.im[i][0];
		c.im += a.im[i][0]*b.re[i][0] - a.re[i][0]*b.im[i][0];
		c.re += a.re[i][1]*b.re[i][1] + a.im[i][1]*b.im[i][1];
		c.im += a.im[i][1]*b.re[i][1] - a.re[i][1]*b.im[i][1];
		c.re += a.re[i][2]*b.re[i][2] + a.im[i][2]*b.im[i][2];
		c.im += a.im[i][2]*b.re[i][2] - a.re[i][2]*b.im[i][2];
		c.re += a.re[i][3]*b.re[i][3] + a.im[i][3]*b.im[i][3];
		c.im += a.im[i][3]*b.re[i][3] - a.re[i][3]*b.im[i][3];
	}

	return c;
}

double inner_product_re(Spinor &a, Spinor &b)
{
	double ret = 0;

	for(int i = 0; i < REPR_DIM; i++)
	{
		ret += a.re[i][0]*b.re[i][0] + a.im[i][0]*b.im[i][0];
		ret += a.re[i][1]*b.re[i][1] + a.im[i][1]*b.im[i][1];
		ret += a.re[i][2]*b.re[i][2] + a.im[i][2]*b.im[i][2];
		ret += a.re[i][3]*b.re[i][3] + a.im[i][3]*b.im[i][3];
	}

	return ret;
}
