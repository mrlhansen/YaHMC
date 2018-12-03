#include <suN.h>
#include <random.h>
#include <cstdio>
#include <cmath>

suNg::suNg()
{
	for(int i = 0; i < NC*NC; i++)
	{
		re[i] = 0;
		im[i] = 0;
	}
}

suNg::suNg(suNa &alg)
{
	double v;
	int n = 0;

	// Construct the anti-Hermitian generator
	for(int i = 0; i < NC; i++)
	{
		re[i*NC+i] = 0;
		im[i*NC+i] = 0;
	}

	for(int i = 0; i < NC; i++)
	{
		for(int j = 0; j < NC; j++)
		{
			if(i > j)
			{
				im[i*NC+j] = alg.av[n];
				im[j*NC+i] = alg.av[n];
				n++;
			}
			else if(i < j)
			{
				re[i*NC+j] = alg.av[n];
				re[j*NC+i] = -alg.av[n];
				n++;
			}
			else if(i > 0)
			{
				v = i*(i+1);
				v = sqrt(2.0/v);

				for(int k = 0; k < i; k++)
				{
					im[k*NC+k] += v*alg.av[n];
				}

				im[i*NC+i] -= i*v*alg.av[n];
				n++;
			}
		}
	}
}

suNg operator+(const suNg &lhs, const suNg &rhs)
{
	suNg ret(lhs);
	ret += rhs;
	return ret;
}

suNg& suNg::operator+=(const suNg &rhs)
{
	for(int i = 0; i < NC*NC; i++)
	{
		re[i] += rhs.re[i];
		im[i] += rhs.im[i];
	}
	return *this;
}

suNg operator-(const suNg &lhs, const suNg &rhs)
{
	suNg ret(lhs);
	ret -= rhs;
	return ret;
}

suNg& suNg::operator-=(const suNg &rhs)
{
	for(int i = 0; i < NC*NC; i++)
	{
		re[i] -= rhs.re[i];
		im[i] -= rhs.im[i];
	}
	return *this;
}

suNg operator*(const suNg &lhs, const suNg &rhs)
{
	suNg ret;

	for(int i = 0; i < NC; i++)
	{
		for(int j = 0; j < NC; j++)
		{
			for(int k = 0; k < NC; k++)
			{
				ret.re[i*NC+j] += lhs.re[i*NC+k]*rhs.re[k*NC+j] - lhs.im[i*NC+k]*rhs.im[k*NC+j];
				ret.im[i*NC+j] += lhs.re[i*NC+k]*rhs.im[k*NC+j] + lhs.im[i*NC+k]*rhs.re[k*NC+j];
			}
		}
	}

	return ret;
}

suNg& suNg::operator*=(const suNg &rhs)
{
	suNg ret;

	for(int i = 0; i < NC; i++)
	{
		for(int j = 0; j < NC; j++)
		{
			for(int k = 0; k < NC; k++)
			{
				ret.re[i*NC+j] += re[i*NC+k]*rhs.re[k*NC+j] - im[i*NC+k]*rhs.im[k*NC+j];
				ret.im[i*NC+j] += re[i*NC+k]*rhs.im[k*NC+j] + im[i*NC+k]*rhs.re[k*NC+j];
			}
		}
	}

	*this = ret;
	return *this;
}

suNg operator*(const complex &lhs, const suNg &rhs)
{
	suNg ret(rhs);
	ret *= lhs;
	return ret;
}

suNg& suNg::operator*=(const complex &rhs)
{
	double rtmp;

	for(int i = 0; i < NC*NC; i++)
	{
		rtmp = re[i];
		re[i] = rtmp*rhs.re - im[i]*rhs.im;
		im[i] = rtmp*rhs.im + im[i]*rhs.re;
	}

	return *this;
}

suNg operator*(const double &lhs, const suNg &rhs)
{
	suNg ret(rhs);
	ret *= lhs;
	return ret;
}

suNg& suNg::operator*=(const double &rhs)
{
	for(int i = 0; i < NC*NC; i++)
	{
		re[i] *= rhs;
		im[i] *= rhs;
	}
	return *this;
}

void suNg::reunitarize()
{
	#if (NC == 2)
	double norm;

	// Normalize the first row
	norm  = re[0]*re[0] + im[0]*im[0];
	norm += re[1]*re[1] + im[1]*im[1];
	norm  = sqrt(norm);

	if(norm == 0)
	{
		return;
	}

	re[0] /= norm;
	im[0] /= norm;
	re[1] /= norm;
	im[1] /= norm;

	// Construct the second row
	re[2] = -re[1];
	im[2] =  im[1];
	re[3] =  re[0];
	im[3] = -im[0];

	#elif (NC == 3)
	double norm, rtmp, itmp;

	// Normalize the first row
	norm  = re[0]*re[0] + im[0]*im[0];
	norm += re[1]*re[1] + im[1]*im[1];
	norm += re[2]*re[2] + im[2]*im[2];
	norm  = sqrt(norm);

	if(norm == 0)
	{
		return;
	}

	re[0] /= norm;
	im[0] /= norm;
	re[1] /= norm;
	im[1] /= norm;
	re[2] /= norm;
	im[2] /= norm;

	// Orthogonalize row 1 and 2
	rtmp  = re[0]*re[3] + im[0]*im[3];
	itmp  = re[0]*im[3] - im[0]*re[3];
	rtmp += re[1]*re[4] + im[1]*im[4];
	itmp += re[1]*im[4] - im[1]*re[4];
	rtmp += re[2]*re[5] + im[2]*im[5];
	itmp += re[2]*im[5] - im[2]*re[5];

	re[3] -= rtmp*re[0] - itmp*im[0];
	im[3] -= rtmp*im[0] + itmp*re[0];
	re[4] -= rtmp*re[1] - itmp*im[1];
	im[4] -= rtmp*im[1] + itmp*re[1];
	re[5] -= rtmp*re[2] - itmp*im[2];
	im[5] -= rtmp*im[2] + itmp*re[2];

	// Normalize the second row
	norm  = re[3]*re[3] + im[3]*im[3];
	norm += re[4]*re[4] + im[4]*im[4];
	norm += re[5]*re[5] + im[5]*im[5];
	norm  = sqrt(norm);

	re[3] /= norm;
	im[3] /= norm;
	re[4] /= norm;
	im[4] /= norm;
	re[5] /= norm;
	im[5] /= norm;

	// Construct third row as the conjugate of the cross product between row 1 and 2
	re[6] =  re[1]*re[5] - im[1]*im[5] - re[2]*re[4] + im[2]*im[4];
	im[6] = -re[1]*im[5] - im[1]*re[5] + re[2]*im[4] + im[2]*re[4];
	re[7] = -re[0]*re[5] + im[0]*im[5] + re[2]*re[3] - im[2]*im[3];
	im[7] =  re[0]*im[5] + im[0]*re[5] - re[2]*im[3] - im[2]*re[3];
	re[8] =  re[0]*re[4] - im[0]*im[4] - re[1]*re[3] + im[1]*im[3];
	im[8] = -re[0]*im[4] - im[0]*re[4] + re[1]*im[3] + im[1]*re[3];

	#else
	double norm, rtmp, itmp, dtmp;
	suNg copy;

	// Use Gram-Schmidt to orthogonalize all rows
	for(int i = 0; i < NC; i++)
	{
		norm = 0;
		for(int k = 0; k < NC; k++)
		{
			norm += re[i*NC+k]*re[i*NC+k] + im[i*NC+k]*im[i*NC+k];
		}
		norm = sqrt(norm);

		if(norm == 0)
		{
			return;
		}

		for(int k = 0; k < NC; k++)
		{
			re[i*NC+k] /= norm;
			im[i*NC+k] /= norm;
		}

		for(int j = i+1; j < NC; j++)
		{
			rtmp = 0;
			itmp = 0;

			for(int k = 0; k < NC; k++)
			{
				rtmp += re[i*NC+k]*re[j*NC+k] + im[i*NC+k]*im[j*NC+k];
				itmp += re[i*NC+k]*im[j*NC+k] - im[i*NC+k]*re[j*NC+k];
			}

			for(int k = 0; k < NC; k++)
			{
				re[j*NC+k] -= re[i*NC+k]*rtmp - im[i*NC+k]*itmp;
				im[j*NC+k] -= im[i*NC+k]*rtmp + re[i*NC+k]*itmp;
			}
		}
	}

	// Use LU decomposition to calculate the determinant
	copy = *this;

	for(int i = 0; i < NC; i++)
	{
		for(int j = 0; j <= i; j++)
		{
			rtmp = re[i*NC+j];
			itmp = im[i*NC+j];

			for(int k = 0; k < j; k++)
			{
				rtmp -= re[i*NC+k]*re[k*NC+j] - im[i*NC+k]*im[k*NC+j];
				itmp -= im[i*NC+k]*re[k*NC+j] + re[i*NC+k]*im[k*NC+j];
			}

			re[i*NC+j] = rtmp;
			im[i*NC+j] = itmp;
		}

		for(int j = i+1; j < NC; j++)
		{
			rtmp = re[i*NC+j];
			itmp = im[i*NC+j];

			for(int k = 0; k < i; k++)
			{
				rtmp -= re[i*NC+k]*re[k*NC+j] - im[i*NC+k]*im[k*NC+j];
				itmp -= im[i*NC+k]*re[k*NC+j] + re[i*NC+k]*im[k*NC+j];
			}

			norm = re[i*NC+i]*re[i*NC+i] + im[i*NC+i]*im[i*NC+i];
			re[i*NC+j] = (re[i*NC+i]*rtmp + im[i*NC+i]*itmp) / norm;
			im[i*NC+j] = (re[i*NC+i]*itmp - im[i*NC+i]*rtmp) / norm;
		}
	}

	rtmp = 1;
	itmp = 0;

	for(int i = 0; i < NC; i++)
	{
		dtmp = rtmp;
		rtmp = re[i*NC+i]*dtmp - im[i*NC+i]*itmp;
		itmp = re[i*NC+i]*itmp + im[i*NC+i]*dtmp;
	}

	// Multiply last row by a phase to ensure unit determinant
	*this = copy;

	dtmp = atan2(itmp, rtmp);
	rtmp = cos(dtmp);
	itmp = -sin(dtmp);

	for(int i = 0, j = NC-1; i < NC; i++)
	{
		dtmp = re[j*NC+i];
		re[j*NC+i] = dtmp*rtmp - im[j*NC+i]*itmp;
		im[j*NC+i] = dtmp*itmp + im[j*NC+i]*rtmp;
	}
	#endif
}

void suNg::unit()
{
	for(int i = 0; i < NC*NC; i++)
	{
		re[i] = 0;
		im[i] = 0;
	}

	for(int i = 0; i < NC; i++)
	{
		re[i*NC+i] = 1;
	}
}

void suNg::random()
{
	suNa alg;
	alg.random();
	*this = exp(alg);
}

double suNg::sqnorm()
{
	double ret = 0;
	for(int i = 0; i < NC*NC; i++)
	{
		ret += re[i]*re[i] + im[i]*im[i];
	}
	return ret;
}

int suNg::read(FILE *fp)
{
	int sz = 2*NC*NC*sizeof(double);
	for(int i = 0; i < NC*NC; i++)
	{
		sz -= fread(&re[i], sizeof(double), 1, fp);
		sz -= fread(&im[i], sizeof(double), 1, fp);
	}
	return sz;
}

int suNg::write(FILE *fp)
{
	int sz = 2*NC*NC*sizeof(double);
	for(int i = 0; i < NC*NC; i++)
	{
		sz -= fwrite(&re[i], sizeof(double), 1, fp);
		sz -= fwrite(&im[i], sizeof(double), 1, fp);
	}
	return sz;
}

suNg dagger(suNg g)
{
	suNg ret;
	for(int i = 0; i < NC; i++)
	{
		for(int j = 0; j < NC; j++)
		{
			ret.re[i*NC+j] =  g.re[j*NC+i];
			ret.im[i*NC+j] = -g.im[j*NC+i];
		}
	}
	return ret;
}

suNg transpose(suNg g)
{
	suNg ret;
	for(int i = 0; i < NC; i++)
	{
		for(int j = 0; j < NC; j++)
		{
			ret.re[i*NC+j] = g.re[j*NC+i];
			ret.im[i*NC+j] = g.im[j*NC+i];
		}
	}
	return ret;
}

double trace(suNg g)
{
	double ret = 0;
	for(int i = 0; i < NC; i++)
	{
		ret += g.re[i*NC+i];
	}
	return ret;
}

double trace(suNg &lhs, suNg &rhs)
{
	double re = 0;
	for(int i = 0; i < NC; i++)
	{
		for(int k = 0; k < NC; k++)
		{
			re += lhs.re[i*NC+k]*rhs.re[k*NC+i] - lhs.im[i*NC+k]*rhs.im[k*NC+i];
		}
	}
	return re;
}

complex ctrace(suNg g)
{
	complex ret;
	for(int i = 0; i < NC; i++)
	{
		ret.re += g.re[i*NC+i];
		ret.im += g.im[i*NC+i];
	}
	return ret;
}

suNg exp(suNg &X)
{
	suNg ret, Xk;
	double error;
	int k = 1;

	error = 1;
	Xk.unit();
	ret.unit();

	while(error > 1e-24)
	{
		Xk = (1./k)*Xk*X;
		ret = ret + Xk;
		error = Xk.sqnorm();
		k++;
	}

	return ret;
}
