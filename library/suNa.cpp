#include <global.h>
#include <suN.h>
#include <random.h>

suNa::suNa()
{
	for(int i = 0; i < NG; i++)
	{
		av[i] = 0;
	}
}

suNa::suNa(int n)
{
	for(int i = 0; i < NG; i++)
	{
		av[i] = 0;
	}
	if(n < NG)
	{
		av[n] = 1;
	}
}

suNa operator+(const suNa &lhs, const suNa &rhs)
{
	suNa ret(lhs);
	ret += rhs;
	return ret;
}

suNa& suNa::operator+=(const suNa &rhs)
{
	for(int i = 0; i < NG; i++)
	{
		av[i] += rhs.av[i];
	}
	return *this;
}

suNa operator-(const suNa &lhs, const suNa &rhs)
{
	suNa ret(lhs);
	ret -= rhs;
	return ret;
}

suNa& suNa::operator-=(const suNa &rhs)
{
	for(int i = 0; i < NG; i++)
	{
		av[i] -= rhs.av[i];
	}
	return *this;
}

suNa operator*(const double &lhs, const suNa &rhs)
{
	suNa ret(rhs);
	ret *= lhs;
	return ret;
}

suNa& suNa::operator*=(const double &rhs)
{
	for(int i = 0; i < NG; i++)
	{
		av[i] *= rhs;
	}
	return *this;
}

double suNa::sqnorm()
{
	double ret = 0;
	for(int i = 0; i < NG; i++)
	{
		ret += av[i]*av[i];
	}
	return ret;
}

void suNa::random()
{
	for(int i = 0; i < NG; i++)
	{
		av[i] = rand_gaussian();
	}
}

suNa algebra_project(suNg &m, double coeff)
{
	suNa ret;
	coeff /= -TF;
	for(int a = 0; a < NG; a++)
	{
		ret.av[a] = coeff*trace(iTfund[a], m);
	}
	return ret;
}

suNa algebra_project(suNf &m, double coeff)
{
	suNa ret;
	coeff /= -TF;
	for(int a = 0; a < NG; a++)
	{
		ret.av[a] = coeff*trace(iTrepr[a], m);
	}
	return ret;
}

suNg exp(suNa &m, double coeff)
{
	suNg X(m);
	suNg ret, Xk;
	double error;
	int k = 1;

	error = 1;
	Xk.unit();
	ret.unit();

	while(error > 1e-24)
	{
		Xk = (coeff/k)*Xk*X;
		ret = ret + Xk;
		error = Xk.sqnorm();
		k++;
	}

	ret.reunitarize();
	return ret;
}
