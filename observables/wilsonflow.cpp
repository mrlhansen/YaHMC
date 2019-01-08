#include <wilsonflow.h>
#include <wilsonloop.h>
#include <global.h>
#include <staples.h>
#include <bcs.h>
#include <cmath>

void wf_zeta(double alpha)
{
	suNg stmp;
	suNa atmp;

	#pragma omp parallel for private(stmp,atmp)
	sites_for(id)
	{
		for(int mu = 0; mu < 4; mu++)
		{
			stmp = link(id,mu) * staples_wilson(id,mu);
			atmp = algebra_project(stmp, alpha);
			momentum(id,mu) -= atmp;
		}
	}
}

void wf_update(double val)
{
	suNg tmp;

	#pragma omp parallel for private(tmp)
	sites_for(id)
	{
		for(int mu = 0; mu < 4; mu++)
		{
			tmp = exp(momentum(id,mu));
			link(id,mu) = tmp * link(id,mu);
			momentum(id,mu) *= val;
		}
	}

	apply_bcs_on_gauge_field();
	mp_transfer_links();
}

void wf_iterate(double epsilon)
{
	wf_zeta(epsilon/4.0);
	wf_update(-17.0/9.0);

	wf_zeta(8.0*epsilon/9.0);
	wf_update(-1.0);

	wf_zeta(3.0*epsilon/4.0);
	wf_update(0);
}

double wf_e_plaq()
{
	double sum = 0;

	#pragma omp parallel for reduction(+:sum)
	sites_for(id)
	{
		for(int mu = 1; mu < 4; mu++)
		{
			for(int nu = 0; nu < mu; nu++)
			{
				sum += NC - wilson_loop(id, mu, nu, 1, 1);
			}
		}
	}

	mp_global_sum(&sum, 1);
	return 2*sum/global.vol4;
}

double wf_e_sym()
{
	double sum, norm;
	suNg stmp;
	suNa atmp;

	sum = 0;
	norm = sqrt(2.0*TF);

	#pragma omp parallel for private(stmp,atmp) reduction(+:sum)
	sites_for(id)
	{
		for(int mu = 1; mu < 4; mu++)
		{
			for(int nu = 0; nu < mu; nu++)
			{
				stmp = clover_fund(id, mu, nu);
				atmp = algebra_project(stmp, norm);
				sum += atmp.sqnorm();
			}
		}
	}

	mp_global_sum(&sum, 1);
	return (2*sum)/(16*4*global.vol4);
}

double wf_charge()
{
	double tc, norm;
	suNg v1, v2;
	suNa a1, a2;

	tc = 0;
	norm = sqrt(2.0*TF);

	#pragma omp parallel for private(v1,v2,a1,a2) reduction(+:tc)
	sites_for(id)
	{
		v1 = clover_fund(id, 0, 3);
		v2 = clover_fund(id, 1, 2);
		a1 = algebra_project(v1, norm);
		a2 = algebra_project(v2, norm);
		for(int i = 0; i < NG; i++)
		{
			tc += a1.av[i] * a2.av[i];
		}

		v1 = clover_fund(id, 0, 2);
		v2 = clover_fund(id, 3, 1);
		a1 = algebra_project(v1, norm);
		a2 = algebra_project(v2, norm);
		for(int i = 0; i < NG; i++)
		{
			tc += a1.av[i] * a2.av[i];
		}

		v1 = clover_fund(id, 0, 1);
		v2 = clover_fund(id, 2, 3);
		a1 = algebra_project(v1, norm);
		a2 = algebra_project(v2, norm);
		for(int i = 0; i < NG; i++)
		{
			tc += a1.av[i] * a2.av[i];
		}
	}

	mp_global_sum(&tc, 1);
	return (2*tc)/(16*16*M_PI*M_PI);
}

void wf_init()
{
	suNa zero;
	sites_for(id)
	{
		for(int mu = 0; mu < 4; mu++)
		{
			momentum(id,mu) = zero;
		}
	}
}
