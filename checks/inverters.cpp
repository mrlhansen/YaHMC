#include <linalg.h>
#include <dirac.h>
#include <inverter.h>
#include <logger.h>
#include <cmath>

static double prec = 1e-28;
static double mass = -0.2;

static int check_bicgstab()
{
	SpinorField u, Mu, v;
	double r1, norm;
	inv_par par;
	int it;

	// Allocate memory
	spinor_allocate(u);
	spinor_allocate(v);
	spinor_allocate(Mu);

	// Random spinor fields
	spinor_random(u);
	norm = spinor_sqnorm(u);

	// Apply dirac operator
	par.prec = prec;
	par.mass = mass;
	par.log = 0;
	par.mvm = Dphi;
	it = bicgstab_inverter(par, Mu, u);

	// Check result
	Dphi(mass, v, Mu);
	spinor_sub_assign(v,u);
	r1 = spinor_sqnorm(v) / norm;

	// Evaluate result
	if(r1 < prec)
	{
		lprintf("CHECK", INFO, "BiCGstab inverter: PASSED [diff = %1.6e, iter = %d]", r1, it);
		return 0;
	}
	else
	{
		lprintf("CHECK", INFO, "BiCGstab inverter: FAILED [diff = %1.6e, iter = %d]", r1, it);
		return 1;
	}
}

static int check_minres()
{
	SpinorField u, Mu, v;
	double r1, norm;
	inv_par par;
	int it;

	// Allocate memory
	spinor_allocate(u);
	spinor_allocate(v);
	spinor_allocate(Mu);

	// Random spinor fields
	spinor_random(u);
	norm = spinor_sqnorm(u);

	// Apply dirac operator
	par.prec = prec;
	par.mass = mass;
	par.log = 0;
	par.mvm = Hphi;
	it = minres_inverter(par, Mu, u);

	// Check result
	Hphi(mass, v, Mu);
	spinor_sub_assign(v,u);
	r1 = spinor_sqnorm(v) / norm;

	// Evaluate result
	if(r1 < prec)
	{
		lprintf("CHECK", INFO, "MINRES inverter: PASSED [diff = %1.6e, iter = %d]", r1, it);
		return 0;
	}
	else
	{
		lprintf("CHECK", INFO, "MINRES inverter: FAILED [diff = %1.6e, iter = %d]", r1, it);
		return 1;
	}
}

static int check_cg()
{
	SpinorField u, Mu, v;
	double r1, norm;
	inv_par par;
	int it;

	// Allocate memory
	spinor_allocate(u);
	spinor_allocate(v);
	spinor_allocate(Mu);

	// Random spinor fields
	spinor_random(u);
	norm = spinor_sqnorm(u);

	// Apply dirac operator
	par.prec = prec;
	par.mass = mass;
	par.log = 0;
	par.mvm = Hphi_sq;
	it = cg_inverter(par, Mu, u);

	// Check result
	Hphi_sq(mass, v, Mu);
	spinor_sub_assign(v,u);
	r1 = spinor_sqnorm(v) / norm;

	// Evaluate result
	if(r1 < prec)
	{
		lprintf("CHECK", INFO, "CG inverter: PASSED [diff = %1.6e, iter = %d]", r1, it);
		return 0;
	}
	else
	{
		lprintf("CHECK", INFO, "CG inverter: FAILED [diff = %1.6e, iter = %d]", r1, it);
		return 1;
	}
}

static int check_cg_mshift()
{
	SpinorField u, Mu[4], v;
	double r1, norm;
	double s[4];
	inv_par par;
	int ret = 0;
	int it;

	// Shifts
	s[0] = 0.50;
	s[1] = 0.10;
	s[2] = 0.05;
	s[3] = 1.05;

	// Allocate memory
	spinor_allocate(u);
	spinor_allocate(v);
	spinor_allocate(Mu[0]);
	spinor_allocate(Mu[1]);
	spinor_allocate(Mu[2]);
	spinor_allocate(Mu[3]);

	// Random spinor fields
	spinor_random(u);
	norm = spinor_sqnorm(u);

	// Apply dirac operator
	par.prec = prec;
	par.mass = mass;
	par.log = 0;
	par.num = 4;
	par.shift = s;
	par.mvm = Hphi_sq;
	it = cg_mshift(par, Mu, u);

	// Check result
	for(int i = 0; i < par.num; i++)
	{
		Hphi_sq(mass, v, Mu[i]);
		spinor_mulr_add_assign(v, s[i], Mu[i]);
		spinor_sub_assign(v,u);
		r1 = spinor_sqnorm(v) / norm;

		// Evaluate result
		if(r1 < prec)
		{
			lprintf("CHECK", INFO, "CG_mshift inverter (shift %d): PASSED [diff = %1.6e, iter = %d]", i+1, r1, it);
		}
		else
		{
			lprintf("CHECK", INFO, "CG_mshift inverter (shift %d): FAILED [diff = %1.6e, iter = %d]", i+1, r1, it);
			ret = 1;
		}
	}

	// Return result
	return ret;
}

int check_inverters()
{
	int ret = 0;
	ret += check_bicgstab();
	ret += check_minres();
	ret += check_cg();
	ret += check_cg_mshift();
	return ret;
}
