#include <inverter.h>
#include <logger.h>
#include <linalg.h>
#include <complex.h>

static SpinorField r, r0, v, p, s, t;
static int init = 0;

int bicgstab_inverter(inv_par &par, SpinorField &dptr, SpinorField &sptr)
{
	int it = 0;
	double norm, residue;
	complex alpha, beta;
	complex rho0, rho1;
	complex omega;

	// Check operator
	if(par.mvm == 0)
	{
		lprintf("INVERTER", CRITICAL, "bicgstab: no MVM operator specified");
	}

	// Allocate spinor fields
	if(init == 0)
	{
		spinor_allocate(r);
		spinor_allocate(r0);
		spinor_allocate(v);
		spinor_allocate(p);
		spinor_allocate(s);
		spinor_allocate(t);
		init = 1;
	}

	// Initial residue
	par.mvm(par.mass, p, dptr);
	spinor_copy(r, sptr);
	spinor_sub_assign(r, p);

	// Initialize spinors
	spinor_copy(r0, r);
	spinor_copy(p, r);

	// Initialize variables
	norm = spinor_sqnorm(sptr);
	rho0 = spinor_product(r, r0);

	// Check initial residue
	residue = rho0.re;
	if(residue < norm*par.prec)
	{
		it = -1;
	}

	// Perform the inversion
	while(++it)
	{
		par.mvm(par.mass, v, p);
		alpha = spinor_product(v, r0);
		alpha = rho0 / alpha;

		spinor_mulc_sub_assign(r, alpha, v);
		par.mvm(par.mass, t, r);

		omega = spinor_product(r, t);
		omega /= spinor_sqnorm(t);
		spinor_mulc_add_assign(dptr, alpha, p);
		spinor_mulc_add_assign(dptr, omega, r);

		spinor_mulc_sub_assign(r, omega, t);
		rho1 = spinor_product(r, r0);
		beta = (rho1 / rho0) * (alpha / omega);

		spinor_mulc_assign(p, beta);
		spinor_add_assign(p, r);
		spinor_mulc_sub_assign(p, beta*omega, v);

		residue = spinor_sqnorm(r);
		if(residue < norm*par.prec)
		{
			break;
		}

		rho0 = rho1;
	}

	// Each iteration requires two calls to the Dirac operator
	it *= 2;

	// Print log info
	if(par.log)
	{
		lprintf("INVERTER", NOTICE, "bicgstab: iterations = %d, precision = %1.4e", it, residue/norm);
	}

	return it;
}
