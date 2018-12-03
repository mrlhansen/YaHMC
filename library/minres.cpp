#include <inverter.h>
#include <logger.h>
#include <linalg.h>
#include <dirac.h>
#include <cmath>

static SpinorField v0, v1;
static SpinorField w0, w1;
static SpinorField tmp;
static int init = 0;

int minres_inverter(inv_par &par, SpinorField &dptr, SpinorField &sptr)
{
	int it = 0;
	double beta0, beta1, eta;
	double g0 = 1, g1 = 1;
	double s0 = 0, s1 = 0;
	double alpha, delta;
	double rho1, rho2, rho3;
	double norm, residue;
	Spinor *ptr;

	// Check operator
	if(par.mvm == 0)
	{
		lprintf("INVERTER", CRITICAL, "minres: no MVM operator specified");
	}

	// Allocate spinor fields
	if(init == 0)
	{
		spinor_allocate(v0);
		spinor_allocate(v1);
		spinor_allocate(w0);
		spinor_allocate(w1);
		spinor_allocate(tmp);
		init = 1;
	}

	// Initial residue
	par.mvm(par.mass, tmp, dptr);
	spinor_copy(v1, sptr);
	spinor_sub_assign(v1, tmp);

	// Zero spinors
	spinor_zero(v0);
	spinor_zero(w0);
	spinor_zero(w1);

	// Initial residue norm
	norm = spinor_sqnorm(sptr);
	residue = spinor_sqnorm(v1);
	beta0 = sqrt(residue);
	eta = beta0;

	// Check initial residue
	if(residue < norm*par.prec)
	{
		it = -1;
	}

	// Perform the inversion
	while(++it)
	{
		spinor_mulr_assign(v1, 1.0/beta0);
		par.mvm(par.mass, tmp, v1);
		alpha = spinor_product_re(v1, tmp);

		spinor_mulr_sub_assign(tmp, alpha, v1);
		spinor_mulr_sub_assign(tmp, beta0, v0);

		ptr = v0.ptr;
		v0.ptr = v1.ptr;
		v1.ptr = tmp.ptr;
		tmp.ptr = ptr;

		beta1 = spinor_sqnorm(v1);
		beta1 = sqrt(beta1);

		delta = g1*alpha - g0*s1*beta0;
		rho1 = sqrt(delta*delta + beta1*beta1);
		rho2 = s1*alpha + g0*g1*beta0;
		rho3 = s0*beta0;

		g0 = g1;
		g1 = delta/rho1;
		s0 = s1;
		s1 = beta1/rho1;

		spinor_mulr(tmp, 1.0/rho1, v0);
		spinor_mulr_sub_assign(tmp, rho3/rho1, w0);
		spinor_mulr_sub_assign(tmp, rho2/rho1, w1);

		ptr = w0.ptr;
		w0.ptr = w1.ptr;
		w1.ptr = tmp.ptr;
		tmp.ptr = ptr;

		spinor_mulr_add_assign(dptr, g1*eta, w1);
		eta = -s1*eta;
		beta0 = beta1;

		residue = eta*eta;
		if(residue < norm*par.prec)
		{
			break;
		}
	}

	// Print log info
	if(par.log)
	{
		lprintf("INVERTER", NOTICE, "minres: iterations = %d, precision = %1.4e", it, residue/norm);
	}

	return it;
}
