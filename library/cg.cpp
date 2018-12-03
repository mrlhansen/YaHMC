#include <inverter.h>
#include <logger.h>
#include <linalg.h>

#define MAX 50

static SpinorField r, p, qqp;
static SpinorField k[MAX];
static int init = 0;

int cg_mshift(inv_par &par, SpinorField *dptr, SpinorField &sptr)
{
	int it = 0;
	double c1[MAX], c2[MAX], c3[MAX], cflag[MAX];
	double alpha, beta, beta_old;
	double atmp, btmp, norm;
	int notconverged;

	// Check operator
	if(par.mvm == 0)
	{
		lprintf("INVERTER", CRITICAL, "cg_mshift: no MVM operator specified");
	}

	// Allocate spinor fields
	if(init == 0)
	{
		spinor_allocate(r);
		spinor_allocate(p);
		spinor_allocate(qqp);
		init = 1;
	}

	// Initial guess must be zero for multi-shift inversion
	spinor_copy(r, sptr);
	spinor_copy(p, sptr);

	for(int n = 0; n < par.num; n++)
	{
		spinor_zero(dptr[n]);
		spinor_allocate(k[n]);
		spinor_copy(k[n], sptr);
	}

	// Initial value for the coefficients
	for(int n = 0; n < par.num; n++)
	{
		c1[n] = 1;
		c2[n] = 1;
		cflag[n] = 0;
	}

	beta = 1;
	alpha = 0;
	norm = spinor_sqnorm(sptr);

	// Perform the inversion
	while(++it)
	{
		par.mvm(par.mass, qqp, p);

		atmp = spinor_sqnorm(r);
		btmp = spinor_product_re(p, qqp);

		beta_old = beta;
		beta = -(atmp / btmp);

		for(int n = 0; n < par.num; n++)
		{
			// Continue if converged
			if(cflag[n]) continue;

			// Zeta coefficient
			c3[n] = c2[n] * c1[n] * beta_old;
			c3[n] /= (beta * alpha * (c1[n] - c2[n])) + (c1[n] * beta_old * (1 - par.shift[n] * beta));

			// Beta coefficient
			btmp = beta * (c3[n] / c2[n]);

			// Update spinor
			spinor_mulr_sub_assign(dptr[n], btmp, k[n]);
		}

		spinor_mulr_add_assign(r, beta, qqp);
		btmp = spinor_sqnorm(r);
		alpha = (btmp / atmp);
		notconverged = 0;

		for(int n = 0; n < par.num; n++)
		{
			// Continue if converged
			if(cflag[n]) continue;

			// Alpha coefficient
			atmp = alpha * (c3[n] * c3[n]) / (c2[n] * c2[n]);

			// Update spinor
			spinor_mulr_assign(k[n], atmp);
			spinor_mulr_add_assign(k[n], c3[n], r);

			if(btmp*c3[n]*c3[n] < norm*par.prec)
			{
				cflag[n] = 1;
			}
			else
			{
				notconverged++;
			}

			c1[n] = c2[n];
			c2[n] = c3[n];
		}

		if(notconverged == 0)
		{
			break;
		}

		spinor_mulr_assign(p, alpha);
		spinor_add_assign(p, r);
	}

	// Print log info
	if(par.log)
	{
		lprintf("INVERTER", NOTICE, "cg_mshift: iterations = %d, precision = %1.4e", it, btmp/norm);
	}

	return it;
}

int cg_inverter(inv_par &par, SpinorField &dptr, SpinorField &sptr)
{
	int it = 0;
	double alpha, beta, norm;
	double atmp, btmp;

	// Check operator
	if(par.mvm == 0)
	{
		lprintf("INVERTER", CRITICAL, "cg: no MVM operator specified");
	}

	// Allocate spinor fields
	if(init == 0)
	{
		spinor_allocate(r);
		spinor_allocate(p);
		spinor_allocate(qqp);
		init = 1;
	}

	// Calculate initial residue
	spinor_copy(p, dptr);
	par.mvm(par.mass, qqp, p);

	spinor_copy(r, sptr);
	spinor_sub_assign(r, qqp);
	spinor_copy(p, r);

	// Norm of input vector
	norm = spinor_sqnorm(sptr);
	btmp = spinor_sqnorm(r);

	// Check initial residue
	if(btmp < norm*par.prec)
	{
		it = -1;
	}

	// Perform the inversion
	while(++it)
	{
		par.mvm(par.mass, qqp, p);

		atmp = spinor_sqnorm(r);
		btmp = spinor_product_re(p, qqp);
		alpha = (atmp / btmp);

		spinor_mulr_add_assign(dptr, alpha, p);
		spinor_mulr_sub_assign(r, alpha, qqp);

		btmp = spinor_sqnorm(r);
		beta = (btmp / atmp);

		if(btmp < norm*par.prec)
		{
			break;
		}

		spinor_mulr_assign(p, beta);
		spinor_add_assign(p, r);
	}

	// Print log info
	if(par.log)
	{
		lprintf("INVERTER", NOTICE, "cg: iterations = %d, precision = %1.4e", it, btmp/norm);
	}

	return it;
}
