#include <mre.h>
#include <logger.h>
#include <linalg.h>
#include <solver.h>
#include <cmath>

// Auxillary spinor field
static SpinorField out;

// Variables used in the LU solver
static complex A[MRE_MAX*MRE_MAX];
static complex b[MRE_MAX];
static complex x[MRE_MAX];

void gram_schmidt(mre_par &par, int max)
{
	complex rij;
	double rii;

	for(int i = 0; i < max; i++)
	{
		rii = spinor_sqnorm(par.s[i]);
		rii = 1.0/sqrt(rii);
		spinor_mulr_assign(par.s[i], rii);

		for(int j = i+1; j < max; j++)
		{
			rij = spinor_product(par.s[i], par.s[j]);
			spinor_mulc_sub_assign(par.s[j], rij, par.s[i]);
		}
	}
}

void mre_store(mre_par &par, SpinorField &sptr)
{
	if(par.init == 0)
	{
		return;
	}

	int last = par.max - 1;
	Spinor *tmp = par.s[last].ptr;

	for(int i = last; i > 0; i--)
	{
		par.s[i].ptr = par.s[i-1].ptr;
	}

	par.s[0].ptr = tmp;
	spinor_copy(par.s[0], sptr);
	par.num++;
}

void mre_guess(mre_par &par, SpinorField &dptr, SpinorField &sptr)
{
	if(par.init == 0)
	{
		return;
	}

	int max = (par.num > par.max) ? par.max : par.num;
	spinor_zero(dptr);
	gram_schmidt(par, max);

	for(int i = 0; i < max; i++)
	{
		par.mvm(par.mass, out, par.s[i]);

		for(int j = 0; j < max; j++)
		{
			A[j*max+i] = spinor_product(out, par.s[j]);
		}

		b[i] = spinor_product(sptr, par.s[i]);
	}

	lu(max, A, x, b);

	for(int k = 0; k < max ; k++)
	{
		spinor_mulc_add_assign(dptr, x[k], par.s[k]);
	}
}

void mre_init(mre_par &par, int max, mvm_operator mvm, double mass)
{
	if(max == 0)
	{
		par.max = 0;
		par.mvm = 0;
		par.init = 0;
		par.mass = 0;
		return;
	}
	else
	{
		par.max = (max < MRE_MAX) ? max : MRE_MAX;
		par.mvm = mvm;
		par.init = 1;
		par.mass = mass;
	}

	// Allocate memory
	spinor_allocate(out);
	for(int n = 0; n < par.max; n++)
	{
		spinor_allocate(par.s[n]);
	}

	// Log info
	lprintf("MRE", NOTICE, "Enabled MRE with %d past solutions", par.max);
}

void mre_reset(mre_par &par)
{
	par.num = 0;
}
