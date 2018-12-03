#include <monomials.h>
#include <linalg.h>
#include <dirac.h>
#include <mre.h>
#include <inverter.h>
#include <logger.h>
#include <settings.h>
#include <clover.h>

MonomialHMC::MonomialHMC(double m, double s, double p, int n)
{
	omass = m;
	smass = m + s;
	precision = p;

	mre_init(mpar1, n, Dphi, smass);
	mre_init(mpar2, n, Dphi, smass);

	spinor_allocate(pf);
	spinor_allocate(atmp);
	spinor_allocate(btmp);

	lprintf("MONOMIAL", INFO, "HMC: mass = %1.6f, precision = %1.1e", smass, p);
}

void MonomialHMC::reset()
{
	spinor_random(atmp);
	la = spinor_sqnorm(atmp);
	Dphi(smass, pf, atmp);

#ifdef CLOVER_TERM_EO
	la -= 2*clover_logdet(omass);
#endif

	mre_reset(mpar1);
	mre_reset(mpar2);
}

double MonomialHMC::action()
{
	inv_par par;
	double value;

	par.prec = precision;
	par.mass = smass;
	par.log = 1;
	par.mvm = Dphi;

	if(la != 0)
	{
		value = la;
		la = 0;
	}
	else
	{
		auto_inverter(par, atmp, pf);
		value = spinor_product_re(atmp, atmp);
#ifdef CLOVER_TERM_EO
		value -= 2*clover_logdet(omass);
#endif
	}

	return value;
}

void MonomialHMC::update(double dt)
{
	inv_par par;
	par.prec = precision;
	par.mass = smass;
	par.log = 1;
	par.mvm = Dphi;

	mre_guess(mpar1, atmp, pf);
	if(auto_inverter(par, atmp, pf))
	{
		mre_store(mpar1, atmp);
	}

	spinor_g5(atmp);
	mre_guess(mpar2, btmp, atmp);
	if(auto_inverter(par, btmp, atmp))
	{
		mre_store(mpar2, btmp);
	}

	wilson_fermion_force_q(dt, smass, 1, btmp);

#ifdef CLOVER_TERM_EO
	clover_logdet_force(omass, 2);
#endif
}
