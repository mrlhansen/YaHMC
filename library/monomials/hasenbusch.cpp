#include <monomials.h>
#include <linalg.h>
#include <dirac.h>
#include <inverter.h>
#include <logger.h>
#include <clover.h>

MonomialHB::MonomialHB(double m, double s, double p)
{
	omass = m;
	smass = m + s;
	precision = p;

	spinor_allocate(pf);
	spinor_allocate(atmp);
	spinor_allocate(btmp);
	spinor_allocate(ctmp);

#ifdef UPDATE_EO
	dm = (4.0+smass)*(4.0+smass) - (4.0+omass)*(4.0+omass);
#else
	dm = s;
#endif

	lprintf("MONOMIAL", INFO, "Hasenbusch: mass = %1.6f, shift = %1.4f, prec = %1.1e", m, s, p);
}

void MonomialHB::reset()
{
	inv_par par;
	par.prec = precision;
	par.mass = smass;
	par.log = 1;
	par.mvm = Hphi;

	spinor_random(atmp);
	la = spinor_sqnorm(atmp);
	Hphi(omass, btmp, atmp);
	auto_inverter(par, pf, btmp);
}

double MonomialHB::action()
{
	inv_par par;
	double value;

	par.prec = precision;
	par.mass = omass;
	par.log = 1;
	par.mvm = Hphi;

	if(la > 0)
	{
		value = la;
		la = -1;
	}
	else
	{
		Hphi(smass, atmp, pf);
		auto_inverter(par, btmp, atmp);
		value = spinor_product_re(btmp, btmp);
	}

	return value;
}

void MonomialHB::update(double dt)
{
#ifdef CLOVER_TERM_EO
	inv_par par;
	par.prec = precision;
	par.mass = omass;
	par.log = 1;
	par.mvm = Hphi_sq;

	Hphi(smass, atmp, pf);
	auto_inverter(par, btmp, atmp);

	wilson_fermion_force_q(dt, omass, 1, btmp);
	wilson_fermion_force(dt, smass, -1, pf, btmp);
#else
	inv_par par;
	par.prec = precision;
	par.mass = omass;
	par.log = 1;
	par.mvm = Dphi;

	auto_inverter(par, atmp, pf);
	spinor_mulr(btmp, dm, atmp);
	spinor_add_assign(btmp, pf);
	spinor_g5(btmp);
	auto_inverter(par, ctmp, btmp);

	wilson_fermion_force(dt, omass, dm, atmp, ctmp);
#endif
}
