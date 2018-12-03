#include <monomials.h>
#include <evlimit.h>
#include <linalg.h>
#include <dirac.h>
#include <cstdlib>
#include <inverter.h>
#include <logger.h>
#include <clover.h>

MonomialRHMC::MonomialRHMC(double m, double rap, double inv)
{
	mass = m;
	inv_precision = inv;
	rap_precision = rap;

	n_fc = -1;
	d_fc = 2;
	n_pf = 1;
	d_pf = 4;

	spinor_allocate(pf);
	spinor_allocate(tmp);
	out = 0;
	order = 0;
	la = 0;

	lprintf("MONOMIAL", INFO, "RHMC: mass = %1.6f, rap prec = %1.1e, inv prec = %1.1e, frac = %d/%d", m, rap, inv, n_fc, d_fc);
}

void MonomialRHMC::reset()
{
	double evmin, evmax;

	evlimits(Hphi_sq, mass, &evmin, &evmax);
	evmin *= 0.95;
	evmax *= 1.05;

	rational_find(rpf, n_pf, d_pf, rap_precision, evmin, evmax);
	rational_find(rfc, n_fc, d_fc, rap_precision, evmin, evmax);

	spinor_random(tmp);
	la = spinor_sqnorm(tmp);
	rational_calc(rpf, Hphi_sq, mass, inv_precision, pf, tmp);

#ifdef CLOVER_TERM_EO
	la -= 1.0*clover_logdet(mass);
#endif

	if(rfc.order != order)
	{
		if(out)
		{
			delete[] out;
		}

		order = rfc.order;
		out = new SpinorField[order];
		for(int i = 0; i < order; i++)
		{
			spinor_allocate(out[i]);
		}
	}
}

double MonomialRHMC::action()
{
	double value;

	if(la != 0)
	{
		value = la;
		la = 0;
	}
	else
	{
		rational_calc(rfc, Hphi_sq, mass, inv_precision, tmp, pf);
		value = spinor_product_re(pf, tmp);
#ifdef CLOVER_TERM_EO
		value -= 1.0*clover_logdet(mass);
#endif
	}

	return value;
}

void MonomialRHMC::update(double dt)
{
	inv_par par;
	par.mass = mass;
	par.prec = inv_precision;
	par.num = rfc.order;
	par.shift = rfc.pole;
	par.log = 1;
	par.mvm = Hphi_sq;

	cg_mshift(par, out, pf);

	for(int n = 0; n < rfc.order; n++)
	{
		wilson_fermion_force_q(dt, mass, rfc.residue[n], out[n]);
	}

#ifdef CLOVER_TERM_EO
	clover_logdet_force(mass, 1.0);
#endif
}
