#include <spectrum.h>
#include <global.h>
#include <inverter.h>
#include <linalg.h>
#include <random.h>
#include <dirac.h>

static int init = 0;
static SpinorField tmp, eo_source;

int random_timeslice()
{
#if (defined(BC_OPEN) || defined(BC_SF))
	return 1;
#else
	int tau = rand_uniform()*global.dim_t;
	mp_broadcast(&tau);
	return tau;
#endif
}

void create_semwall_source(int n, int tau, SpinorField &source)
{
	// Manually construct the wall source
	if(n == 0)
	{
		sites_for(id)
		{
			source[id] = Spinor();

			if(global_time(id) != tau)
			{
				continue;
			}

			for(int i = 0; i < REPR_DIM; i++)
			{
				source[id].re[i][0] = rand_z2();
				source[id].im[i][0] = rand_z2();
			}
		}
	}
	else
	{
		sites_for(id)
		{
			for(int i = 0; i < REPR_DIM; i++)
			{
				source[id].re[i][n] = source[id].re[i][n-1];
				source[id].im[i][n] = source[id].im[i][n-1];
				source[id].re[i][n-1] = 0;
				source[id].im[i][n-1] = 0;
			}
		}
	}
}

void create_point_source(int n, int tau, SpinorField &source)
{
	int id, c, s;

	// n = s*REPR_DIM + c
	s = n/REPR_DIM;
	c = n%REPR_DIM;

	// Manually construct the point source
	id = extended_index_global(tau, 0, 0, 0);
	spinor_zero(source);
	if(id >= 0)
	{
		source[id].re[c][s] = 1;
	}
}

void calc_propagator_eo(double mass, double prec, SpinorField &dptr, SpinorField &sptr)
{
	if(init == 0)
	{
		spinor_allocate(tmp);
		spinor_allocate(eo_source);
		init = 1;
	}

	// Inverter settings
	inv_par par;
	par.prec = prec;
	par.mass = mass;
	par.log = 1;
	par.mvm = Dphi;

#ifdef CLOVER_TERM
	// eo_source = b_e - D_eo D_oo^-1 b_o
	Dphi_oo_inv(mass, tmp, sptr);
	Dphi_eo(tmp, tmp);
	spinor_mulr(eo_source, -1.0, tmp);
	spinor_add_assign(eo_source, sptr);

	// Call inverter
	auto_inverter(par, dptr, eo_source);

	// tmp = D_oe x_e
	Dphi_oe(tmp, dptr);

	// x_o = D_oo^-1(b_o - D_oe x_e)
	spinor_select_odd(sptr);
	spinor_select_odd(dptr);
	spinor_select_odd(tmp);

	spinor_mulr(dptr, 1.0, sptr);
	spinor_sub_assign(dptr, tmp);
	Dphi_oo_inv(mass, dptr, dptr);

	spinor_select_default(sptr);
	spinor_select_default(dptr);
	spinor_select_default(tmp);
#else
	// Adjust mass
	mass = 4.0 + mass;

	// eo_source = D_oo b_e - D_eo b_o
	Dphi_eo(tmp, sptr);
	spinor_mulr(eo_source, mass, sptr);
	spinor_sub_assign(eo_source, tmp);

	// Call inverter
	auto_inverter(par, dptr, eo_source);

	// tmp = D_oe x_e
	Dphi_oe(tmp, dptr);

	// x_o = D_oo^-1(b_o - D_oe x_e)
	spinor_select_odd(sptr);
	spinor_select_odd(dptr);
	spinor_select_odd(tmp);

	spinor_mulr(dptr, 1./mass, sptr);
	spinor_mulr_sub_assign(dptr, 1./mass, tmp);

	spinor_select_default(sptr);
	spinor_select_default(dptr);
	spinor_select_default(tmp);
#endif
}

void calc_propagator_std(double mass, double prec, SpinorField &dptr, SpinorField &sptr)
{
	inv_par par;

	// Inverter settings
	par.prec = prec;
	par.mass = mass;
	par.log = 1;
	par.mvm = Dphi;

	// Invert dirac operator
	auto_inverter(par, dptr, sptr);
}

void calc_propagator(double mass, double prec, SpinorField &dptr, SpinorField &sptr)
{
	#ifdef UPDATE_EO
	calc_propagator_eo(mass, prec, dptr, sptr);
	#else
	calc_propagator_std(mass, prec, dptr, sptr);
	#endif
}
