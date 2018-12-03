#include <global.h>
#include <cmath>
#include <linalg.h>
#include <spectrum.h>
#include <cmath>
#include <cstring>

static SpinorField source;
static SpinorField prop[2*REPR_DIM];
static int init = 0;

static void sf_mpcac_init()
{
	if(init == 0)
	{
		spinor_allocate(source);
		for(int n = 0; n < 2*REPR_DIM; n++)
		{
			spinor_allocate(prop[n]);
		}
		init = 1;
	}
}

static void sf_mpcac_print(const char *ch, double *v)
{
	string str;
	char buf[32];

	for(int i = 0; i < global.dim_t-3; i++)
	{
		sprintf(buf, "%1.8e ", v[i]);
		str.append(buf);
	}

	lprintf("SF", INFO, "%s: %s", ch, str.c_str());
}

void sf_mpcac_measure(double mass)
{
	double prec, f;
	Spinor s1, s2;
	int t, bk;

	double fp[global.dim_t];
	double fa[global.dim_t];
	double gp[global.dim_t];
	double ga[global.dim_t];

	memset(fp, 0, sizeof(fp));
	memset(fa, 0, sizeof(fa));
	memset(gp, 0, sizeof(gp));
	memset(ga, 0, sizeof(ga));

	sf_mpcac_init();
	f = 2./sqrt(global.vol3);
	prec = 1e-20;

	// Lower boundary
	for(int n = 0; n < 2*REPR_DIM; n++)
	{
		spinor_zero(source);
		sites_for(id)
		{
			if(global_time(id) == 1)
			{
				bk = bk_index(id,0);
				source[id].re[n%REPR_DIM][n/REPR_DIM] = f;
				source[id] = dagger(fermion_link(bk,0))*source[id];
				source[id].p_plus();
			}
		}
		calc_propagator(mass, prec, prop[n], source);
	}

	sites_for(id)
	{
		t = global_time(id);
		if(t == 0) continue;
		t = t-1;

		for(int n = 0; n < 2*REPR_DIM; n++)
		{
			s1 = prop[n].at(id);
			s2 = s1;
			fp[t] += 0.5*inner_product_re(s1,s2);
			s2.gamma0();
			fa[t] -= 0.5*inner_product_re(s1,s2);
		}
	}

	// Upper boundary
	for(int n = 0; n < 2*REPR_DIM; n++)
	{
		spinor_zero(source);
		sites_for(id)
		{
			if(global_time(id) == global.dim_t-3)
			{
				source[id].re[n%REPR_DIM][n/REPR_DIM] = f;
				source[id] = fermion_link(id,0)*source[id];
				source[id].p_minus();
			}
		}
		calc_propagator(mass, prec, prop[n], source);
	}

	sites_for(id)
	{
		t = global_time(id);
		if(t == 0) continue;
		t = t-1;

		for(int n = 0; n < 2*REPR_DIM; n++)
		{
			Spinor s1 = prop[n].at(id);
			Spinor s2 = s1;
			gp[t] += 0.5*inner_product_re(s1,s2);
			s2.gamma0();
			ga[t] += 0.5*inner_product_re(s1,s2);
		}
	}

	// Print correlators
	sf_mpcac_print("fp", fp);
	sf_mpcac_print("fa", fa);
	sf_mpcac_print("gp", gp);
	sf_mpcac_print("ga", ga);
}
