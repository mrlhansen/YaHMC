#include <global.h>
#include <linalg.h>
#include <dirac.h>
#include <logger.h>
#include <cmath>

static int check_hermitian()
{
	SpinorField u, v, Mu, Mv;
	complex r1, r2;
	double norm, r3;

	// Allocate memory
	spinor_allocate(u);
	spinor_allocate(v);
	spinor_allocate(Mu);
	spinor_allocate(Mv);

	// Random spinor fields
	spinor_random(u);
	spinor_random(v);

	// Apply dirac operator
	Hphi(-0.5, Mu, u);
	Hphi(-0.5, Mv, v);

	// Calculate norm
	norm = spinor_sqnorm(u);
	norm *= spinor_sqnorm(v);
	norm = sqrt(norm);

	// Calculate difference
	r1 = spinor_product(u, Mv);
	r2 = spinor_product(Mu, v);

	r1 = r1 - r2;
	r3 = cabs(r1);
	r3 /= norm;

	// Evaluate result
	if(r3 < 1e-14)
	{
		lprintf("CHECK", INFO, "Hermiticity of the Dirac operator: PASSED [diff = %1.6e]", r3);
		return 0;
	}
	else
	{
		lprintf("CHECK", INFO, "Hermiticity of the Dirac operator: FAILED [diff = %1.6e]", r3);
		return 1;
	}
}

#ifdef CLOVER_TERM
static int check_clover_inv()
{
	SpinorField u, v;
	double r1, norm;

	// Allocate memory
	spinor_allocate(u);
	spinor_allocate(v);
	spinor_type_odd(u);
	spinor_type_odd(v);

	// Random spinor field
	spinor_random(u);
	norm = spinor_sqnorm(u);

	// Apply dirac operator
	Dphi_oo(-0.5, v, u);
	Dphi_oo_inv(-0.5, v, v);

	// Calculate difference
	spinor_sub_assign(u, v);
	r1 = spinor_sqnorm(u);
	r1 = sqrt(r1/norm);

	// Evaluate result
	if(r1 < 1e-14)
	{
		lprintf("CHECK", INFO, "Inversion of the clover term: PASSED [diff = %1.6e]", r1);
		return 0;
	}
	else
	{
		lprintf("CHECK", INFO, "Inversion of the clover term: FAILED [diff = %1.6e]", r1);
		return 1;
	}
}
#endif

static double lcg_rand()
{
	static long max = (1L << 32);
	static long state = 1234;
	state = (1664525*state + 1013904223) % max;
	return (double)state/(double)max;
}

static void random_gtrans()
{
	suNa alg;
	int id;

	for(int t = 0; t < global.dim_t; t++)
	for(int x = 0; x < global.dim_x; x++)
	for(int y = 0; y < global.dim_y; y++)
	for(int z = 0; z < global.dim_z; z++)
	{
		id = extended_index_global(t, x, y, z);

		for(int a = 0; a < NG; a++)
		{
			alg.av[a] = lcg_rand();
		}

		if(id >= 0)
		{
			gfield_copy[id] = exp(alg);
		}
	}
}

static void transform_spinor(SpinorField &sptr)
{
	suNf tmp;
	spinor_for(id,sptr)
	{
		represent_links(&tmp, &gfield_copy[id], 1);
		sptr[id] = tmp*sptr[id];
	}
}

static void transform_gfield()
{
	sites_for(id)
	{
		for(int mu = 0; mu < 4; mu++)
		{
			int fw = fw_index(id,mu);
			link(id,mu) = gfield_copy[id]*link(id,mu)*dagger(gfield_copy[fw]);
		}
	}
	represent_gauge_field();
}

static int check_covariance()
{
	double r1, norm;
	SpinorField s0, s1, s2;

	// Allocate spinor fields
	spinor_allocate(s0);
	spinor_allocate(s1);
	spinor_allocate(s2);

	// Random fields
	spinor_random(s0);
	random_gtrans();
	norm = spinor_sqnorm(s0);

	// Dphi
	Dphi(-0.5, s1, s0);

	// Transform fields
	transform_spinor(s1);
	transform_spinor(s0);
	transform_gfield();

	// Dphi
	Dphi(-0.5, s2, s0);

	// Calculate difference
	spinor_sub_assign(s1, s2);
	r1 = spinor_sqnorm(s1);
	r1 = sqrt(r1/norm);

	// Evaluate result
	if(r1 < 1e-14)
	{
		lprintf("CHECK", INFO, "Covariance of the Dirac operator: PASSED [diff = %1.6e]", r1);
		return 0;
	}
	else
	{
		lprintf("CHECK", INFO, "Covariance of the Dirac operator: FAILED [diff = %1.6e]", r1);
		return 1;
	}
}

int check_dirac()
{
	int ret = 0;
	ret += check_hermitian();
	ret += check_covariance();
#ifdef CLOVER_TERM
	ret += check_clover_inv();
#endif
	return ret;
}
