#include <global.h>
#include <linalg.h>
#include <cstring>
#include <bcs.h>
#include <timing.h>

#define check_spinors(s1,s2) \
	if(s1.sites != s2.sites) lprintf("LINALG", CRITICAL, "Spinors do not match");

void spinor_allocate(SpinorField &sf)
{
	int memsz = (sizeof(Spinor)*outer.vol4)/1024;
	if(sf.allocated == false)
	{
		try
		{
			sf.ptr = new Spinor[outer.vol4];
			sf.allocated = true;
			spinor_type_default(sf);
			lprintf("LINALG", DEBUG, "Allocated %lu kb of memory", memsz);
		}
		catch(bad_alloc&)
		{
			lprintf("LINALG", CRITICAL, "Failed to allocate %lu kb of memory", memsz);
		}
	}
}

void spinor_zero(SpinorField &sf)
{
	timing_start(tm_linalg);

	#pragma omp parallel for
	sites_for(id)
	{
		sf[id] = Spinor();
	}

	timing_end(tm_linalg);
}

void spinor_random(SpinorField &sf)
{
	timing_start(tm_linalg);

	spinor_for(id,sf)
	{
		sf[id].randomize();
	}

	apply_bcs_on_spinor_field(sf);
	timing_end(tm_linalg);
}

double spinor_sqnorm(SpinorField &sf)
{
	double res = 0;
	timing_start(tm_linalg);

	#pragma omp parallel for reduction(+:res)
	spinor_for(id,sf)
	{
		res += sf[id].sqnorm();
	}

	timing_end(tm_linalg);
	mp_global_sum(&res, 1);
	return res;
}

void spinor_copy(SpinorField &lhs, SpinorField &rhs)
{
	timing_start(tm_linalg);
	memcpy(lhs.ptr, rhs.ptr, sizeof(Spinor)*outer.vol4);
	timing_end(tm_linalg);
}

void spinor_g5(SpinorField &sf)
{
	timing_start(tm_linalg);

	#pragma omp parallel for
	spinor_for(id,sf)
	{
		sf[id].gamma5();
	}

	timing_end(tm_linalg);
}

void spinor_add_assign(SpinorField &lhs, SpinorField &rhs)
{
	check_spinors(lhs, rhs);
	timing_start(tm_linalg);

	#pragma omp parallel for
	spinor_for(id,lhs)
	{
		lhs[id] += rhs[id];
	}

	timing_end(tm_linalg);
}

void spinor_sub_assign(SpinorField &lhs, SpinorField &rhs)
{
	check_spinors(lhs, rhs);
	timing_start(tm_linalg);

	#pragma omp parallel for
	spinor_for(id,lhs)
	{
		lhs[id] -= rhs[id];
	}

	timing_end(tm_linalg);
}

void spinor_mulr(SpinorField &lhs, double v, SpinorField &rhs)
{
	check_spinors(lhs, rhs);
	timing_start(tm_linalg);

	#pragma omp parallel for
	spinor_for(id,lhs)
	{
		lhs[id] = v*rhs[id];
	}

	timing_end(tm_linalg);
}

void spinor_mulc(SpinorField &lhs, complex v, SpinorField &rhs)
{
	check_spinors(lhs, rhs);
	timing_start(tm_linalg);

	#pragma omp parallel for
	spinor_for(id,lhs)
	{
		lhs[id] = v*rhs[id];
	}

	timing_end(tm_linalg);
}

void spinor_mulr_assign(SpinorField &lhs, double v)
{
	timing_start(tm_linalg);

	#pragma omp parallel for
	spinor_for(id,lhs)
	{
		lhs[id] *= v;
	}

	timing_end(tm_linalg);
}

void spinor_mulc_assign(SpinorField &lhs, complex v)
{
	timing_start(tm_linalg);

	#pragma omp parallel for
	spinor_for(id,lhs)
	{
		lhs[id] *= v;
	}

	timing_end(tm_linalg);
}

void spinor_mulr_add_assign(SpinorField &lhs, double v, SpinorField &rhs)
{
	check_spinors(lhs, rhs);
	timing_start(tm_linalg);

	#pragma omp parallel for
	spinor_for(id,lhs)
	{
		lhs[id] += v*rhs[id];
	}

	timing_end(tm_linalg);
}

void spinor_mulc_add_assign(SpinorField &lhs, complex v, SpinorField &rhs)
{
	check_spinors(lhs, rhs);
	timing_start(tm_linalg);

	#pragma omp parallel for
	spinor_for(id,lhs)
	{
		lhs[id] += v*rhs[id];
	}

	timing_end(tm_linalg);
}

void spinor_mulr_sub_assign(SpinorField &lhs, double v, SpinorField &rhs)
{
	check_spinors(lhs, rhs);
	timing_start(tm_linalg);

	#pragma omp parallel for
	spinor_for(id,lhs)
	{
		lhs[id] -= v*rhs[id];
	}

	timing_end(tm_linalg);
}

void spinor_mulc_sub_assign(SpinorField &lhs, complex v, SpinorField &rhs)
{
	check_spinors(lhs, rhs);
	timing_start(tm_linalg);

	#pragma omp parallel for
	spinor_for(id,lhs)
	{
		lhs[id] -= v*rhs[id];
	}

	timing_end(tm_linalg);
}

complex spinor_product(SpinorField &sf_a, SpinorField &sf_b)
{
	complex sum;
	double re = 0;
	double im = 0;

	check_spinors(sf_a, sf_b);
	timing_start(tm_linalg);

	#pragma omp parallel for private(sum) reduction(+:re,im)
	spinor_for(id,sf_a)
	{
		sum = inner_product(sf_a[id], sf_b[id]);
		re += sum.re;
		im += sum.im;
	}

	sum.re = re;
	sum.im = im;

	timing_end(tm_linalg);
	mp_global_sum(&sum, 2);
	return sum;
}

double spinor_product_re(SpinorField &sf_a, SpinorField &sf_b)
{
	double sum = 0;

	check_spinors(sf_a, sf_b);
	timing_start(tm_linalg);

	#pragma omp parallel for reduction(+:sum)
	spinor_for(id,sf_a)
	{
		sum += inner_product_re(sf_a[id], sf_b[id]);
	}

	timing_end(tm_linalg);
	mp_global_sum(&sum, 1);
	return sum;
}

// Spinor select
void spinor_select_even(SpinorField &sf)
{
	sf.offset = 0;
	sf.sites = (local.vol4 / 2);
}

void spinor_select_odd(SpinorField &sf)
{
	sf.offset = (local.vol4 / 2);
	sf.sites = (local.vol4 / 2);
}

void spinor_select_full(SpinorField &sf)
{
	sf.offset = 0;
	sf.sites = local.vol4;
}

void spinor_select_default(SpinorField &sf)
{
	sf.offset = sf.default_offset;
	sf.sites = sf.default_sites;
}

// Spinor type
void spinor_type_even(SpinorField &sf)
{
	sf.default_offset = 0;
	sf.default_sites = (local.vol4 / 2);
	sf.offset = sf.default_offset;
	sf.sites = sf.default_sites;
}

void spinor_type_odd(SpinorField &sf)
{
	sf.default_offset = (local.vol4 / 2);
	sf.default_sites = (local.vol4 / 2);
	sf.offset = sf.default_offset;
	sf.sites = sf.default_sites;
}

void spinor_type_full(SpinorField &sf)
{
	sf.default_offset = 0;
	sf.default_sites = local.vol4;
	sf.offset = sf.default_offset;
	sf.sites = sf.default_sites;
}

void spinor_type_default(SpinorField &sf)
{
#ifdef UPDATE_EO
	spinor_type_even(sf);
#else
	spinor_type_full(sf);
#endif
}
