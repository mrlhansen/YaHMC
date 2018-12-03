#include <dirac.h>
#include <global.h>
#include <linalg.h>
#include <simd.h>
#include <clover.h>
#include <bcs.h>
#include <timing.h>

static int mvm_counter = 0;

int mvm_get()
{
#ifdef UPDATE_EO
	return (mvm_counter >> 1);
#else
	return mvm_counter;
#endif
}

void mvm_reset()
{
	mvm_counter = 0;
}

static void clover_kernel(double mass, SpinorField &dst, SpinorField &src, int part, bool assign)
{
	Spinor stmp;
	mass = 4.0+mass;

	// Timing
	timing_start(tm_dirac);

	// Select spinor block
	switch(part)
	{
		case EVEN:
			spinor_select_even(src);
			break;
		case ODD:
			spinor_select_odd(src);
			break;
		case BOTH:
			spinor_select_full(src);
			break;
		default:
			break;
	}

	#pragma omp parallel for private(stmp)
	spinor_for(id,src)
	{
		double cf = 0;

#ifdef DIRAC_BOUNDARY_TERM
		if(global_time(id) == bcs.Ta || global_time(id) == bcs.Tb)
		{
			cf = (bcs.cf-1.0);
		}
#endif

		stmp = (mass+cf)*src[id];

		for(int i = 0; i < REPR_DIM; i++)
		{
			for(int j = 0; j < REPR_DIM; j++)
			{
				int ij = i*REPR_DIM+j;
				int ji = j*REPR_DIM+i;

				stmp.re[i][0] += clover_re(id,0,ij) * src[id].re[j][0] - clover_im(id,0,ij) * src[id].im[j][0];
				stmp.im[i][0] += clover_re(id,0,ij) * src[id].im[j][0] + clover_im(id,0,ij) * src[id].re[j][0];
				stmp.re[i][0] += clover_re(id,1,ij) * src[id].re[j][1] - clover_im(id,1,ij) * src[id].im[j][1];
				stmp.im[i][0] += clover_re(id,1,ij) * src[id].im[j][1] + clover_im(id,1,ij) * src[id].re[j][1];

				stmp.re[i][1] += clover_re(id,1,ji) * src[id].re[j][0] + clover_im(id,1,ji) * src[id].im[j][0];
				stmp.im[i][1] += clover_re(id,1,ji) * src[id].im[j][0] - clover_im(id,1,ji) * src[id].re[j][0];
				stmp.re[i][1] -= clover_re(id,0,ij) * src[id].re[j][1] - clover_im(id,0,ij) * src[id].im[j][1];
				stmp.im[i][1] -= clover_re(id,0,ij) * src[id].im[j][1] + clover_im(id,0,ij) * src[id].re[j][1];

				stmp.re[i][2] += clover_re(id,2,ij) * src[id].re[j][2] - clover_im(id,2,ij) * src[id].im[j][2];
				stmp.im[i][2] += clover_re(id,2,ij) * src[id].im[j][2] + clover_im(id,2,ij) * src[id].re[j][2];
				stmp.re[i][2] += clover_re(id,3,ij) * src[id].re[j][3] - clover_im(id,3,ij) * src[id].im[j][3];
				stmp.im[i][2] += clover_re(id,3,ij) * src[id].im[j][3] + clover_im(id,3,ij) * src[id].re[j][3];

				stmp.re[i][3] += clover_re(id,3,ji) * src[id].re[j][2] + clover_im(id,3,ji) * src[id].im[j][2];
				stmp.im[i][3] += clover_re(id,3,ji) * src[id].im[j][2] - clover_im(id,3,ji) * src[id].re[j][2];
				stmp.re[i][3] -= clover_re(id,2,ij) * src[id].re[j][3] - clover_im(id,2,ij) * src[id].im[j][3];
				stmp.im[i][3] -= clover_re(id,2,ij) * src[id].im[j][3] + clover_im(id,2,ij) * src[id].re[j][3];
			}
		}

		if(assign)
		{
			dst[id] += stmp;
		}
		else
		{
			dst[id] = stmp;
		}
	}

	// Select default spinor block
	spinor_select_default(src);

	// Timing
	timing_end(tm_dirac);
}

static void clover_inv_kernel(double mass, SpinorField &dst, SpinorField &src, int part, bool assign)
{
	const int N = 2*REPR_DIM;
	complex x[N];
	Spinor stmp;
	mass = 4.0+mass;

	// Timing
	timing_start(tm_dirac);

	// Update LDL decomposition
	compute_ldl_decomp(mass);

	// Select spinor block
	switch(part)
	{
		case EVEN:
			spinor_select_even(src);
			break;
		case ODD:
			spinor_select_odd(src);
			break;
		case BOTH:
			spinor_select_full(src);
			break;
		default:
			break;
	}

	#pragma omp parallel for private(x,stmp)
	spinor_for(id,src)
	{
		// Upper two spin components
		for(int i = 0; i < REPR_DIM; i++)
		{
			x[i].re = src[id].re[i][0];
			x[i].im = src[id].im[i][0];
			x[i+REPR_DIM].re = src[id].re[i][1];
			x[i+REPR_DIM].im = src[id].im[i][1];
		}

		for(int i = 1; i < N; i++)
		{
			for(int k = 0; k < i; k++)
			{
				x[i] -= upper_ldl(id,i,k) * x[k];
			}
		}

		for(int i = N-1; i >= 0; i--)
		{
			x[i] /= upper_ldl(id,i,i).re;
			for(int k = i+1; k < N; k++)
			{
				x[i] -= conj(upper_ldl(id,k,i)) * x[k];
			}
		}

		for(int i = 0; i < REPR_DIM; i++)
		{
			stmp.re[i][0] = x[i].re;
			stmp.im[i][0] = x[i].im;
			stmp.re[i][1] = x[i+REPR_DIM].re;
			stmp.im[i][1] = x[i+REPR_DIM].im;
		}

		// Lower two spin components
		for(int i = 0; i < REPR_DIM; i++)
		{
			x[i].re = src[id].re[i][2];
			x[i].im = src[id].im[i][2];
			x[i+REPR_DIM].re = src[id].re[i][3];
			x[i+REPR_DIM].im = src[id].im[i][3];
		}

		for(int i = 1; i < N; i++)
		{
			for(int k = 0; k < i; k++)
			{
				x[i] -= lower_ldl(id,i,k) * x[k];
			}
		}

		for(int i = N-1; i >= 0; i--)
		{
			x[i] /= lower_ldl(id,i,i).re;
			for(int k = i+1; k < N; k++)
			{
				x[i] -= conj(lower_ldl(id,k,i)) * x[k];
			}
		}

		for(int i = 0; i < REPR_DIM; i++)
		{
			stmp.re[i][2] = x[i].re;
			stmp.im[i][2] = x[i].im;
			stmp.re[i][3] = x[i+REPR_DIM].re;
			stmp.im[i][3] = x[i+REPR_DIM].im;
		}

		// Store result
		if(assign)
		{
			dst[id] += stmp;
		}
		else
		{
			dst[id] = stmp;
		}
	}

	// Select default spinor block
	spinor_select_default(src);

	// Timing
	timing_end(tm_dirac);
}

static void dirac_kernel(SpinorField &dst, SpinorField &src, int part, double sign, bool assign)
{
	int fw, bk;
	Spinor sum, stmp;
	suNf *link;

	mvm_counter++;
	sign = -0.5*sign;

	// Select spinor block
	switch(part)
	{
		case EVEN:
			mp_transfer_spinor_odd(src);
			spinor_select_even(src);
			break;
		case ODD:
			mp_transfer_spinor_even(src);
			spinor_select_odd(src);
			break;
		case BOTH:
			mp_transfer_spinor(src);
			spinor_select_full(src);
			break;
		default:
			break;
	}

	// Timing
	timing_start(tm_dirac);

	#pragma omp parallel for private(stmp,sum,link,fw,bk) schedule(guided)
	spinor_for(id,src)
	{
		sse_vector res_re;
		sse_vector res_im;
		sse_vector rhs_re[REPR_DIM];
		sse_vector rhs_im[REPR_DIM];
		sse_vector lhs_re;
		sse_vector lhs_im;
		sse_vector tmp;

		// Forwards 0
		fw = fw_index(id,0);
		link = &fermion_link(id,0);

		for(int i = 0; i < REPR_DIM; i++)
		{
			stmp.re[i][0] = src[fw].re[i][0] + src[fw].re[i][2];
			stmp.im[i][0] = src[fw].im[i][0] + src[fw].im[i][2];
			stmp.re[i][1] = src[fw].re[i][1] + src[fw].re[i][3];
			stmp.im[i][1] = src[fw].im[i][1] + src[fw].im[i][3];

			sse_load(rhs_re[i], stmp.re[i][0]);
			sse_load(rhs_im[i], stmp.im[i][0]);
		}

		for(int i = 0; i < REPR_DIM; i++)
		{
			sse_set_zero(res_re);
			sse_set_zero(res_im);

			for(int j = 0; j < REPR_DIM; j++)
			{
				sse_set_dbl(lhs_re, link->re[i*REPR_DIM+j]);
				sse_set_dbl(lhs_im, link->im[i*REPR_DIM+j]);

				sse_mul(tmp, lhs_re, rhs_re[j]);
				sse_add_assign(res_re, tmp);
				sse_mul(tmp, lhs_im, rhs_im[j]);
				sse_sub_assign(res_re, tmp);

				sse_mul(tmp, lhs_re, rhs_im[j]);
				sse_add_assign(res_im, tmp);
				sse_mul(tmp, lhs_im, rhs_re[j]);
				sse_add_assign(res_im, tmp);
			}

			sse_store(stmp.re[i][2], res_re);
			sse_store(stmp.im[i][2], res_im);
		}

		for(int i = 0; i < REPR_DIM; i++)
		{
			sum.re[i][0] = stmp.re[i][2];
			sum.im[i][0] = stmp.im[i][2];
			sum.re[i][1] = stmp.re[i][3];
			sum.im[i][1] = stmp.im[i][3];
			sum.re[i][2] = stmp.re[i][2];
			sum.im[i][2] = stmp.im[i][2];
			sum.re[i][3] = stmp.re[i][3];
			sum.im[i][3] = stmp.im[i][3];
		}

		// Backwards 0
		bk = bk_index(id,0);
		link = &fermion_link(bk,0);

		for(int i = 0; i < REPR_DIM; i++)
		{
			stmp.re[i][0] = src[bk].re[i][2] - src[bk].re[i][0];
			stmp.im[i][0] = src[bk].im[i][2] - src[bk].im[i][0];
			stmp.re[i][1] = src[bk].re[i][3] - src[bk].re[i][1];
			stmp.im[i][1] = src[bk].im[i][3] - src[bk].im[i][1];

			sse_load(rhs_re[i], stmp.re[i][0]);
			sse_load(rhs_im[i], stmp.im[i][0]);
		}

		for(int i = 0; i < REPR_DIM; i++)
		{
			sse_set_zero(res_re);
			sse_set_zero(res_im);

			for(int j = 0; j < REPR_DIM; j++)
			{
				sse_set_dbl(lhs_re, link->re[j*REPR_DIM+i]);
				sse_set_dbl(lhs_im, link->im[j*REPR_DIM+i]);

				sse_mul(tmp, lhs_re, rhs_re[j]);
				sse_add_assign(res_re, tmp);
				sse_mul(tmp, lhs_im, rhs_im[j]);
				sse_add_assign(res_re, tmp);

				sse_mul(tmp, lhs_re, rhs_im[j]);
				sse_add_assign(res_im, tmp);
				sse_mul(tmp, lhs_im, rhs_re[j]);
				sse_sub_assign(res_im, tmp);
			}

			sse_store(stmp.re[i][2], res_re);
			sse_store(stmp.im[i][2], res_im);
		}

		for(int i = 0; i < REPR_DIM; i++)
		{
			sum.re[i][0] -= stmp.re[i][2];
			sum.im[i][0] -= stmp.im[i][2];
			sum.re[i][1] -= stmp.re[i][3];
			sum.im[i][1] -= stmp.im[i][3];
			sum.re[i][2] += stmp.re[i][2];
			sum.im[i][2] += stmp.im[i][2];
			sum.re[i][3] += stmp.re[i][3];
			sum.im[i][3] += stmp.im[i][3];
		}

		// Forwards 1
		fw = fw_index(id,1);
		link = &fermion_link(id,1);

		for(int i = 0; i < REPR_DIM; i++)
		{
			stmp.re[i][0] = src[fw].re[i][2] + src[fw].im[i][1];
			stmp.im[i][0] = src[fw].im[i][2] - src[fw].re[i][1];
			stmp.re[i][1] = src[fw].re[i][3] + src[fw].im[i][0];
			stmp.im[i][1] = src[fw].im[i][3] - src[fw].re[i][0];

			sse_load(rhs_re[i], stmp.re[i][0]);
			sse_load(rhs_im[i], stmp.im[i][0]);
		}

		for(int i = 0; i < REPR_DIM; i++)
		{
			sse_set_zero(res_re);
			sse_set_zero(res_im);

			for(int j = 0; j < REPR_DIM; j++)
			{
				sse_set_dbl(lhs_re, link->re[i*REPR_DIM+j]);
				sse_set_dbl(lhs_im, link->im[i*REPR_DIM+j]);

				sse_mul(tmp, lhs_re, rhs_re[j]);
				sse_add_assign(res_re, tmp);
				sse_mul(tmp, lhs_im, rhs_im[j]);
				sse_sub_assign(res_re, tmp);

				sse_mul(tmp, lhs_re, rhs_im[j]);
				sse_add_assign(res_im, tmp);
				sse_mul(tmp, lhs_im, rhs_re[j]);
				sse_add_assign(res_im, tmp);
			}

			sse_store(stmp.re[i][2], res_re);
			sse_store(stmp.im[i][2], res_im);
		}

		for(int i = 0; i < REPR_DIM; i++)
		{
			sum.re[i][0] -= stmp.im[i][3];
			sum.im[i][0] += stmp.re[i][3];
			sum.re[i][1] -= stmp.im[i][2];
			sum.im[i][1] += stmp.re[i][2];
			sum.re[i][2] += stmp.re[i][2];
			sum.im[i][2] += stmp.im[i][2];
			sum.re[i][3] += stmp.re[i][3];
			sum.im[i][3] += stmp.im[i][3];
		}

		// Backwards 1
		bk = bk_index(id,1);
		link = &fermion_link(bk,1);

		for(int i = 0; i < REPR_DIM; i++)
		{
			stmp.re[i][0] = src[bk].re[i][2] - src[bk].im[i][1];
			stmp.im[i][0] = src[bk].im[i][2] + src[bk].re[i][1];
			stmp.re[i][1] = src[bk].re[i][3] - src[bk].im[i][0];
			stmp.im[i][1] = src[bk].im[i][3] + src[bk].re[i][0];

			sse_load(rhs_re[i], stmp.re[i][0]);
			sse_load(rhs_im[i], stmp.im[i][0]);
		}

		for(int i = 0; i < REPR_DIM; i++)
		{
			sse_set_zero(res_re);
			sse_set_zero(res_im);

			for(int j = 0; j < REPR_DIM; j++)
			{
				sse_set_dbl(lhs_re, link->re[j*REPR_DIM+i]);
				sse_set_dbl(lhs_im, link->im[j*REPR_DIM+i]);

				sse_mul(tmp, lhs_re, rhs_re[j]);
				sse_add_assign(res_re, tmp);
				sse_mul(tmp, lhs_im, rhs_im[j]);
				sse_add_assign(res_re, tmp);

				sse_mul(tmp, lhs_re, rhs_im[j]);
				sse_add_assign(res_im, tmp);
				sse_mul(tmp, lhs_im, rhs_re[j]);
				sse_sub_assign(res_im, tmp);
			}

			sse_store(stmp.re[i][2], res_re);
			sse_store(stmp.im[i][2], res_im);
		}

		for(int i = 0; i < REPR_DIM; i++)
		{
			sum.re[i][0] += stmp.im[i][3];
			sum.im[i][0] -= stmp.re[i][3];
			sum.re[i][1] += stmp.im[i][2];
			sum.im[i][1] -= stmp.re[i][2];
			sum.re[i][2] += stmp.re[i][2];
			sum.im[i][2] += stmp.im[i][2];
			sum.re[i][3] += stmp.re[i][3];
			sum.im[i][3] += stmp.im[i][3];
		}

		// Forwards 2
		fw = fw_index(id,2);
		link = &fermion_link(id,2);

		for(int i = 0; i < REPR_DIM; i++)
		{
			stmp.re[i][0] = src[fw].re[i][2] - src[fw].re[i][1];
			stmp.im[i][0] = src[fw].im[i][2] - src[fw].im[i][1];
			stmp.re[i][1] = src[fw].re[i][3] + src[fw].re[i][0];
			stmp.im[i][1] = src[fw].im[i][3] + src[fw].im[i][0];

			sse_load(rhs_re[i], stmp.re[i][0]);
			sse_load(rhs_im[i], stmp.im[i][0]);
		}

		for(int i = 0; i < REPR_DIM; i++)
		{
			sse_set_zero(res_re);
			sse_set_zero(res_im);

			for(int j = 0; j < REPR_DIM; j++)
			{
				sse_set_dbl(lhs_re, link->re[i*REPR_DIM+j]);
				sse_set_dbl(lhs_im, link->im[i*REPR_DIM+j]);

				sse_mul(tmp, lhs_re, rhs_re[j]);
				sse_add_assign(res_re, tmp);
				sse_mul(tmp, lhs_im, rhs_im[j]);
				sse_sub_assign(res_re, tmp);

				sse_mul(tmp, lhs_re, rhs_im[j]);
				sse_add_assign(res_im, tmp);
				sse_mul(tmp, lhs_im, rhs_re[j]);
				sse_add_assign(res_im, tmp);
			}

			sse_store(stmp.re[i][2], res_re);
			sse_store(stmp.im[i][2], res_im);
		}

		for(int i = 0; i < REPR_DIM; i++)
		{
			sum.re[i][0] += stmp.re[i][3];
			sum.im[i][0] += stmp.im[i][3];
			sum.re[i][1] -= stmp.re[i][2];
			sum.im[i][1] -= stmp.im[i][2];
			sum.re[i][2] += stmp.re[i][2];
			sum.im[i][2] += stmp.im[i][2];
			sum.re[i][3] += stmp.re[i][3];
			sum.im[i][3] += stmp.im[i][3];
		}

		// Backwards 2
		bk = bk_index(id,2);
		link = &fermion_link(bk,2);

		for(int i = 0; i < REPR_DIM; i++)
		{
			stmp.re[i][0] = src[bk].re[i][1] + src[bk].re[i][2];
			stmp.im[i][0] = src[bk].im[i][1] + src[bk].im[i][2];
			stmp.re[i][1] = src[bk].re[i][3] - src[bk].re[i][0];
			stmp.im[i][1] = src[bk].im[i][3] - src[bk].im[i][0];

			sse_load(rhs_re[i], stmp.re[i][0]);
			sse_load(rhs_im[i], stmp.im[i][0]);
		}

		for(int i = 0; i < REPR_DIM; i++)
		{
			sse_set_zero(res_re);
			sse_set_zero(res_im);

			for(int j = 0; j < REPR_DIM; j++)
			{
				sse_set_dbl(lhs_re, link->re[j*REPR_DIM+i]);
				sse_set_dbl(lhs_im, link->im[j*REPR_DIM+i]);

				sse_mul(tmp, lhs_re, rhs_re[j]);
				sse_add_assign(res_re, tmp);
				sse_mul(tmp, lhs_im, rhs_im[j]);
				sse_add_assign(res_re, tmp);

				sse_mul(tmp, lhs_re, rhs_im[j]);
				sse_add_assign(res_im, tmp);
				sse_mul(tmp, lhs_im, rhs_re[j]);
				sse_sub_assign(res_im, tmp);
			}

			sse_store(stmp.re[i][2], res_re);
			sse_store(stmp.im[i][2], res_im);
		}

		for(int i = 0; i < REPR_DIM; i++)
		{
			sum.re[i][0] -= stmp.re[i][3];
			sum.im[i][0] -= stmp.im[i][3];
			sum.re[i][1] += stmp.re[i][2];
			sum.im[i][1] += stmp.im[i][2];
			sum.re[i][2] += stmp.re[i][2];
			sum.im[i][2] += stmp.im[i][2];
			sum.re[i][3] += stmp.re[i][3];
			sum.im[i][3] += stmp.im[i][3];
		}

		// Forwards 3
		fw = fw_index(id,3);
		link = &fermion_link(id,3);

		for(int i = 0; i < REPR_DIM; i++)
		{
			stmp.re[i][0] = src[fw].re[i][2] + src[fw].im[i][0];
			stmp.im[i][0] = src[fw].im[i][2] - src[fw].re[i][0];
			stmp.re[i][1] = src[fw].re[i][3] - src[fw].im[i][1];
			stmp.im[i][1] = src[fw].im[i][3] + src[fw].re[i][1];

			sse_load(rhs_re[i], stmp.re[i][0]);
			sse_load(rhs_im[i], stmp.im[i][0]);
		}

		for(int i = 0; i < REPR_DIM; i++)
		{
			sse_set_zero(res_re);
			sse_set_zero(res_im);

			for(int j = 0; j < REPR_DIM; j++)
			{
				sse_set_dbl(lhs_re, link->re[i*REPR_DIM+j]);
				sse_set_dbl(lhs_im, link->im[i*REPR_DIM+j]);

				sse_mul(tmp, lhs_re, rhs_re[j]);
				sse_add_assign(res_re, tmp);
				sse_mul(tmp, lhs_im, rhs_im[j]);
				sse_sub_assign(res_re, tmp);

				sse_mul(tmp, lhs_re, rhs_im[j]);
				sse_add_assign(res_im, tmp);
				sse_mul(tmp, lhs_im, rhs_re[j]);
				sse_add_assign(res_im, tmp);
			}

			sse_store(stmp.re[i][2], res_re);
			sse_store(stmp.im[i][2], res_im);
		}

		for(int i = 0; i < REPR_DIM; i++)
		{
			sum.re[i][0] -= stmp.im[i][2];
			sum.im[i][0] += stmp.re[i][2];
			sum.re[i][1] += stmp.im[i][3];
			sum.im[i][1] -= stmp.re[i][3];
			sum.re[i][2] += stmp.re[i][2];
			sum.im[i][2] += stmp.im[i][2];
			sum.re[i][3] += stmp.re[i][3];
			sum.im[i][3] += stmp.im[i][3];
		}

		// Backwards 3
		bk = bk_index(id,3);
		link = &fermion_link(bk,3);

		for(int i = 0; i < REPR_DIM; i++)
		{
			stmp.re[i][0] = src[bk].re[i][2] - src[bk].im[i][0];
			stmp.im[i][0] = src[bk].re[i][0] + src[bk].im[i][2];
			stmp.re[i][1] = src[bk].re[i][3] + src[bk].im[i][1];
			stmp.im[i][1] = src[bk].im[i][3] - src[bk].re[i][1];

			sse_load(rhs_re[i], stmp.re[i][0]);
			sse_load(rhs_im[i], stmp.im[i][0]);
		}

		for(int i = 0; i < REPR_DIM; i++)
		{
			sse_set_zero(res_re);
			sse_set_zero(res_im);

			for(int j = 0; j < REPR_DIM; j++)
			{
				sse_set_dbl(lhs_re, link->re[j*REPR_DIM+i]);
				sse_set_dbl(lhs_im, link->im[j*REPR_DIM+i]);

				sse_mul(tmp, lhs_re, rhs_re[j]);
				sse_add_assign(res_re, tmp);
				sse_mul(tmp, lhs_im, rhs_im[j]);
				sse_add_assign(res_re, tmp);

				sse_mul(tmp, lhs_re, rhs_im[j]);
				sse_add_assign(res_im, tmp);
				sse_mul(tmp, lhs_im, rhs_re[j]);
				sse_sub_assign(res_im, tmp);
			}

			sse_store(stmp.re[i][2], res_re);
			sse_store(stmp.im[i][2], res_im);
		}

		for(int i = 0; i < REPR_DIM; i++)
		{
			sum.re[i][0] += stmp.im[i][2];
			sum.im[i][0] -= stmp.re[i][2];
			sum.re[i][1] -= stmp.im[i][3];
			sum.im[i][1] += stmp.re[i][3];
			sum.re[i][2] += stmp.re[i][2];
			sum.im[i][2] += stmp.im[i][2];
			sum.re[i][3] += stmp.re[i][3];
			sum.im[i][3] += stmp.im[i][3];
		}

		// Store result
		if(assign)
		{
			dst[id] += sign*sum;
		}
		else
		{
			dst[id] = sign*sum;
		}
	}

	// Timing
	timing_end(tm_dirac);

	// Boundary conditions
	apply_bcs_on_spinor_field(dst);

	// Default spinor
	spinor_select_default(src);
}

void Dphi_oe(SpinorField &dptr, SpinorField &sptr)
{
	// dptr_o = D_oe * sptr_e
	dirac_kernel(dptr, sptr, ODD, 1, false);
}

void Dphi_eo(SpinorField &dptr, SpinorField &sptr)
{
	// dptr_e = D_eo * sptr_o
	dirac_kernel(dptr, sptr, EVEN, 1, false);
}

void Dphi_ee(double mass, SpinorField &dptr, SpinorField &sptr)
{
	// dptr_e = D_ee * sptr_e
	clover_kernel(mass, dptr, sptr, EVEN, false);
}

void Dphi_oo(double mass, SpinorField &dptr, SpinorField &sptr)
{
	// dptr_o = D_oo * sptr_o
	clover_kernel(mass, dptr, sptr, ODD, false);
}

void Dphi_oo_inv(double mass, SpinorField &dptr, SpinorField &sptr)
{
	// dptr_o = D_oo^-1 * sptr_o
	clover_inv_kernel(mass, dptr, sptr, ODD, false);
}

#ifdef CLOVER_TERM

void Dphi(double mass, SpinorField &dptr, SpinorField &sptr)
{
#ifdef UPDATE_EO
	static SpinorField tmp;
	spinor_allocate(tmp);
	dirac_kernel(sptr, sptr, ODD, -1, false);
	clover_inv_kernel(mass, tmp, sptr, ODD, false);
	dirac_kernel(dptr, tmp, EVEN, 1, false);
	clover_kernel(mass, dptr, sptr, EVEN, true);
#else
	dirac_kernel(dptr, sptr, BOTH, 1, false);
	clover_kernel(mass, dptr, sptr, BOTH, true);
#endif
}

#else

void Dphi(double mass, SpinorField &dptr, SpinorField &sptr)
{
#ifdef UPDATE_EO
	dirac_kernel(sptr, sptr, ODD, -1, false);
	dirac_kernel(dptr, sptr, EVEN, 1, false);
	spinor_mulr_add_assign(dptr, (4.0+mass)*(4.0+mass), sptr);
#else
	dirac_kernel(dptr, sptr, BOTH, 1, false);
	spinor_mulr_add_assign(dptr, (4.0+mass), sptr);
#endif
}

#endif

void Hphi(double mass, SpinorField &dptr, SpinorField &sptr)
{
	Dphi(mass, dptr, sptr);
	spinor_g5(dptr);
}

void Hphi_sq(double mass, SpinorField &dptr, SpinorField &sptr)
{
	static SpinorField tmp;
	spinor_allocate(tmp);
	Hphi(mass, tmp, sptr);
	Hphi(mass, dptr, tmp);
}
