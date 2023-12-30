#include <clover.h>
#include <global.h>
#include <solver.h>
#include <bcs.h>
#include <memory.h>
#include <wilsonloop.h>
#include <timing.h>
#include <cmath>

#define N (2*REPR_DIM)

static double cl_sigma;
static double cl_csw;
static int init = 0;

suNf *cl_force;
suNf *cl_term;
ldl_t *cl_ldl;

double get_csw()
{
	return cl_csw;
}

static void _compute_ldl_decomp(int id)
{
	int o1, o2, o3, ndx;
	complex A[N*N];
	complex B[N*N];
	double cf = 0;

#ifdef DIRAC_BOUNDARY_TERM
	if(global_time(id) == bcs.Ta || global_time(id) == bcs.Tb)
	{
		cf = (bcs.cf-1.0);
	}
#endif

	// Construct matrices
	for(int i = 0; i < REPR_DIM; i++)
	{
		o1 = N * i;
		o2 = N * (i + REPR_DIM);
		o3 = o2 + REPR_DIM;

		for(int j = 0; j < REPR_DIM; j++)
		{
			ndx = j*REPR_DIM+i;
			A[o2+j].re =  clover_re(id,1,ndx);
			A[o2+j].im = -clover_im(id,1,ndx);
			B[o2+j].re =  clover_re(id,3,ndx);
			B[o2+j].im = -clover_im(id,3,ndx);

			if(i >= j)
			{
				ndx = i*REPR_DIM+j;
				A[o1+j].re =  clover_re(id,0,ndx);
				A[o1+j].im =  clover_im(id,0,ndx);
				A[o3+j].re = -clover_re(id,0,ndx);
				A[o3+j].im = -clover_im(id,0,ndx);
				B[o1+j].re =  clover_re(id,2,ndx);
				B[o1+j].im =  clover_im(id,2,ndx);
				B[o3+j].re = -clover_re(id,2,ndx);
				B[o3+j].im = -clover_im(id,2,ndx);

				if(i == j)
				{
					A[o1+j].re += cl_sigma+cf;
					A[o3+j].re += cl_sigma+cf;
					B[o1+j].re += cl_sigma+cf;
					B[o3+j].re += cl_sigma+cf;
				}
			}
		}
	}

	// LDL factorization
	ldl(N, A);
	ldl(N, B);

	// Store result
	for(int i = 0; i < N; i++)
	{
		for(int j = 0; j <= i; j++)
		{
			o2 = i*N+j;
			upper_ldl(id,i,j) = A[o2];
			lower_ldl(id,i,j) = B[o2];
		}
	}
}

static int _bcs_clover_term(int id)
{
	suNf u;
	int t, cond;

	t = global_time(id);
	cond = 0;

#ifdef BC_SF
	cond = (t == 0 || t == global.dim_t-1 || t == global.dim_t-2);
#endif

#ifdef BC_OPEN
	cond = (t == 0 || t == global.dim_t-1);
#endif

	if(cond)
	{
		clover_term(id,0) = u;
		clover_term(id,1) = u;
		clover_term(id,2) = u;
		clover_term(id,3) = u;
	}

	return cond;
}

static void _compute_clover_term(int id)
{
	suNf tmp[6];
	double csw;
	double atmp_re, atmp_im;
	double btmp_re, btmp_im;
	double ctmp_re, ctmp_im;
	double dtmp_re, dtmp_im;

	if(_bcs_clover_term(id))
	{
		return;
	}

	csw = -cl_csw/16.0;
	tmp[0] = clover_repr(id, 0, 1);
	tmp[1] = clover_repr(id, 0, 2);
	tmp[2] = clover_repr(id, 0, 3);
	tmp[3] = clover_repr(id, 1, 2);
	tmp[4] = clover_repr(id, 1, 3);
	tmp[5] = clover_repr(id, 2, 3);

	for(int i = 0; i < REPR_DIM; i++)
	{
		for(int j = 0; j < REPR_DIM; j++)
		{
			int ij = i*REPR_DIM+j;
			int ji = j*REPR_DIM+i;

			atmp_re = tmp[2].im[ij] + tmp[2].im[ji];
			atmp_im = tmp[2].re[ji] - tmp[2].re[ij];
			btmp_re = tmp[3].im[ij] + tmp[3].im[ji];
			btmp_im = tmp[3].re[ij] - tmp[3].re[ji];
			ctmp_re = tmp[0].im[ij] + tmp[0].im[ji] - tmp[1].re[ij] + tmp[1].re[ji];
			ctmp_im = tmp[0].re[ji] - tmp[0].re[ij] - tmp[1].im[ij] - tmp[1].im[ji];
			dtmp_re = tmp[4].re[ij] - tmp[4].re[ji] + tmp[5].im[ij] + tmp[5].im[ji];
			dtmp_im = tmp[4].im[ij] + tmp[4].im[ji] - tmp[5].re[ij] + tmp[5].re[ji];

			clover_re(id,0,ij) =  csw * (atmp_re - btmp_re);
			clover_im(id,0,ij) =  csw * (atmp_im + btmp_im);
			clover_re(id,1,ij) =  csw * (ctmp_re - dtmp_re);
			clover_im(id,1,ij) =  csw * (ctmp_im - dtmp_im);
			clover_re(id,2,ij) = -csw * (atmp_re + btmp_re);
			clover_im(id,2,ij) =  csw * (btmp_im - atmp_im);
			clover_re(id,3,ij) = -csw * (ctmp_re + dtmp_re);
			clover_im(id,3,ij) = -csw * (ctmp_im + dtmp_im);
		}
	}
}

static void _compute_clover_force(int id, double coeff)
{
	complex A[N][N];
	complex B[N][N];

	double a11_re, a11_im;
	double a12_re, a12_im;
	double a21_re, a21_im;
	double a22_re, a22_im;
	double a33_re, a33_im;
	double a34_re, a34_im;
	double a43_re, a43_im;
	double a44_re, a44_im;

	// Calculate inverse from LDL
	for(int n = 0; n < N; n++)
	{
		A[n][n].re = coeff;
		B[n][n].re = coeff;

		for(int i = n; i < N; i++)
		{
			for(int k = 0; k < i; k++)
			{
				A[i][n] -= upper_ldl(id,i,k) * A[k][n];
				B[i][n] -= lower_ldl(id,i,k) * B[k][n];
			}
		}

		for(int i = N-1; i >= n; i--)
		{
			A[i][n] /= upper_ldl(id,i,i).re;
			B[i][n] /= lower_ldl(id,i,i).re;
			for(int k = i+1; k < N; k++)
			{
				A[i][n] -= conj(upper_ldl(id,k,i)) * A[k][n];
				B[i][n] -= conj(lower_ldl(id,k,i)) * B[k][n];
			}
		}
	}

	// Construct force matrices
	for(int i = 0; i < REPR_DIM; i++)
	{
		for(int j = 0; j < REPR_DIM; j++)
		{
			int ij = i*REPR_DIM+j;
			a21_re =  A[i+REPR_DIM][j].re;
			a21_im =  A[i+REPR_DIM][j].im;
			a12_re =  A[j+REPR_DIM][i].re;
			a12_im = -A[j+REPR_DIM][i].im;
			a43_re =  B[i+REPR_DIM][j].re;
			a43_im =  B[i+REPR_DIM][j].im;
			a34_re =  B[j+REPR_DIM][i].re;
			a34_im = -B[j+REPR_DIM][i].im;

			if(i < j)
			{
				a11_re =  A[j][i].re;
				a11_im = -A[j][i].im;
				a22_re =  A[j+REPR_DIM][i+REPR_DIM].re;
				a22_im = -A[j+REPR_DIM][i+REPR_DIM].im;
				a33_re =  B[j][i].re;
				a33_im = -B[j][i].im;
				a44_re =  B[j+REPR_DIM][i+REPR_DIM].re;
				a44_im = -B[j+REPR_DIM][i+REPR_DIM].im;
			}
			else
			{
				a11_re = A[i][j].re;
				a11_im = A[i][j].im;
				a22_re = A[i+REPR_DIM][j+REPR_DIM].re;
				a22_im = A[i+REPR_DIM][j+REPR_DIM].im;
				a33_re = B[i][j].re;
				a33_im = B[i][j].im;
				a44_re = B[i+REPR_DIM][j+REPR_DIM].re;
				a44_im = B[i+REPR_DIM][j+REPR_DIM].im;
			}

			clover_force(id,0).re[ij] +=  a12_im + a21_im - a34_im - a43_im; // X_01
			clover_force(id,0).im[ij] += -a12_re - a21_re + a34_re + a43_re;
			clover_force(id,1).re[ij] +=  a12_re - a21_re + a43_re - a34_re; // X_02
			clover_force(id,1).im[ij] +=  a12_im - a21_im + a43_im - a34_im;
			clover_force(id,2).re[ij] +=  a22_im - a11_im + a44_im - a33_im; // X_12
			clover_force(id,2).im[ij] +=  a11_re - a22_re + a33_re - a44_re;
			clover_force(id,3).re[ij] +=  a11_im - a22_im + a44_im - a33_im; // X_03
			clover_force(id,3).im[ij] +=  a22_re - a11_re + a33_re - a44_re;
			clover_force(id,4).re[ij] +=  a12_re - a21_re + a34_re - a43_re; // X_13
			clover_force(id,4).im[ij] +=  a12_im - a21_im + a34_im - a43_im;
			clover_force(id,5).re[ij] += -a12_im - a21_im - a34_im - a43_im; // X_23
			clover_force(id,5).im[ij] +=  a12_re + a21_re + a34_re + a43_re;
		}
	}
}

void compute_clover_force(double mass, double coeff)
{
	if(init == 0)
	{
		return;
	}

	compute_ldl_decomp(4.0+mass);

	#pragma omp parallel for
	odd_sites_for(id)
	{
		_compute_clover_force(id, coeff);
	}
}

double clover_logdet(double mass)
{
	if(init == 0)
	{
		return 0;
	}

	double det = 0;
	compute_ldl_decomp(4.0+mass);

	#pragma omp parallel for reduction(+:det)
	odd_sites_for(id)
	{
		double prod = 1;
		for(int n = 0; n < N; n++)
		{
			prod *= upper_ldl(id,n,n).re;
			prod *= lower_ldl(id,n,n).re;
		}
		det += log(prod);
	}

	mp_global_sum(&det, 1);
	return det;
}

void compute_ldl_decomp(double sigma)
{
	if(init == 0)
	{
		return;
	}

	if(sigma == cl_sigma)
	{
		return;
	}
	else
	{
		cl_sigma = sigma;
	}

	#pragma omp parallel for
	odd_sites_for(id)
	{
		_compute_ldl_decomp(id);
	}
}

void compute_clover_term()
{
	if(init == 0)
	{
		return;
	}

	cl_sigma = 0xF00F;
	timing_start(tm_clover_term);

	#pragma omp parallel for
	sites_for(id)
	{
		_compute_clover_term(id);
	}

	timing_end(tm_clover_term);
}

void clover_init(double csw)
{
#ifdef CLOVER_TERM
	init = 1;
	cl_csw = csw;

	cl_force = suNf_allocate(6);
	cl_term = suNf_allocate(4);
	cl_ldl = ldl_allocate(1);

	lprintf("CLOVER", INFO, "Coefficient: csw = %1.4f", cl_csw);
#endif
}
