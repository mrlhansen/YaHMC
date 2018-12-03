#include <repr.h>
#include <clover.h>
#include <global.h>
#include <bcs.h>
#include <cmath>

#define FL_re(idx) repr_link[id].re[idx]
#define FL_im(idx) repr_link[id].im[idx]
#define GL_re(i,j) fund_link[id].re[(i)*NC+j]
#define GL_im(i,j) fund_link[id].im[(i)*NC+j]
#define SQRT2 1.41421356237309504880168872421

suNg basis[REPR_DIM];
suNg iTfund[NG];
suNf iTrepr[NG];

void repr_init()
{
	suNa p;
	complex ctmp;
	suNg tmp;

	lprintf("REPR", INFO, "Representation: %s SU(%d) [dim = %d]", REPR_NAME, NC, REPR_DIM);

	#ifdef REPR_SYMMETRIC
	for(int i = 0, n = 0; i < NC; i++)
	{
		for(int j = 0; j <= i; j++, n++)
		{
			if(i == j)
			{
				basis[n].re[i*NC+i] = 1;
			}
			else
			{
				basis[n].re[i*NC+j] = 1.0/SQRT2;
				basis[n].re[j*NC+i] = 1.0/SQRT2;
			}
		}
	}
	#endif

	#ifdef REPR_ANTISYMMETRIC
	for(int i = 0, n = 0; i < NC; i++)
	{
		for(int j = 0; j < i; j++, n++)
		{
			basis[n].re[j*NC+i] = 1.0/SQRT2;
			basis[n].re[i*NC+j] = -1.0/SQRT2;
		}
	}
	#endif

	#ifdef REPR_ADJOINT
	ctmp = complex(0, -1.0/SQRT2);
	for(int n = 0; n < NG; n++)
	{
		p = suNa(n);
		basis[n] = ctmp*suNg(p);
	}
	#endif

	// Generators for the fundamental representation
	for(int n = 0; n < NG; n++)
	{
		p = suNa(n);
		iTfund[n] = suNg(p);
	}

	// Generators for the fermionic representation
	for(int k = 0; k < NG; k++)
	{
		for(int i = 0; i < REPR_DIM; i++)
		{
			for(int j = 0; j < REPR_DIM; j++)
			{
				#ifdef REPR_FUNDAMENTAL
				iTrepr[k].re[i*REPR_DIM+j] = iTfund[k].re[i*REPR_DIM+j];
				iTrepr[k].im[i*REPR_DIM+j] = iTfund[k].im[i*REPR_DIM+j];
				#endif

				#ifdef REPR_ADJOINT
				tmp = basis[j] * basis[i] * iTfund[k] + basis[i] * basis[j] * dagger(iTfund[k]);
				ctmp = ctrace(tmp);
				iTrepr[k].re[i*REPR_DIM+j] = ctmp.re;
				iTrepr[k].im[i*REPR_DIM+j] = ctmp.im;
				#endif

				#ifdef REPR_SYMMETRIC
				tmp = basis[j] * basis[i] * iTfund[k] + basis[i] * basis[j] * transpose(iTfund[k]);
				ctmp = ctrace(tmp);
				iTrepr[k].re[i*REPR_DIM+j] = ctmp.re;
				iTrepr[k].im[i*REPR_DIM+j] = ctmp.im;
				#endif

				#ifdef REPR_ANTISYMMETRIC
				tmp = basis[j] * transpose(basis[i]) * iTfund[k] + transpose(basis[i]) * basis[j] * transpose(iTfund[k]);
				ctmp = ctrace(tmp);
				iTrepr[k].re[i*REPR_DIM+j] = ctmp.re;
				iTrepr[k].im[i*REPR_DIM+j] = ctmp.im;
				#endif
			}
		}
	}
}

void represent_links(suNf *repr_link, suNg *fund_link, int count)
{
	for(int id = 0; id < count; id++)
	{
		#ifdef REPR_FUNDAMENTAL
		repr_link[id] = *(suNf*)&fund_link[id];
		#endif

		#ifdef REPR_ANTISYMMETRIC
		int n = 0;
		for(int i = 1; i < NC; i++)
		for(int j = 0; j < i; j++)
		for(int k = 1; k < NC; k++)
		for(int l = 0; l < k; l++)
		{
			FL_re(n)  = GL_re(i,k) * GL_re(j,l) - GL_im(i,k) * GL_im(j,l);
			FL_im(n)  = GL_re(i,k) * GL_im(j,l) + GL_im(i,k) * GL_re(j,l);
			FL_re(n) -= GL_re(i,l) * GL_re(j,k) - GL_im(i,l) * GL_im(j,k);
			FL_im(n) -= GL_re(i,l) * GL_im(j,k) + GL_im(i,l) * GL_re(j,k);
			n++;
		}
		#endif

		#ifdef REPR_SYMMETRIC
		int n = 0;
		for(int i = 0; i < NC; i++)
		for(int j = 0; j <= i; j++)
		for(int k = 0; k < NC; k++)
		for(int l = 0; l <= k; l++)
		{
			if(i == j && k == l)
			{
				FL_re(n) = GL_re(i,k) * GL_re(j,l) - GL_im(i,k) * GL_im(j,l);
				FL_im(n) = GL_re(i,k) * GL_im(j,l) + GL_im(i,k) * GL_re(j,l);
			}
			else if(i == j || k == l)
			{
				FL_re(n) = GL_re(i,k) * GL_re(j,l) - GL_im(i,k) * GL_im(j,l);
				FL_re(n) *= SQRT2;
				FL_im(n) = GL_re(i,k) * GL_im(j,l) + GL_im(i,k) * GL_re(j,l);
				FL_im(n) *= SQRT2;
			}
			else
			{
				FL_re(n)  = GL_re(i,k) * GL_re(j,l) - GL_im(i,k) * GL_im(j,l);
				FL_im(n)  = GL_re(i,k) * GL_im(j,l) + GL_im(i,k) * GL_re(j,l);
				FL_re(n) += GL_re(i,l) * GL_re(j,k) - GL_im(i,l) * GL_im(j,k);
				FL_im(n) += GL_re(i,l) * GL_im(j,k) + GL_im(i,l) * GL_re(j,k);
			}
			n++;
		}
		#endif

		#ifdef REPR_ADJOINT
		int n = 0;
		int m = 0;

		for(int i = 0; i < NC; i++)
		for(int j = 0; j < NC; j++)
		{
			if(i > j)
			{
				m = 0;
				for(int a = 0; a < NC; a++)
				for(int b = 0; b < NC; b++)
				{
					if(a > b)
					{
						FL_re(n*REPR_DIM+m) = GL_re(j,b)*GL_re(i,a) + GL_re(j,a)*GL_re(i,b) + GL_im(j,b)*GL_im(i,a) + GL_im(j,a)*GL_im(i,b);
						m++;
					}
					else if(a < b)
					{
						FL_re(n*REPR_DIM+m) = GL_re(i,b)*GL_im(j,a) + GL_re(j,b)*GL_im(i,a) - GL_re(i,a)*GL_im(j,b) - GL_re(j,a)*GL_im(i,b);
						m++;
					}
					else if(a > 0)
					{
						FL_re(n*REPR_DIM+m) = 0;
						for(int k = 0; k <= a; k++)
						{
							FL_re(n*REPR_DIM+m) += SQRT2 * basis[m].re[k*NC+k] * (GL_re(j,k)*GL_re(i,k) + GL_im(j,k)*GL_im(i,k));
						}
						m++;
					}
				}
				n++;
			}
			else if(i < j)
			{
				m = 0;
				for(int a = 0; a < NC; a++)
				for(int b = 0; b < NC; b++)
				{
					if(a > b)
					{
						FL_re(n*REPR_DIM+m) = GL_re(i,b)*GL_im(j,a) - GL_re(j,b)*GL_im(i,a) + GL_re(i,a)*GL_im(j,b) - GL_re(j,a)*GL_im(i,b);
						m++;
					}
					else if(a < b)
					{
						FL_re(n*REPR_DIM+m) = GL_re(j,b)*GL_re(i,a) - GL_re(j,a)*GL_re(i,b) + GL_im(j,b)*GL_im(i,a) - GL_im(j,a)*GL_im(i,b);
						m++;
					}
					else if(a > 0)
					{
						FL_re(n*REPR_DIM+m) = 0;
						for(int k = 0; k <= a; k++)
						{
							FL_re(n*REPR_DIM+m) += SQRT2 * basis[m].re[k*NC+k] * (GL_re(i,k)*GL_im(j,k) - GL_re(j,k)*GL_im(i,k));
						}
						m++;
					}
				}
				n++;
			}
			else if(i > 0)
			{
				m = 0;
				for(int a = 0; a < NC; a++)
				for(int b = 0; b < NC; b++)
				{
					if(a > b)
					{
						FL_re(n*REPR_DIM+m) = 0;
						for(int k = 0; k <= i; k++)
						{
							FL_re(n*REPR_DIM+m) += SQRT2 * basis[n].re[k*NC+k] * (GL_re(k,b)*GL_re(k,a) + GL_im(k,b)*GL_im(k,a));
						}
						m++;
					}
					else if(a < b)
					{
						FL_re(n*REPR_DIM+m) = 0;
						for(int k = 0; k <= i; k++)
						{
							FL_re(n*REPR_DIM+m) += SQRT2 * basis[n].re[k*NC+k] * (GL_re(k,b)*GL_im(k,a) - GL_im(k,b)*GL_re(k,a));
						}
						m++;
					}
					else if(a > 0)
					{
						FL_re(n*REPR_DIM+m) = 0;
						for(int x = 0; x <= i; x++)
						for(int y = 0; y <= a; y++)
						{
							FL_re(n*REPR_DIM+m) += basis[n].re[x*NC+x] * basis[m].re[y*NC+y] * (GL_re(x,y)*GL_re(x,y) + GL_im(x,y)*GL_im(x,y));
						}
						m++;
					}
				}
				n++;
			}
		}
		#endif
	}
}

void represent_gauge_field()
{
	mp_transfer_links();

	#pragma omp parallel for
	extended_sites_for(id)
	{
		represent_links(&fermion_link(id,0), &link(id,0), 4);
	}

	// Calculate clover term
	compute_clover_term();

	// Boundary conditions
	apply_bcs_on_represented_gauge_field();
}
