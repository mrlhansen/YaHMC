#ifndef CLOVER_H
#define CLOVER_H

#include <repr.h>
#include <complex.h>
#include <suN.h>

#if (defined(UPDATE_EO) && defined(CLOVER_TERM))
#define CLOVER_TERM_EO
#endif

#define clover_term(id,mu) \
	cl_term[4*(id)+mu]

#define clover_force(id,mu) \
	cl_force[6*(id)+mu]

#define clover_re(id,mu,ndx) \
	cl_term[4*(id)+mu].re[ndx]

#define clover_im(id,mu,ndx) \
	cl_term[4*(id)+mu].im[ndx]

#define upper_ldl(id,i,j) \
	cl_ldl[id].upper[(i)*(i+1)/2+j]

#define lower_ldl(id,i,j) \
	cl_ldl[id].lower[(i)*(i+1)/2+j]

typedef struct {
	complex upper[REPR_DIM*(2*REPR_DIM+1)];
	complex lower[REPR_DIM*(2*REPR_DIM+1)];
} ldl_t;

extern suNf *cl_force;
extern suNf *cl_term;
extern ldl_t *cl_ldl;

double get_csw();
void compute_clover_force(double, double);
double clover_logdet(double);
void compute_ldl_decomp(double);
void compute_clover_term();
void clover_init(double);

#endif
