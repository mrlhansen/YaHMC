#ifndef BCS_H

#include <spinorfield.h>
#define BC_DEFAULT

#ifdef BC_OPEN
#undef BC_DEFAULT
#define PLAQ_WEIGHTS
#define RECT_WEIGHTS
#define DIRAC_BOUNDARY_TERM
#endif

#ifdef BC_ANTIPERIODIC
#undef BC_DEFAULT
#endif

#ifdef BC_PERIODIC
#undef BC_DEFAULT
#endif

#ifdef BC_SF
#undef BC_DEFAULT
#define PLAQ_WEIGHTS
#define RECT_WEIGHTS
#define DIRAC_BOUNDARY_TERM
#endif

#ifdef BC_DEFAULT
#undef BC_DEFAULT
#define BC_PERIODIC
#endif

#ifdef PLAQ_WEIGHTS
#define plaq_weight(id,mu,nu) plaq_weights[16*(id)+4*(mu)+nu]
#else
#define plaq_weight(id,mu,nu) 1.f
#endif

#ifdef RECT_WEIGHTS
#define rect_weight(id,mu,nu) rect_weights[16*(id)+4*(mu)+nu]
#else
#define rect_weight(id,mu,nu) 1.f
#endif

typedef struct {
	double cf;
	double Ta;
	double Tb;
} bc_par;

extern double *plaq_weights;
extern double *rect_weights;
extern bc_par bcs;

void bcs_init(double);
void apply_bcs_on_spinor_field(SpinorField&);
void apply_bcs_on_momentum_field();
void apply_bcs_on_gauge_field();
void apply_bcs_on_represented_gauge_field();

#endif
