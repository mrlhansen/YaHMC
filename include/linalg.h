#ifndef LINALG_H
#define LINALG_H

#include <complex.h>
#include <spinorfield.h>

void    spinor_allocate(SpinorField&);
void    spinor_zero(SpinorField&);
void    spinor_random(SpinorField&);
double  spinor_sqnorm(SpinorField&);
void    spinor_copy(SpinorField&, SpinorField&);
void    spinor_g5(SpinorField&);
void    spinor_add_assign(SpinorField&, SpinorField&);
void    spinor_sub_assign(SpinorField&, SpinorField&);
void    spinor_mulr(SpinorField&, double, SpinorField&);
void    spinor_mulc(SpinorField&, complex, SpinorField&);
void    spinor_mulr_assign(SpinorField&, double);
void    spinor_mulc_assign(SpinorField&, complex);
void    spinor_mulr_add_assign(SpinorField&, double, SpinorField&);
void    spinor_mulc_add_assign(SpinorField&, complex, SpinorField&);
void    spinor_mulr_sub_assign(SpinorField&, double, SpinorField&);
void    spinor_mulc_sub_assign(SpinorField&, complex, SpinorField&);
complex spinor_product(SpinorField&, SpinorField&);
double  spinor_product_re(SpinorField&, SpinorField&);

void    spinor_select_even(SpinorField&);
void    spinor_select_odd(SpinorField&);
void    spinor_select_full(SpinorField&);
void    spinor_select_default(SpinorField&);

void    spinor_type_even(SpinorField&);
void    spinor_type_odd(SpinorField&);
void    spinor_type_full(SpinorField&);
void    spinor_type_default(SpinorField&);

#endif
