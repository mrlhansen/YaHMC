#ifndef DIRAC_H
#define DIRAC_H

#include <spinorfield.h>

void mvm_reset();
int  mvm_get();

void Dphi_oe(SpinorField&, SpinorField&);
void Dphi_eo(SpinorField&, SpinorField&);
void Dphi_ee(double, SpinorField&, SpinorField&);
void Dphi_oo(double, SpinorField&, SpinorField&);
void Dphi_oo_inv(double, SpinorField&, SpinorField&);
void Dphi(double, SpinorField&, SpinorField&);
void Hphi(double, SpinorField&, SpinorField&);
void Hphi_sq(double, SpinorField&, SpinorField&);

#endif
