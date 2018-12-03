#ifndef MRE_H
#define MRE_H

#include <spinorfield.h>
#include <inverter.h>

#define MRE_MAX 15

typedef struct {
	SpinorField s[MRE_MAX];
	mvm_operator mvm;
	double mass;
	int num;
	int max;
	int init;
} mre_par;

void mre_guess(mre_par&, SpinorField&, SpinorField&);
void mre_store(mre_par&, SpinorField&);
void mre_reset(mre_par&);
void mre_init(mre_par&, int, mvm_operator, double);

#endif
