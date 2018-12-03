#ifndef RANDOM_H
#define RANDOM_H

#include <gsl/gsl_rng.h>
#include <gsl/gsl_randist.h>

extern gsl_rng *rng_state;

void rand_init(int);
int rand_store(void*);
void rand_load(void*, int);

#define SQRT2INV 0.70710678118654752440084436210
#define rand_uniform() gsl_rng_uniform_pos(rng_state)
#define rand_gaussian() gsl_ran_gaussian_ziggurat(rng_state, SQRT2INV)
#define rand_z2() ((rand_uniform() < 0.5) ? SQRT2INV : -SQRT2INV)

#endif
