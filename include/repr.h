#ifndef REPR_H
#define REPR_H

// The generators of the fundamental representation
// must have trace normalization 2
#define TF 2
#define NG (NC*NC - 1)

// Fundamental
#ifdef REPR_FUNDAMENTAL
#define REPR_DIM NC
#define REPR_NAME "Fundamental"
#define REPR_TRACE 2
#define REPR_ID 0x00
#endif

// Adjoint
#ifdef REPR_ADJOINT
#define REPR_DIM (NC*NC-1)
#define REPR_NAME "Adjoint"
#define REPR_TRACE (4*NC)
#define REPR_ID 0x01
#endif

// Symmetric
#ifdef REPR_SYMMETRIC
#define REPR_DIM (NC*(NC+1)/2)
#define REPR_NAME "Symmetric"
#define REPR_TRACE (2*(NC+2))
#define REPR_ID 0x02
#endif

// Antisymmetric
#ifdef REPR_ANTISYMMETRIC
#define REPR_DIM (NC*(NC-1)/2)
#define REPR_NAME "Antisymmetric"
#define REPR_TRACE (2*(NC-2))
#define REPR_ID 0x03
#endif

#ifndef REPR_NAME
#error Unknown or unspecified representation
#endif

void repr_init();
void represent_gauge_field();

#endif
