#ifndef MP_H
#define MP_H

#include <spinorfield.h>

#ifdef _OPENMP
#define ENABLE_OMP
#endif

void mp_init(int, char**);
void mp_setup();

#ifdef ENABLE_MPI
void mp_finalize();
void mp_transfer_spinor(SpinorField&);
void mp_transfer_spinor_even(SpinorField&);
void mp_transfer_spinor_odd(SpinorField&);
void mp_transfer_links();
void mp_transfer_clover_force();
void mp_global_sum(void*, int);
void mp_broadcast(int*);
void mp_barrier();
void mp_gather(void*, void*, int);
void mp_scatter(void*, void*, int);
#else
#define mp_transfer_spinor(a)
#define mp_transfer_spinor_even(a)
#define mp_transfer_spinor_odd(a)
#define mp_transfer_links()
#define mp_transfer_clover_force()
#define mp_global_sum(a,b)
#define mp_broadcast(a)
#define mp_barrier()
#define mp_gather(a,b,c)
#define mp_scatter(a,b,c)
#endif

extern int mpi_size;
extern int mpi_rank;

#endif
