#ifndef TIMING_H
#define TIMING_H

typedef enum {
	tm_mpi = 0,
	tm_dirac,
	tm_gforce,
	tm_wforce,
	tm_cforce,
	tm_linalg,
	tm_num_types,
} tm_type;

double timestamp();
void timing_reset();
void timing_start(tm_type);
void timing_end(tm_type);
void timing_print();

#endif
