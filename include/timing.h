#ifndef TIMING_H
#define TIMING_H

typedef enum {
	tm_linalg = 0,
	tm_hopping_term,
	tm_clover_term,
	tm_gauge_force,
	tm_hopping_force,
	tm_clover_force,
	tm_mpi,
	tm_num_types,
} tm_type;

double timestamp();
void timing_reset();
void timing_start(tm_type);
void timing_end(tm_type);
void timing_print();

#endif
