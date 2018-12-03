#include <timing.h>
#include <logger.h>
#include <sys/time.h>

static double tm[tm_num_types];
static double tm_total;

double timestamp()
{
	struct timeval tm;
	gettimeofday(&tm, NULL);
	return tm.tv_sec + 1.0e-6 * tm.tv_usec;
}

void timing_reset()
{
	for(int n = 0; n < tm_num_types; n++)
	{
		tm[n] = 0;
	}
	tm_total = timestamp();
}

void timing_start(tm_type t)
{
	tm[t] -= timestamp();
}

void timing_end(tm_type t)
{
	tm[t] += timestamp();
}

void timing_print()
{
	double total;
	total = timestamp() - tm_total;

	lprintf("TIMING", DEBUG, "total time: %1.4f seconds", total);
	lprintf("TIMING", DEBUG, "linear algebra: %1.4f seconds", tm[tm_linalg]);
	lprintf("TIMING", DEBUG, "dirac operator: %1.4f seconds", tm[tm_dirac]);
	lprintf("TIMING", DEBUG, "gauge force: %1.4f seconds", tm[tm_gforce]);
	lprintf("TIMING", DEBUG, "hopping force: %1.4f seconds", tm[tm_wforce]);
	lprintf("TIMING", DEBUG, "clover force: %1.4f seconds", tm[tm_cforce]);
	lprintf("TIMING", DEBUG, "communication: %1.4f seconds", tm[tm_mpi]);
}
