#include <timing.h>
#include <logger.h>
#include <sys/time.h>

static double tm[tm_num_types];
static double tm_total;

static const char *tm_names[tm_num_types] = {
	"linear algebra",
	"hopping term",
	"clover term",
	"gauge force",
	"hopping force",
	"clover force",
	"communication",
};

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
	double total, other;
	total = timestamp() - tm_total;
	other = total;

	lprintf("TIMING", DEBUG, "total time: %1.4f seconds", total);
	for(int i = 0; i < tm_num_types; i++)
	{
		other -= tm[i];
		lprintf("TIMING", DEBUG, "%s: %1.4f seconds (%1.1f%%)", tm_names[i], tm[i], 100*tm[i]/total);
	}
	lprintf("TIMING", DEBUG, "other: %1.4f seconds (%1.1f%%)", other, 100*other/total);
}
