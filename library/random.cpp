#include <random.h>
#include <global.h>
#include <cstring>
#include <x86intrin.h>

gsl_rng *rng_state;

static int random_seed()
{
	unsigned long seed;
	seed = __rdtsc();
	return (seed % 1048576);
}

static void rand_exit()
{
	gsl_rng_free(rng_state);
}

void rand_init(int seed)
{
	gsl_rng_env_setup();

	if(seed == 0)
	{
		seed = random_seed();
	}

	rng_state = gsl_rng_alloc(gsl_rng_ranlxd1);
	gsl_rng_set(rng_state, seed + mpi_rank);
	lprintf("RANDOM", INFO, "Random generator initialized with seed %d", seed);

	atexit(rand_exit);
}

int rand_store(void *ptr)
{
	int sz = gsl_rng_size(rng_state);
	memcpy(ptr, gsl_rng_state(rng_state), sz);
	return sz;
}

void rand_load(void *ptr, int sz)
{
	memcpy(gsl_rng_state(rng_state), ptr, sz);
	lprintf("RANDOM", INFO, "Random generator continued from previous state");
}
