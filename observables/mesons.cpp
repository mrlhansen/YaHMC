#include <spectrum.h>
#include <global.h>
#include <linalg.h>
#include <random.h>
#include <cstring>
#include <string>

static SpinorField source;
static SpinorField prop[4*REPR_DIM];

// The gamma matrices below are in fact g5*GAMMA
static dirac_t ID = {{0,1,2,3}, {complex(+1,0), complex(+1,0), complex(-1,0), complex(-1,0)}};
static dirac_t G1 = {{3,2,1,0}, {complex(0,-1), complex(0,-1), complex(0,-1), complex(0,-1)}};
static dirac_t G2 = {{3,2,1,0}, {complex(-1,0), complex(+1,0), complex(-1,0), complex(+1,0)}};
static dirac_t G3 = {{2,3,0,1}, {complex(0,-1), complex(0,+1), complex(0,-1), complex(0,+1)}};
static dirac_t G5 = {{0,1,2,3}, {complex(+1,0), complex(+1,0), complex(+1,0), complex(+1,0)}};
static dirac_t G5G1 = {{3,2,1,0}, {complex(0,-1), complex(0,-1), complex(0,+1), complex(0,+1)}};
static dirac_t G5G2 = {{3,2,1,0}, {complex(-1,0), complex(+1,0), complex(+1,0), complex(-1,0)}};
static dirac_t G5G3 = {{2,3,0,1}, {complex(0,-1), complex(0,+1), complex(0,+1), complex(0,-1)}};
static dirac_t G0G1 = {{1,0,3,2}, {complex(0,-1), complex(0,-1), complex(0,-1), complex(0,-1)}};
static dirac_t G0G2 = {{1,0,3,2}, {complex(-1,0), complex(+1,0), complex(-1,0), complex(+1,0)}};
static dirac_t G0G3 = {{0,1,2,3}, {complex(0,-1), complex(0,+1), complex(0,-1), complex(0,+1)}};
static dirac_t G0G5 = {{2,3,0,1}, {complex(+1,0), complex(+1,0), complex(+1,0), complex(+1,0)}};
static dirac_t G0G5G1 = {{1,0,3,2}, {complex(0,+1), complex(0,+1), complex(0,-1), complex(0,-1)}};
static dirac_t G0G5G2 = {{1,0,3,2}, {complex(+1,0), complex(-1,0), complex(-1,0), complex(+1,0)}};
static dirac_t G0G5G3 = {{0,1,2,3}, {complex(0,+1), complex(0,-1), complex(0,-1), complex(0,+1)}};

// Correlators
static double *corr_id = 0;
static double *corr_g1 = 0;
static double *corr_g2 = 0;
static double *corr_g3 = 0;
static double *corr_g5 = 0;
static double *corr_g5g1 = 0;
static double *corr_g5g2 = 0;
static double *corr_g5g3 = 0;
static double *corr_g0g1 = 0;
static double *corr_g0g2 = 0;
static double *corr_g0g3 = 0;
static double *corr_g0g5 = 0;
static double *corr_g0g5g1 = 0;
static double *corr_g0g5g2 = 0;
static double *corr_g0g5g3 = 0;
static double *corr_g5_g0g5 = 0;

// Settings
static double mes_prec;   // Inverter precision
static double mes_mass;   // Inverter mass
static int    mes_length; // Correlator length
static int    mes_method; // Calculation method
static int    mes_hits;   // Number of sources
static int    mes_norm;   // Correlator normalization
static int    mes_ninv;   // Number of inversions needed for the source

static double contract_at_index(int id, dirac_t &g1, dirac_t &g2)
{
	complex sum, gtmp;
	complex ctmp1, ctmp2;

	if(mes_method == SEMWALL)
	{
		for(int a = 0; a < 4; a++)
		{
			for(int A = 0; A < 4; A++)
			{
				int b = g1.col[a];
				int B = g2.col[A];
				gtmp = g1.val[a] * g2.val[A];

				for(int x = 0; x < REPR_DIM; x++)
				{
					ctmp1 = prop[A].at(id).get(b,x);
					ctmp2 = prop[B].at(id).get(a,x);
					sum += gtmp * ctmp1 * conj(ctmp2);
				}
			}
		}
	}

	if(mes_method == POINT)
	{
		for(int a = 0; a < 4; a++)
		{
			for(int c = 0; c < 4; c++)
			{
				int b = g1.col[a];
				int d = g2.col[c];
				gtmp = g1.val[a] * g2.val[c];

				for(int x = 0; x < REPR_DIM; x++)
				for(int y = 0; y < REPR_DIM; y++)
				{
					ctmp1 = prop[c*REPR_DIM+x].at(id).get(b,y);
					ctmp2 = prop[d*REPR_DIM+x].at(id).get(a,y);
					sum += gtmp * ctmp1 * conj(ctmp2);
				}
			}
		}
	}

	return sum.re;
}

static void contract_all_correlators(int tau)
{
	sites_for(id)
	{
		int t = global_time(id);
		t = (t - tau + global.dim_t) % global.dim_t;

		corr_id[t] += contract_at_index(id, ID, ID);
		corr_g1[t] += contract_at_index(id, G1, G1);
		corr_g2[t] += contract_at_index(id, G2, G2);
		corr_g3[t] += contract_at_index(id, G3, G3);
		corr_g5[t] += contract_at_index(id, G5, G5);
		corr_g5g1[t] += contract_at_index(id, G5G1, G5G1);
		corr_g5g2[t] += contract_at_index(id, G5G2, G5G2);
		corr_g5g3[t] += contract_at_index(id, G5G3, G5G3);
		corr_g0g1[t] += contract_at_index(id, G0G1, G0G1);
		corr_g0g2[t] += contract_at_index(id, G0G2, G0G2);
		corr_g0g3[t] += contract_at_index(id, G0G3, G0G3);
		corr_g0g5[t] += contract_at_index(id, G0G5, G0G5);
		corr_g0g5g1[t] += contract_at_index(id, G0G5G1, G0G5G1);
		corr_g0g5g2[t] += contract_at_index(id, G0G5G2, G0G5G2);
		corr_g0g5g3[t] += contract_at_index(id, G0G5G3, G0G5G3);
		corr_g5_g0g5[t] += contract_at_index(id, G0G5, G5);
	}
}

static void print_correlator(const char *ch, double *v, double sign)
{
	string str;
	char buf[32];
	double sf;

	// Sign and normalization
	sf = sign/(mes_norm*mes_hits);

	// Implode to string
	for(int i = 0; i < mes_length; i++)
	{
		sprintf(buf, "%1.6e ", sf*v[i]);
		str.append(buf);
	}

	// Print correlator
	lprintf("MESONS", INFO, "%s: %s", ch, str.c_str());
}

void meson_init(double mass, double prec, int hits, int method)
{
	// Settings
	mes_ninv = 0;
	mes_mass = mass;
	mes_prec = prec;
	mes_hits = hits;
	mes_method = method;
	mes_length = global.dim_t;

	// SEMWALL source
	if(mes_method == SEMWALL)
	{
		lprintf("MESONS", INFO, "Calculation method: SEMWALL");
		mes_norm = global.vol3*global.vol3;
		mes_ninv = 4;
	}
	else
	{
		mes_method = POINT;
	}

	// Point source
	if(mes_method == POINT)
	{
		lprintf("MESONS", INFO, "Calculation method: point sources");
		mes_norm = global.vol3;
		mes_ninv = 4*REPR_DIM;
	}

	// Allocate source
	spinor_allocate(source);

	// Allocate propagator
	for(int n = 0; n < mes_ninv; n++)
	{
		spinor_allocate(prop[n]);
	}

	// Allocate correlators
	corr_id = new double[mes_length];
	corr_g1 = new double[mes_length];
	corr_g2 = new double[mes_length];
	corr_g3 = new double[mes_length];
	corr_g5 = new double[mes_length];
	corr_g5g1 = new double[mes_length];
	corr_g5g2 = new double[mes_length];
	corr_g5g3 = new double[mes_length];
	corr_g0g1 = new double[mes_length];
	corr_g0g2 = new double[mes_length];
	corr_g0g3 = new double[mes_length];
	corr_g0g5 = new double[mes_length];
	corr_g0g5g1 = new double[mes_length];
	corr_g0g5g2 = new double[mes_length];
	corr_g0g5g3 = new double[mes_length];
	corr_g5_g0g5 = new double[mes_length];

	// Print log info
	lprintf("MESONS", INFO, "Averaging over %d sources", mes_hits);
}

void meson_measure()
{
	int tau, memsz;

	// Reset correlators
	memsz = mes_length * sizeof(double);
	memset(corr_id, 0, memsz);
	memset(corr_g1, 0, memsz);
	memset(corr_g2, 0, memsz);
	memset(corr_g3, 0, memsz);
	memset(corr_g5, 0, memsz);
	memset(corr_g5g1, 0, memsz);
	memset(corr_g5g2, 0, memsz);
	memset(corr_g5g3, 0, memsz);
	memset(corr_g0g1, 0, memsz);
	memset(corr_g0g2, 0, memsz);
	memset(corr_g0g3, 0, memsz);
	memset(corr_g0g5, 0, memsz);
	memset(corr_g0g5g1, 0, memsz);
	memset(corr_g0g5g2, 0, memsz);
	memset(corr_g0g5g3, 0, memsz);
	memset(corr_g5_g0g5, 0, memsz);

	// Calculate the propagators and perform contractions
	for(int n = 0; n < mes_hits; n++)
	{
		tau = random_timeslice();

		if(mes_method == SEMWALL)
		{
			lprintf("MESONS", INFO, "Using stochastic source on timeslice %d", tau);
		}
		else
		{
			lprintf("MESONS", INFO, "Using point source on timeslice %d", tau);
		}

		for(int n = 0; n < mes_ninv; n++)
		{
			if(mes_method == SEMWALL)
			{
				create_semwall_source(n, tau, source);
			}
			else
			{
				create_point_source(n, tau, source);
			}
			calc_propagator(mes_mass, mes_prec, prop[n], source);
		}

		contract_all_correlators(tau);
	}

	// Collect results
	mp_global_sum(corr_id, mes_length);
	mp_global_sum(corr_g1, mes_length);
	mp_global_sum(corr_g2, mes_length);
	mp_global_sum(corr_g3, mes_length);
	mp_global_sum(corr_g5, mes_length);
	mp_global_sum(corr_g5g1, mes_length);
	mp_global_sum(corr_g5g2, mes_length);
	mp_global_sum(corr_g5g3, mes_length);
	mp_global_sum(corr_g0g1, mes_length);
	mp_global_sum(corr_g0g2, mes_length);
	mp_global_sum(corr_g0g3, mes_length);
	mp_global_sum(corr_g0g5, mes_length);
	mp_global_sum(corr_g0g5g1, mes_length);
	mp_global_sum(corr_g0g5g2, mes_length);
	mp_global_sum(corr_g0g5g3, mes_length);
	mp_global_sum(corr_g5_g0g5, mes_length);

	// Print correlators
	print_correlator("id", corr_id, -1.0);
	print_correlator("g1", corr_g1, -1.0);
	print_correlator("g2", corr_g2, -1.0);
	print_correlator("g3", corr_g3, -1.0);
	print_correlator("g5", corr_g5, +1.0);
	print_correlator("g5g1", corr_g5g1, -1.0);
	print_correlator("g5g2", corr_g5g2, -1.0);
	print_correlator("g5g3", corr_g5g3, -1.0);
	print_correlator("g0g1", corr_g0g1, -1.0);
	print_correlator("g0g2", corr_g0g2, -1.0);
	print_correlator("g0g3", corr_g0g3, -1.0);
	print_correlator("g0g5", corr_g0g5, +1.0);
	print_correlator("g0g5g1", corr_g0g5g2, +1.0);
	print_correlator("g0g5g2", corr_g0g5g3, +1.0);
	print_correlator("g0g5g3", corr_g0g5g3, +1.0);
	print_correlator("g5_g0g5", corr_g5_g0g5, +1.0);
}
