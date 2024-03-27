#include <global.h>
#include <random.h>
#include <repr.h>
#include <bcs.h>
#include <wilsonflow.h>
#include <cstring>
#include <fstream>

string ifn;
string lfn;
string ofn;

void parse_args(int argc, char *argv[])
{
	for(int i = 1; i < argc; i++)
	{
		if(strcmp(argv[i], "-i") == 0)
		{
			ifn = argv[++i];
			continue;
		}
		if(strcmp(argv[i], "-l") == 0)
		{
			lfn = argv[++i];
			continue;
		}
		if(strcmp(argv[i], "-o") == 0)
		{
			ofn = argv[++i];
			continue;
		}
	}

	if(ifn.empty())
	{
		ifn = "wilsonflow.cfg";
	}
}

int main(int argc, char *argv[])
{
	double tmax, tcur;
	double prec, tmeas;
	int outer, inner, ncfg;
	ifstream f;
	string s;

	// Handle arguments
	parse_args(argc, argv);

	// Initialize MPI
	mp_init(argc, argv);

	// Read settings
	var_init(ifn);

	// Enable logger
	logger_init(ofn, var_int("log", "level"));

	// Setup geometry
	geometry_init(1);

	// Setup communications
	mp_setup();

	// Initialize random generator
	rand_init(var_int("rand", "seed"));

	// Setup representation
	repr_init();

	// Initialize boundary conditions
	bcs_init(1.0);

	// Initialize gauge field
	gauge_field_init();

	// Check list file
	if(lfn.empty())
	{
		lprintf("WILSONFLOW", CRITICAL, "No configuration list specified (use -l)");
	}

	// Open file
	f.open(lfn);

	if(f.is_open() == false)
	{
		lprintf("WILSONFLOW", CRITICAL, "Unable to open file: %s", lfn.c_str());
	}

	// Settings
	tmax  = var_dbl("wilsonflow", "tmax");
	prec  = var_dbl("wilsonflow", "prec");
	tmeas = var_dbl("wilsonflow", "tmeas");
	outer = tmax / tmeas;
	inner = tmeas / prec;
	ncfg = 0;

	// Print log info
	lprintf("WILSONFLOW", INFO, "Settings: tmax = %g, tmeas = %g, prec = %g", tmax, tmeas, prec);

	// Analyse configurations
	while(getline(f, s))
	{
		gauge_field_load(s.c_str());
		wf_init();

		tcur = 0;
		ncfg++;

		for(int n = 0; n <= outer; n++)
		{
			lprintf("WILSONFLOW", INFO, "(ncfg,t,E,Esym,TC) = %d %1.3e %1.8e %1.8e %1.8e", ncfg, tcur, wf_e_plaq(), wf_e_sym(), wf_charge());
			for(int m = 0; m < inner; m++)
			{
				wf_iterate(prec);
				tcur += prec;
			}
		}
	}

	return 0;
}
