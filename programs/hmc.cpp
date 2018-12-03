// Martin Hansen (martin.hansen@roma2.infn.it)
// December 2, 2018
// Copyright (C), all rights reserved.

#include <global.h>
#include <random.h>
#include <staples.h>
#include <update.h>
#include <wilsonloop.h>
#include <polyakov.h>
#include <spectrum.h>
#include <checks.h>
#include <bcs.h>
#include <clover.h>
#include <cstring>
#include <cmath>

static string ifn;
static string ofn;
static int checks = 0;

static void parse_args(int argc, char *argv[])
{
	for(int i = 1; i < argc; i++)
	{
		if(strcmp(argv[i], "-i") == 0)
		{
			ifn = argv[++i];
			continue;
		}
		if(strcmp(argv[i], "-o") == 0)
		{
			ofn = argv[++i];
			continue;
		}
		if(strcmp(argv[i], "--run-checks") == 0)
		{
			checks = 1;
		}
	}

	if(ifn.empty())
	{
		ifn = "hmc.cfg";
	}
}

static void load_configuration()
{
	char filename[256];
	sprintf(filename, "%s/%s", var_cstr("save:dir"), var_cstr("run:start"));
	gauge_field_load(filename);
}

static void save_configuration()
{
	char filename[256];
	sprintf(filename, "%s/%s_%dx%dx%dx%d_nc%db%1.4fm%1.4fn%04d",
			  var_cstr("save:dir"),
			  var_cstr("save:name"),
			  global.dim_t,
			  global.dim_x,
			  global.dim_y,
			  global.dim_z,
			  NC,
			  var_dbl("run:beta"),
			  fabs(var_dbl("run:mass")),
			  num_cnfg);
	gauge_field_save(filename);
}

int main(int argc, char *argv[])
{
	complex val;
	int border;

	// Handle arguments
	parse_args(argc, argv);

	// Initialize MPI
	mp_init(argc, argv);

	// Read settings
	var_init(ifn);

	// Enable logger
	logger_init(ofn, var_int("log:level"));

	// Determine border size
	if(var_dbl("run:c0") == 1.0)
	{
		border = 1;
	}
	else
	{
		border = 2;
	}

	// Setup geometry
	geometry_init(border);

	// Setup communications
	mp_setup();

	// Initialize random generator
	rand_init(var_int("rand:seed"));

	// Setup representation
	repr_init();

	// Initialize boundary conditions
	bcs_init(var_dbl("bc:cf"));

	// Initialize clover term
	clover_init(var_dbl("run:csw"));

	// Initialize gauge field
	gauge_field_init();

	// Load or set gauge field
	if(var_str("run:start").compare("random") == 0)
	{
		gauge_field_set_random();
	}
	else if(var_str("run:start").compare("unit") == 0)
	{
		gauge_field_set_unit();
	}
	else
	{
		load_configuration();
	}

	// Test suite
	if(checks)
	{
		run_checks();
		return 0;
	}

	// Perform staple test
	staples_test();

	// Initialize algorithm
	hmc_init();

	// Initialize mesonic correlators
	if(var_int("mes:freq"))
	{
		meson_init(var_dbl("run:mass"), var_dbl("mes:prec"), var_int("mes:hits"), var_int("mes:method"));
	}

	// Generate configurations
	while(num_cnfg < var_int("run:last"))
	{
		update();

		if(var_int("plaq:freq"))
		{
			if((num_cnfg % var_int("plaq:freq")) == 0)
			{
				lprintf("HMC", INFO, "Plaquette: %1.6e", avr_wilson(1,1));
			}
		}

		if(var_int("poly:freq"))
		{
			if((num_cnfg % var_int("poly:freq")) == 0)
			{
				val = avr_polyakov(0);
				lprintf("HMC", INFO, "Polyakov, mu = 0: %+1.6e %+1.6e", val.re, val.im);
				val = avr_polyakov(1);
				lprintf("HMC", INFO, "Polyakov, mu = 1: %+1.6e %+1.6e", val.re, val.im);
				val = avr_polyakov(2);
				lprintf("HMC", INFO, "Polyakov, mu = 2: %+1.6e %+1.6e", val.re, val.im);
				val = avr_polyakov(3);
				lprintf("HMC", INFO, "Polyakov, mu = 3: %+1.6e %+1.6e", val.re, val.im);
			}
		}

		if(var_int("mes:freq"))
		{
			if((num_cnfg % var_int("mes:freq")) == 0)
			{
				meson_measure();
			}
		}

		if(var_int("save:freq"))
		{
			if((num_cnfg % var_int("save:freq")) == 0)
			{
				save_configuration();
			}
		}
	}

	return 0;
}
