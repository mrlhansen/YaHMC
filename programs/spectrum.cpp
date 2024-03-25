// Martin Hansen (martin.hansen@roma2.infn.it)
// December 2, 2018
// Copyright (C), all rights reserved.

#include <spectrum.h>
#include <clover.h>
#include <global.h>
#include <random.h>
#include <repr.h>
#include <bcs.h>
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
		ifn = "spectrum.cfg";
	}
}

int main(int argc, char *argv[])
{
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
	bcs_init(var_dbl("clover", "cf"));

	// Initialize clover term
	clover_init(var_dbl("clover", "csw"));

	// Initialize gauge field
	gauge_field_init();

	// Initialize mesons
	meson_init(var_dbl("mesons", "mass"), var_dbl("mesons", "prec"), var_int("mesons", "hits"), var_int("mesons", "method"));

	// Check list file
	if(lfn.empty())
	{
		lprintf("SPECTRUM", CRITICAL, "No configuration list specified (use -l)");
	}

	// Open file
	f.open(lfn);

	if(f.is_open() == false)
	{
		lprintf("SPECTRUM", CRITICAL, "Unable to open file: %s", lfn.c_str());
	}

	// Analyse configurations
	while(getline(f, s))
	{
		gauge_field_load(s.c_str());
		meson_measure();
	}

	return 0;
}
