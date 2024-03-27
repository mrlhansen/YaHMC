#include <update.h>
#include <global.h>
#include <integrator.h>
#include <monomials.h>
#include <random.h>
#include <timing.h>
#include <repr.h>
#include <dirac.h>
#include <bcs.h>
#include <cmath>

static MonomialList monomials;
static int_par *ipar;

int num_cnfg = 0;
int num_accept = 0;

bool detect_improved_gauge()
{
	char section[32];
	double c0;
	string s;

	for(int i = 0; i < var_num_mon; i++)
	{
		sprintf(section, "monomial%d", i);
		s = var_str(section, "type");
		if(s.compare("gauge"))
		{
			continue;
		}

		c0 = var_dbl(section, "c0");
		c0 = fabs(c0 - 1.0);
		return (c0 > 1e-8);
	}

	return false;
}

void hmc_init()
{
	int *m = new int[var_num_int];
	Monomial *mono = 0;
	char section[32];
	int level, steps;
	const char *sc;
	string s;

	// Not found
	if(var_num_int == 0)
	{
		lprintf("UPDATE", CRITICAL, "Integrator: not found");
	}

	if(var_num_mon == 0)
	{
		lprintf("UPDATE", CRITICAL, "Monomials: not found");
	}

	// Mapping
	for(int i = 0; i < var_num_int; i++)
	{
		m[i] = -1;
	}

	for(int i = 0; i < var_num_int; i++)
	{
		sprintf(section, "integrator%d", i);
		level = var_int(section, "level");

		if(level < 0 || level >= var_num_int)
		{
			lprintf("UPDATE", CRITICAL, "Integrator: invalid level: %d", level);
		}
		if(m[level] >= 0)
		{
			lprintf("UPDATE", CRITICAL, "Integrator: duplicated level: %d", level);
		}

		m[level] = i;
	}

	// Integrator
	ipar = new int_par;
	int_par *ip = ipar;

	for(int i = var_num_int-1; i >= 0; i--)
	{
		sprintf(section, "integrator%d", i);
		steps = var_int(section, "steps");
		s = var_str(section, "type");
		sc = s.c_str();

		if(steps <= 0)
		{
			lprintf("UPDATE", CRITICAL, "Integrator: invalid steps: %d", steps);
		}

		ip->level = i;
		ip->steps = steps;

		if(s.compare("o2lf") == 0)
		{
			ip->integrator = &o2lf_multilevel;
		}
		else if(s.compare("o2mn") == 0)
		{
			ip->integrator = &o2mn_multilevel;
		}
		else if(s.compare("o4mn") == 0)
		{
			ip->integrator = &o4mn_multilevel;
		}
		else
		{
			lprintf("UPDATE", CRITICAL, "Integrator: invalid type: %s", sc);
		}

		if(i > 0)
		{
			ip->next = new int_par;
			ip = ip->next;
		}
		else
		{
			ip->next = 0;
		}

		lprintf("UPDATE", INFO, "Integrator: level = %d, type = %s, steps = %d", i, sc, steps);
	}

	// Monomials
	for(int i = 0; i < var_num_int; i++)
	{
		m[i] = 0;
	}

	for(int i = 0; i < var_num_mon; i++)
	{
		sprintf(section, "monomial%d", i);
		level = var_int(section, "level");
		s = var_str(section, "type");
		sc = s.c_str();

		if(level < 0 || level >= var_num_int)
		{
			lprintf("UPDATE", CRITICAL, "Monomial: invalid integrator level: %d", level);
		}

		if(s.compare("gauge") == 0)
		{
			if(level > 0)
			{
				lprintf("UPDATE", CRITICAL, "Monomial: type '%s' must be on integrator level 0", sc);
			}
			mono = new MonomialGauge(var_dbl(section, "beta"), var_dbl(section, "c0"));
		}
		else if(s.compare("hmc") == 0)
		{
			if(level == 0)
			{
				lprintf("UPDATE", CRITICAL, "Monomial: type '%s' cannot be on integrator level 0", sc);
			}
			mono = new MonomialHMC(var_dbl(section, "mass"), var_dbl(section, "dm"), var_dbl(section, "prec"), var_int(section, "mre_past"));
		}
		else if(s.compare("hasenbusch") == 0)
		{
			if(level == 0)
			{
				lprintf("UPDATE", CRITICAL, "Monomial: type '%s' cannot be on integrator level 0", sc);
			}
			mono = new MonomialHB(var_dbl(section, "mass"), var_dbl(section, "dm"), var_dbl(section, "prec"));
		}
		else if(s.compare("rhmc") == 0)
		{
			if(level == 0)
			{
				lprintf("UPDATE", CRITICAL, "Monomial: type '%s' cannot be on integrator level 0", sc);
			}
			mono = new MonomialRHMC(var_dbl(section, "mass"), var_dbl(section, "rprec"), var_dbl(section, "prec"));
		}
		else
		{
			lprintf("UPDATE", CRITICAL, "Monomial: invalid type: %s", sc);
		}

		m[level]++;
		monomials.add(mono, level);
	}

	for(int i = 0; i < var_num_int; i++)
	{
		if(m[i] == 0)
		{
			lprintf("UPDATE", CRITICAL, "Integrator: no monomials on level: %d", i);
		}
	}
}

double hamiltonian()
{
	double value = 0;

	// Boundary conditions
	apply_bcs_on_momentum_field();

	// Contribution from the momenta
	#pragma omp parallel for reduction(+:value)
	sites_for(id)
	{
		value += momentum(id,0).sqnorm();
		value += momentum(id,1).sqnorm();
		value += momentum(id,2).sqnorm();
		value += momentum(id,3).sqnorm();
	}

	mp_global_sum(&value, 1);
	value *= (TF / 2.0);

	// Contribution from monomials
	value += monomials.action();

	// Return action
	return value;
}

void update_links(double eps)
{
	// Boundary conditions
	apply_bcs_on_momentum_field();

	#pragma omp parallel for
	sites_for(id)
	{
		for(int mu = 0; mu < 4; mu++)
		{
			link(id,mu) = exp(momentum(id,mu),eps)*link(id,mu);
		}
	}

	// Boundary conditions
	apply_bcs_on_gauge_field();

	// Communication
	mp_transfer_links();
}

void update_momenta(double dt, int level)
{
	monomials.update(dt, level);
}

void reset_variables()
{
	// Reset momenta
	sites_for(id)
	{
		momentum(id,0).random();
		momentum(id,1).random();
		momentum(id,2).random();
		momentum(id,3).random();
	}

	// Reset monomials
	monomials.reset();
}

void reunitarize()
{
	#pragma omp parallel for
	sites_for(id)
	{
		link(id,0).reunitarize();
		link(id,1).reunitarize();
		link(id,2).reunitarize();
		link(id,3).reunitarize();
	}
}

void update(double length)
{
	double dH, p, ts;
	int accepted;

	// Backup gauge field
	gauge_field_backup();

	// Trajectory starts
	num_cnfg++;
	timing_reset();
	lprintf("UPDATE", INFO, "Trajectory #%d starting", num_cnfg);

	ts = timestamp();   // Current timestamp
	reset_variables();  // Generate new momenta and pseudofermions
	dH = hamiltonian(); // Calculate Hamiltonian value

	// Reset MVM counter
	mvm_reset();

	// Perform integration
	ipar->integrator(length, ipar);

	// Calculate Hamiltonian and time difference
	dH = hamiltonian() - dH;
	ts = timestamp() - ts;
	timing_print();

	// Acceptance probability
	p = (dH < 0) ? 1 : exp(-dH);

	// Trajectory ends
	lprintf("UPDATE", INFO, "Trajectory finished in %1.2f seconds (MVM = %d)", ts, mvm_get());
	lprintf("UPDATE", INFO, "Hamiltonian: dH = %1.6e, exp(-dH) = %1.6e", dH, exp(-dH));

	if(mpi_rank == 0)
	{
		accepted = (rand_uniform() < p);
		mp_broadcast(&accepted);
	}
	else
	{
		accepted = 0;
		mp_broadcast(&accepted);
	}

	if(accepted)
	{
		num_accept++;
		reunitarize();
		lprintf("UPDATE", INFO, "Configuration accepted: %d/%d (%1.2f%%)", num_accept, num_cnfg, 100.0*num_accept/num_cnfg);
	}
	else
	{
		gauge_field_restore();
		represent_gauge_field();
		lprintf("UPDATE", INFO, "Configuration rejected: %d/%d (%1.2f%%)", num_accept, num_cnfg, 100.0*num_accept/num_cnfg);
	}
}
