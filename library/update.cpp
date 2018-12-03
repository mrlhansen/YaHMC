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

void hmc_init()
{
	Monomial *mono;
	ipar = new int_par;
	int_par *ip = ipar;

	// Add monomials and integrators
	if(var_int("act:puregauge"))
	{
		// Gauge
		mono = new MonomialGauge(var_dbl("run:beta"), var_dbl("run:c0"));
		monomials.add(mono, 0);

		// Integrator
		ip->level = 0;
		ip->steps = var_int("int:gsteps");
		ip->integrator = &o2mn_multilevel;
		ip->next = 0;
	}
	else if(var_dbl("act:dm"))
	{
		// Gauge
		mono = new MonomialGauge(var_dbl("run:beta"), var_dbl("run:c0"));
		monomials.add(mono, 0);

		// HMC
		mono = new MonomialHMC(var_dbl("run:mass"), var_dbl("act:dm"), var_dbl("inv:prec"), var_int("mre:past"));
		monomials.add(mono, 1);

		// Hasenbusch
		mono = new MonomialHB(var_dbl("run:mass"), var_dbl("act:dm"), var_dbl("inv:prec"));
		monomials.add(mono, 2);

		// Integrator
		ip->level = 2;
		ip->steps = var_int("int:nsteps");
		ip->integrator = &o2mn_multilevel;
		ip->next = new int_par;
		ip = ip->next;

		ip->level = 1;
		ip->steps = var_int("int:hsteps");
		ip->integrator = &o2mn_multilevel;
		ip->next = new int_par;
		ip = ip->next;

		ip->level = 0;
		ip->steps = var_int("int:gsteps");
		ip->integrator = &o2mn_multilevel;
		ip->next = 0;
	}
	else
	{
		// Gauge
		mono = new MonomialGauge(var_dbl("run:beta"), var_dbl("run:c0"));
		monomials.add(mono, 0);

		// HMC
		mono = new MonomialHMC(var_dbl("run:mass"), 0, var_dbl("inv:prec"), var_int("mre:past"));
		monomials.add(mono, 1);

		// Integrator
		ip->level = 1;
		ip->steps = var_int("int:nsteps");
		ip->integrator = &o2mn_multilevel;
		ip->next = new int_par;
		ip = ip->next;

		ip->level = 0;
		ip->steps = var_int("int:gsteps");
		ip->integrator = &o2mn_multilevel;
		ip->next = 0;
	}

	// Log info
	lprintf("UPDATE", INFO, "Integration steps: %d, length: %1.2f", var_int("int:nsteps"), var_dbl("int:length"));
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

	// Only represent gauge field if necessary
	if(var_int("act:puregauge") == 0)
	{
		represent_gauge_field();
	}
	else
	{
		mp_transfer_links();
	}
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

void update()
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
	ipar->integrator(var_dbl("int:length"), ipar);

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
