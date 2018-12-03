#include <evlimit.h>
#include <logger.h>
#include <linalg.h>
#include <cmath>

static SpinorField b;
static SpinorField stmp;

// Uses inverse power iteration
double smallest_eigenvalue(double mass, mvm_operator mvm)
{
	double norm;
	double eig = 0;
	double old = 0;
	double tmp = 1;

	// Allocate memory
	spinor_allocate(b);
	spinor_allocate(stmp);

	// Use gaussian spinor
	spinor_random(b);

	// Inverter settings
	inv_par par;
	par.prec = 1.0e-6;
	par.mass = mass;
	par.log = 0;
	par.mvm = mvm;

	// Calculate eigenvalue
	while(tmp > 1.0e-3)
	{
		auto_inverter(par, stmp, b);

		norm = spinor_sqnorm(stmp);
		norm = sqrt(norm);
		old = eig;
		eig = norm;

		norm = 1.0/norm;
		spinor_mulr(b, norm, stmp);

		tmp = fabs((eig - old) / eig);
	}

	// Return eigenvalue
	return (1.0 / eig);
}

// Uses power iteration
double largest_eigenvalue(double mass, mvm_operator mvm)
{
	double norm;
	double eig = 0;
	double old = 0;
	double tmp = 1;

	// Allocate memory
	spinor_allocate(b);
	spinor_allocate(stmp);

	// Use gaussian spinor
	spinor_random(b);

	// Calculate eigenvalue
	while(tmp > 1.0e-3)
	{
		mvm(mass, stmp, b);

		norm = spinor_sqnorm(stmp);
		norm = sqrt(norm);
		old = eig;
		eig = norm;
 
		norm = 1.0/norm;
		spinor_mulr(b, norm, stmp);

		tmp = fabs((eig - old) / eig);
	}

	// Return eigenvalue
	return eig;
}

void evlimits(mvm_operator mvm, double mass, double *min, double *max)
{
	*max = largest_eigenvalue(mass, mvm);
	*min = smallest_eigenvalue(mass, mvm);
	lprintf("EVLIMIT", INFO, "Spectral range is [%1.4e,%1.4e]", *min, *max);
}
