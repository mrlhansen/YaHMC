#include <rational.h>
#include <global.h>
#include <approx.h>
#include <linalg.h>
#include <inverter.h>
#include <cstring>
#include <cmath>
#include <cstdlib>

double residue(rational_app &data, int n)
{
	double tmp;
	double p = 0;
	double q = 1;
	double x = data.pole[n];

	for(int i = 0; i < data.order; i++)
	{
		q *= (data.root[i] - x);
	}

	for(int i = 0; i < data.order; i++)
	{
		tmp = 1;
		for(int j = 0; j < data.order; j++)
		{
			if(i != j)
			{
				tmp *= (data.pole[j] - x);
			}
		}
		p += tmp;
	}

	return (q / p) * data.a0;
}

void rescale(rational_app &data, double max)
{
	double k = (max / data.max);
	double fexp = data.n;

	fexp /= data.d;
	data.a0 *= pow(k, fexp);

	for(int i = 0; i < data.order; i++)
	{
		data.root[i] *= k;
		data.pole[i] *= k;
	}

	data.min /= k;
	data.max /= k;
}

void rational_find(rational_app &data, int n, int d, double prec, double min, double max)
{
	int sz = sizeof(appdb)/sizeof(double);
	double tprec, tmin, tmax, root, pole;
	int to, tn, td;
	int bo = 0;
	int bi = -1;

	// Search for a valid approximation
	for(int i = 0; i < sz; i += 7+2*to)
	{
		tn    = appdb[i+0]; // Numerator
		td    = appdb[i+1]; // Denominator
		to    = appdb[i+2]; // Polynomial order
		tprec = appdb[i+3]; // Precision
		tmin  = appdb[i+4]; // Interval start
		tmax  = appdb[i+5]; // Interval end

		// Check for a match
		if(tn == abs(n) && td == abs(d) && (tmin/tmax < min/max) && (tprec < prec))
		{
			if(bo == 0 || to < bo)
			{
				bo = to;
				bi = i;
			}
		}
	}

	// If no approximation was found
	if(bi < 0)
	{
		lprintf("RATIONAL", WARNING, "Unable to find rational approximation for x^(%d/%d)", n, d);
		lprintf("RATIONAL", CRITICAL, "Needed precision is %1.6e in [%1.3e,%1.3e]", prec, min, max);
	}

	// Store the best approximation
	data.n     = n;
	data.d     = d;
	data.order = bo;
	data.prec  = appdb[bi+3];
	data.min   = appdb[bi+4];
	data.max   = appdb[bi+5];
	data.a0    = appdb[bi+6];

	// Print info
	lprintf("RATIONAL", INFO, "Best approximation for x^(%+d/%d): order = %d, precision = %1.6e", n, d, bo, data.prec);

	// Store roots and poles
	for(int i = 0; i < bo; i++)
	{
		data.root[i] = appdb[7+bi+i];
		data.pole[i] = appdb[7+bi+bo+i];
	}

	// If the fraction is positive we must swap roots and poles
	if(d*n > 0)
	{
		for(int i = 0; i < bo; i++)
		{
			root = data.root[i];
			pole = data.pole[i];
			data.root[i] = pole;
			data.pole[i] = root;
		}
		data.a0 = (1.0 / data.a0);
	}

	// Rescale to the given interval
	rescale(data, max);

	// Calculate the partial fraction residues
	for(int i = 0; i < bo; i++)
	{
		data.residue[i] = residue(data, i);
	}
}

void rational_calc(rational_app &data, mvm_operator mvm, double mass, double precision, SpinorField &dptr, SpinorField &sptr)
{
	SpinorField *inv_out;
	inv_par par;

	// Inverter parameters
	par.num = data.order;
	par.prec = precision;
	par.shift = data.pole;
	par.mass = mass;
	par.log = 1;
	par.mvm = mvm;

	// Allocate spinor fields
	inv_out = new SpinorField[data.order];
	for(int i = 0; i < data.order; i++)
	{
		spinor_allocate(inv_out[i]);
	}

	// Call inverter
	cg_mshift(par, inv_out, sptr);

	// Calculate solution
	spinor_mulr(dptr, data.a0, sptr);
	for(int i = 0; i < data.order; i++)
	{
		spinor_mulr_add_assign(dptr, data.residue[i], inv_out[i]);
	}

	// Deallocate spinor fields
	delete[] inv_out;
}
