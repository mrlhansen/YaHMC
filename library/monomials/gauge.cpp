#include <monomials.h>
#include <wilsonloop.h>
#include <global.h>
#include <cmath>

MonomialGauge::MonomialGauge(double _beta, double _c0)
{
	beta = _beta;
	c0 = _c0;
	c1 = (1.0-c0)/8.0;
	improved = (fabs(c1) > 1e-8);
	lprintf("MONOMIAL", INFO, "Gauge: beta = %1.6f, c0 = %1.6f, c1 = %1.6f", beta, c0, c1);
}

void MonomialGauge::reset()
{
	return;
}

double MonomialGauge::action()
{
	double value;
	double v0, v1;

	if(improved)
	{
		v0 = avr_wilson(1,1);
		v1 = avr_wilson(2,1) + avr_wilson(1,2);
		value = beta * (c0*v0 + c1*v1) * 6 * volume;
	}
	else
	{
		v0 = avr_wilson(1,1);
		value = beta * v0 * 6 * volume;
	}

	return -value;
}

void MonomialGauge::update(double dt)
{
	if(improved)
	{
		improved_gauge_force(dt, beta, c0, c1);
	}
	else
	{
		wilson_gauge_force(dt, beta);
	}
}
