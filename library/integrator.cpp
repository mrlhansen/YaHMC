#include <integrator.h>
#include <update.h>
#include <repr.h>

// Update dynamical fields
static void update_fields(double dt, int_par *par)
{
	if(par->next == 0)
	{
		update_links(dt);
		set_gauge_field_state(true);
	}
	else
	{
		par->next->integrator(dt, par->next);
		set_gauge_field_state(false);
	}
}

// 2nd order leap-frog integrator
void o2lf_multilevel(double len, int_par *par)
{
	double dt = (len / par->steps);

	for(int n = 0; n < par->steps; n++)
	{
		if(n == 0)
		{
			update_momenta(dt/2, par->level);
		}
		else
		{
			update_momenta(dt, par->level);
		}
		update_fields(dt, par);
	}

	update_momenta(dt/2, par->level);
}

// 2nd order minimal norm integrator
void o2mn_multilevel(double len, int_par *par)
{
	const double a = 0.1931833275037836;
	double dt = (len / par->steps);

	for(int n = 0; n < par->steps; n++)
	{
		if(n == 0)
		{
			update_momenta(a*dt, par->level);
		}
		else
		{
			update_momenta(2*a*dt, par->level);
		}
		update_fields(dt/2, par);
		update_momenta((1-2*a)*dt, par->level);
		update_fields(dt/2, par);
	}

	update_momenta(a*dt, par->level);
}

// 4th order minimal norm integrator
void o4mn_multilevel(double len, int_par *par)
{
	const double a = 0.08398315262876693;
	const double b = 0.2539785108410595;
	const double c = 0.6822365335719091;
	const double d = -0.03230286765269967;
	double dt = (len / par->steps);

	for(int n = 0; n < par->steps; n++)
	{
		if(n == 0)
		{
			update_momenta(a*dt, par->level);
		}
		else
		{
			update_momenta(2*a*dt, par->level);
		}
		update_fields(b*dt, par);
		update_momenta(c*dt, par->level);
		update_fields(d*dt, par);
		update_momenta((0.5-a-c)*dt, par->level);
		update_fields((1-2*(b+d))*dt, par);
		update_momenta((0.5-a-c)*dt, par->level);
		update_fields(d*dt, par);
		update_momenta(c*dt, par->level);
		update_fields(b*dt, par);
	}

	update_momenta(a*dt, par->level);
}
