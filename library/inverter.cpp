#include <inverter.h>
#include <dirac.h>
#include <linalg.h>
#include <logger.h>

static SpinorField tmp;
static int init = 0;

int auto_inverter(inv_par &par, SpinorField &dptr, SpinorField &sptr)
{
	int iter = 0;

	// Allocate memory
	if(init == 0)
	{
		spinor_allocate(tmp);
		init = 1;
	}

	// No support for multishift inversion
	if(par.num > 0)
	{
		lprintf("INVERTER", CRITICAL, "auto: no support for multishift inversion");
	}

	// Find the right inverter
	if(par.mvm == Hphi_sq)
	{
		par.mvm = Dphi;
		spinor_g5(sptr);
		iter += bicgstab_inverter(par, tmp, sptr);
		spinor_g5(tmp);
		iter += bicgstab_inverter(par, dptr, tmp);
		spinor_g5(sptr);
		par.mvm = Hphi_sq;
	}
	else if(par.mvm == Hphi)
	{
		par.mvm = Dphi;
		spinor_g5(sptr);
		iter += bicgstab_inverter(par, dptr, sptr);
		spinor_g5(sptr);
		par.mvm = Hphi;
	}
	else if(par.mvm == Dphi)
	{
		iter += bicgstab_inverter(par, dptr, sptr);
	}
	else
	{
		lprintf("INVERTER", CRITICAL, "auto: MVM operator not supported");
	}

	return iter;
}
