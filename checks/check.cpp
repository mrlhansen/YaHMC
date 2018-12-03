#include <checks.h>
#include <logger.h>

void run_checks()
{
	int failed = 0;

	lprintf("CHECK", INFO, "Running checks now");

	failed += check_group();
	failed += check_dirac();
	failed += check_inverters();

	if(failed == 0)
	{
		lprintf("CHECK", INFO, "Status: PASSED (all checks succeded)");
	}
	else if(failed == 1)
	{
		lprintf("CHECK", INFO, "Status: FAILED (%d check reported an error)", failed);
	}
	else
	{
		lprintf("CHECK", INFO, "Status: FAILED (%d checks reported an error)", failed);
	}
}
