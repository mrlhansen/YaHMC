#include <global.h>
#include <memory.h>

suNg *suNg_allocate(int ndir)
{
	suNg *ptr;
	int memsz;

	ptr = 0;
	memsz = ndir * outer.vol4 * sizeof(suNg);
	memsz /= 1024;

	try
	{
		ptr = new suNg[ndir*outer.vol4];
	}
	catch(std::bad_alloc&)
	{
		lprintf("MEMORY", CRITICAL, "Failed to allocate %lu kb of memory", memsz);
	}

	lprintf("MEMORY", DEBUG, "Allocated %lu kb of memory", memsz);
	return ptr;
}

suNf *suNf_allocate(int ndir)
{
	suNf *ptr;
	int memsz;

	ptr = 0;
	memsz = ndir * outer.vol4 * sizeof(suNf);
	memsz /= 1024;

	try
	{
		ptr = new suNf[ndir*outer.vol4];
	}
	catch(std::bad_alloc&)
	{
		lprintf("MEMORY", CRITICAL, "Failed to allocate %lu kb of memory", memsz);
	}

	lprintf("MEMORY", DEBUG, "Allocated %lu kb of memory", memsz);
	return ptr;
}

suNa *suNa_allocate(int ndir)
{
	suNa *ptr;
	int memsz;

	ptr = 0;
	memsz = ndir * outer.vol4 * sizeof(suNa);
	memsz /= 1024;

	try
	{
		ptr = new suNa[ndir*outer.vol4];
	}
	catch(std::bad_alloc&)
	{
		lprintf("MEMORY", CRITICAL, "Failed to allocate %lu kb of memory", memsz);
	}

	lprintf("MEMORY", DEBUG, "Allocated %lu kb of memory", memsz);
	return ptr;
}

ldl_t *ldl_allocate(int ndir)
{
	ldl_t *ptr;
	int memsz;

	ptr = 0;
	memsz = ndir * outer.vol4 * sizeof(ldl_t);
	memsz /= 1024;

	try
	{
		ptr = new ldl_t[ndir*outer.vol4];
	}
	catch(std::bad_alloc&)
	{
		lprintf("MEMORY", CRITICAL, "Failed to allocate %lu kb of memory", memsz);
	}

	lprintf("MEMORY", DEBUG, "Allocated %lu kb of memory", memsz);
	return ptr;
}
