#ifndef SPINORFIELD_H
#define SPINORFIELD_H

#include <spinor.h>

struct SpinorField
{
	bool allocated;
	long default_offset;
	long default_sites;
	long offset;
	long sites;
	Spinor *ptr;

	SpinorField()
	{
		ptr = 0;
		allocated = false;
	}

	~SpinorField()
	{
		if(allocated)
		{
			delete[] ptr;
			allocated = false;
		}
	}

	Spinor& operator[](const int n)
	{
		return ptr[n];
	}

	Spinor& at(const int n)
	{
		return ptr[n];
	}
};

#endif
