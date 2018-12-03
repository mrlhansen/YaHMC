#include <monomials.h>

void MonomialList::add(Monomial *ptr, int level)
{
	MonomialItem item;
	item.level = level;
	item.ptr = ptr;
	list.push_back(item);
}

void MonomialList::update(double dt, int level)
{
	if(level)
	{
		begin_fermion_force();
	}

	for(auto& entry : list)
	{
		if(entry.level == level)
		{
			entry.ptr->update(dt);
		}
	}

	if(level)
	{
		end_fermion_force(dt);
	}
}

void MonomialList::reset()
{
	for(auto& entry : list)
	{
		entry.ptr->reset();
	}
}

double MonomialList::action()
{
	double value = 0;
	for(auto& entry : list)
	{
		value += entry.ptr->action();
	}
	return value;
}
