#include <settings.h>
#include <logger.h>
#include <fstream>
#include <map>

typedef struct {
	double dbl;
	string str;
	int num;
} value_t;

static map<string,value_t> list;

double var_dbl(string key)
{
	if(list.find(key) == list.end())
	{
		lprintf("SETTINGS", WARNING, "Parameter '%s' not found in settings", key.c_str());
		return 0;
	}
	else
	{
		return list[key].dbl;
	}
}

void var_dbl(string key, double val)
{
	list[key].dbl = val;
}

int var_int(string key)
{
	if(list.find(key) == list.end())
	{
		lprintf("SETTINGS", WARNING, "Parameter '%s' not found in settings", key.c_str());
		return 0;
	}
	else
	{
		return list[key].num;
	}
}

void var_int(string key, int val)
{
	list[key].num = val;
}

string var_str(string key)
{
	if(list.find(key) == list.end())
	{
		lprintf("SETTINGS", WARNING, "Parameter '%s' not found in settings", key.c_str());
		return "";
	}
	else
	{
		return list[key].str;
	}
}

const char* var_cstr(string key)
{
	if(list.find(key) == list.end())
	{
		lprintf("SETTINGS", WARNING, "Parameter '%s' not found in settings", key.c_str());
		return 0;
	}
	else
	{
		return list[key].str.c_str();
	}
}

void var_str(string key, string val)
{
	list[key].str = val;
}

void var_init(string filename)
{
	char key[128];
	char val[128];
	string s;
	ifstream f;

	f.open(filename);

	if(f.is_open() == false)
	{
		lprintf("SETTINGS", CRITICAL, "Unable to open file: %s", filename.c_str());
	}
	else
	{
		lprintf("SETTINGS", INFO, "Reading settings from: %s", filename.c_str());
	}

	while(getline(f, s))
	{
		if(s.find("//") == 0)
		{
			continue;
		}

		if(sscanf(s.c_str(), "%s = %s", key, val) == 2)
		{
			list[key].dbl = atof(val);
			list[key].num = atof(val);
			list[key].str = val;
		}
	}

	f.close();
}
