#include <settings.h>
#include <logger.h>
#include <fstream>
#include <cstring>
#include <map>

typedef struct {
	double dbl;
	string str;
	int num;
} value_t;

static map<string,value_t> list;
int var_num_int = 0;
int var_num_mon = 0;

double var_dbl(string section, string key)
{
	key = section + "." + key;
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

int var_int(string section, string key)
{
	key = section + "." + key;
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

string var_str(string section, string key)
{
	key = section + "." + key;
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

const char* var_cstr(string section, string key)
{
	key = section + "." + key;
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

void var_init(string filename)
{
	char key[128];
	char val[128];
	string s;
	ifstream f;
	string section;

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
		if(s.length() == 0)
		{
			continue;
		}

		if(s.find("//") == 0)
		{
			continue;
		}

		if(s.find("#") == 0)
		{
			continue;
		}

		if(sscanf(s.c_str(), "[%[a-z0-9]]", key) == 1)
		{
			section = key;

			if(strcmp(key, "integrator") == 0)
			{
				sprintf(val, "%s%d", key, var_num_int);
				section = val;
				var_num_int++;
			}

			if(strcmp(key, "monomial") == 0)
			{
				sprintf(val, "%s%d", key, var_num_mon);
				section = val;
				var_num_mon++;
			}

			continue;
		}

		if(sscanf(s.c_str(), "%s = %s", key, val) == 2)
		{
			s = section + "." + key;
			list[s].dbl = atof(val);
			list[s].num = atof(val);
			list[s].str = val;
		}
	}

	f.close();
}
