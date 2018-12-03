#include <global.h>
#include <cstdio>
#include <cstdlib>
#include <cstdarg>
#include <cstring>
#include <vector>

typedef struct {
	char sys[32];
	char msg[256];
	int lvl;
} log_t;

static vector<log_t> buf;
static FILE *fp;
static int init = 0;
static int log_level = 0;

void lprintf(const char *sys, int lvl, const char* format, ...)
{
	log_t entry;

	if(init == 0)
	{
		entry.lvl = lvl;
		strcpy(entry.sys, sys);

		va_list args;
		va_start(args, format);
		vsprintf(entry.msg, format, args);
		va_end(args);

		buf.push_back(entry);
		return;
	}

	if(mpi_rank == 0)
	{
		if(lvl <= log_level)
		{
			va_list args;
			fprintf(fp, "[%s][%02d]", sys, lvl);

			va_start(args, format);
			vfprintf(fp, format, args);
			va_end(args);

			fprintf(fp, "\n");
			fflush(fp);
		}
	}

	if(lvl == CRITICAL)
	{
		exit(0);
	}
}

void logger_init(string fn, int lvl)
{
	// Enable logger
	init = 1;
	log_level = lvl;

	// Only master should write output
	if(mpi_rank == 0)
	{
		if(fn.empty())
		{
			fp = stdout;
		}
		else
		{
			fp = fopen(fn.c_str(), "a+");

			if(fp == NULL)
			{
				fp = stdout;
				lprintf("LOG", CRITICAL, "Unable to open log file: %s\n", fn.c_str());
			}

			atexit(logger_exit);
		}
	}

	// Print buffer
	for(log_t entry : buf)
	{
		lprintf(entry.sys, entry.lvl, entry.msg);
	}

	// Clear buffer
	buf.clear();
}

void logger_exit()
{
	fclose(fp);
}
