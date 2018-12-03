#ifndef GLOBAL_H
#define GLOBAL_H

#include <settings.h>
#include <geometry.h>
#include <gaugefield.h>
#include <logger.h>
#include <mp.h>

// Spinor blocks
#define EVEN 0x01
#define ODD  0x02
#define BOTH 0x03

// Configurations
extern int num_cnfg;
extern int num_accept;

// Generators etc
extern suNg iTfund[NG];
extern suNf iTrepr[NG];
void represent_links(suNf*, suNg*, int);

#endif
