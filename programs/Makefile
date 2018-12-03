# Martin Hansen (martin.hansen@roma2.infn.it)
# December 2, 2018
# Copyright (C), all rights reserved.

# Variables
CC = g++
LDFLAGS = -lm -lgsl
CCFLAGS = -Wall -O3 -mfpmath=sse -msse4 -Wno-unknown-pragmas

# Colors and representation
COLORS = 3
REPR = FUNDAMENTAL
#REPR = SYMMETRIC
#REPR = ANTISYMMETRIC
#REPR = ADJOINT

# Enable OpenMP
LDFLAGS += -fopenmp
CCFLAGS += -fopenmp

# Enable MPI
#CC = mpicxx
#CCFLAGS += -DENABLE_MPI

# Enable clover term
CCFLAGS += -DCLOVER_TERM

# Enable even-odd preconditioning
CCFLAGS += -DUPDATE_EO

# Enable vectorization
#CCFLAGS += -DENABLE_AVX -mavx

# Boundary conditions
CCFLAGS += -DBC_PERIODIC
#CCFLAGS += -DBC_ANTIPERIODIC
#CCFLAGS += -DBC_OPEN
#CCFLAGS += -DBC_SF

# Do not change the following
SOURCES := $(shell find ../library ../checks ../observables -type f -name '*.cpp')
SOURCES += $(wildcard *.cpp)
OBJECTS := $(SOURCES:.cpp=.o)
TARGETS := $(wildcard *.cpp)
TARGETS := $(TARGETS:.cpp=)

CCFLAGS += -I../include
CCFLAGS += -DNC=$(COLORS)
CCFLAGS += -DREPR_$(REPR)

all: $(SOURCES) $(TARGETS)

$(TARGETS): $(OBJECTS)
	@echo "linking $@"
	@$(CC) $(filter-out $(addsuffix .o,$(TARGETS)),$(OBJECTS)) $@.o $(LDFLAGS) -o $@

%.o: %.cpp
	@echo "compiling $<"
	@$(CC) $(CCFLAGS) -c $< -o $@

clean:
	@echo "cleaning up"
	@rm -rf $(OBJECTS) $(TARGETS)
