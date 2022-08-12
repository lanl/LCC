# ------------------------------------------------
# Makefile for LCC program
#
# Author: cnegre@lanl.gov
# Date  : July 2018
#
# ------------------------------------------------

# Compiler and flags:

F90      			= gfortran
F90FLAGS 			= -fopenmp -g
LINKFLAGS 		= -llapack -lblas -fopenmp
MAKELIB       = ON
AR            = /usr/bin/ar cq
RANLIB        = /usr/bin/ranlib
MPI           = OFF

# Links to PRGRESS and BML libraries:

BML_PATH			 =	../bml/install/lib
PROGRESS_PATH	 =	../qmd-progress/install/lib
LIBS 					 = -L$(PROGRESS_PATH) -lprogress -L$(BML_PATH) -lbml_fortran -lbml -llapack -lblas
LIBS 					+= -I$(PROGRESS_PATH)/include/ -L$(PROGRESS_PATH)/lib/ -lprogress 
LIBS 					+= -I$(BML_PATH)/include/ -L$(BML_PATH)/lib/ -lbml_fortran -lbml 
#LIBS 					+= -I$(HOME)/magma/include -L$(HOME)/magma/lib -lmagma 
#LIBS				  += -L${CUDA_CUDART_LIBRARY} -lcublas -lcudart -lcusolver
LIBS 					+= -llapack -lblas
INCLUDES  	 	+= -I$(BML_PATH)/../include -I$(PROGRESS_PATH)/../include
F90FLAGS 			+= -I$(BML_PATH)/../include -I$(PROGRESS_PATH)/../include

# Process preproc flags 
PPFLAGS = -DMPI$(MPI)

# Roules for creating objects:

%.o:%.F90
	$(F90) $(F90FLAGS) $(PPFLAGS) -c $<

# Program executable name:

TARGET = lcc_main

default: $(TARGET)

MAIN = lcc_main.o
OBJECTS = lcc_message_mod.o lcc_parser_mod.o lcc_constants_mod.o
OBJECTS += lcc_structs_mod.o lcc_lattice_mod.o lcc_aux_mod.o 
OBJECTS += lcc_allocation_mod.o lcc_build_mod.o lcc_regular_mod.o
OBJECTS += lcc_check_mod.o lcc_lib.o lcc_radial_mod.o 
OBJECTS += lcc_string_mod.o 

# Dependencies:

lcc_main.o : lcc_message_mod.o lcc_parser_mod.o lcc_constants_mod.o lcc_structs_mod.o lcc_lattice_mod.o lcc_aux_mod.o lcc_build_mod.o lcc_regular_mod.o lcc_message_mod.o lcc_radial_mod.o lcc_string_mod.o
lcc_lib.o : lcc_message_mod.o lcc_parser_mod.o lcc_constants_mod.o lcc_structs_mod.o lcc_lattice_mod.o lcc_aux_mod.o lcc_build_mod.o lcc_regular_mod.o lcc_message_mod.o lcc_radial_mod.o lcc_string_mod.o

lcc_parser_mod.o : lcc_structs_mod.o lcc_message_mod.o lcc_aux_mod.o
lcc_structs_mod.o : lcc_constants_mod.o	
lcc_lattice_mod.o : lcc_structs_mod.o lcc_message_mod.o
lcc_lattice_mod.o : lcc_allocation_mod.o
lcc_aux_mod.o : lcc_allocation_mod.o  
lcc_build_mod.o : lcc_allocation_mod.o
lcc_regular_mod.o : lcc_allocation_mod.o
lcc_message_mod.o : lcc_constants_mod.o
lcc_check_mod.o : lcc_message_mod.o	
lcc_radial_mod.o : lcc_constants_mod.o lcc_string_mod.o

$(TARGET):	$(OBJECTS) $(MAIN)
	(mkdir build)
	$(F90) -o ./build/$(TARGET) $(LINKFLAGS) $(OBJECTS) $(MAIN) $(LIBS)
	$(AR) ./build/liblcc.a lcc_lib.o  $(OBJECTS)
	$(RANLIB) ./build/liblcc.a

test: 
	(./run_test.sh)

clean:
	rm *.o *.mod main



