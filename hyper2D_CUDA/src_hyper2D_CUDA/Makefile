############################################################
#                         MAKEFILE                         #
############################################################

CMPF  = nvfortran -c
LNK   = nvfortran 

# ---- Compile options. cc35 stands for "compute capability 3.5"
# OPTS = -O3 -cuda -gpu=cc35 -r4 # float
OPTS = -O3 -cuda -gpu=cc35 -r8 # double precision
# OPTS = -cuda -gpu=cc35,debug -r4 -g # For debugging, use these options

OBJS = global_module.o  tools.o  pde.o  integration.o

######### Executable generation by the linker

hyper2D.exe:     hyper2D.o $(OBJS) 
	$(LNK) $(OPTS) hyper2D.o $(OBJS) -o hyper2D.exe 

######### Objects generation

hyper2D.o:       hyper2D.cuf  $(OBJS)
	$(CMPF) $(OPTS) hyper2D.cuf

integration.o:   integration.cuf global_module.o pde.o 
	$(CMPF) $(OPTS) integration.cuf

pde.o:           pde.cuf global_module.o
	$(CMPF) $(OPTS) pde.cuf

tools.o:         tools.cuf pde.o global_module.o
	$(CMPF) $(OPTS) tools.cuf

global_module.o: global_module.cuf
	$(CMPF) $(OPTS) global_module.cuf

########## Cleaning commands

# Use "clean" to remove the objects and exectutable files
clean: 
	@echo cleaning objects, modules and executables 
	rm  -f  *.o  *.mod  *.exe  *~

# Use "cleanoutput" to remove the output files, located in the dumps directory
cleanoutput:
	@echo cleaning output and dump files
	rm  -f  dumps/*

