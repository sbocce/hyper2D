############################################################
#                         MAKEFILE                         #
############################################################

CMPF  = gfortran -c
LNK   = gfortran 

OPTS = -O3

OBJS = global_module.o  tools.o  pde.o  integration.o 

######### Executable generation by the linker

hyper2D.exe:     hyper2D.o $(OBJS) 
	$(LNK) $(OPTS) hyper2D.o $(OBJS) -o hyper2D.exe 

######### Objects generation

hyper2D.o:       hyper2D.f03  $(OBJS)
	$(CMPF) $(OPTS) hyper2D.f03

integration.o:   integration.f03 global_module.o pde.o 
	$(CMPF) $(OPTS) integration.f03

pde.o:           pde.f03 global_module.o
	$(CMPF) $(OPTS) pde.f03

tools.o:         tools.f03 pde.o global_module.o
	$(CMPF) $(OPTS) tools.f03

global_module.o: global_module.f03
	$(CMPF) $(OPTS) global_module.f03

########## Cleaning commands

# Use "clean" to remove the objects and exectutable files
clean: 
	@echo cleaning objects, modules and executables 
	rm  -f  *.o  *.mod  *.exe  *~

# Use "cleanoutput" to remove the output files, located in the dumps directory
cleanoutput:
	@echo cleaning output and dump files
	rm  -f  dumps/*

