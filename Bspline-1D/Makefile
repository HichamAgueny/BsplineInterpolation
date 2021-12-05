FC = nvfortran
FCFLAGS = -O3

EXEC = bspline_test-1d.exe

OBJECTS = bspline_kinds_module.o bspline_sub_module.o \
          bspline_module.o bspline_test.o
MODULES = bspline_kinds_module.mod bspline_sub_module.mod \
          bspline_module.mod

$(EXEC): $(OBJECTS) $(MODULES)
	$(FC) $(FCFLAGS) $(OBJECTS) -o $(EXEC)
##
%.o: %.f90 
	$(FC) $(FCFLAGS) -c $<
%.mod: %.f90
	$(FC) $(FCFLAGS) -c $<

#Cleaning
.PHONY: clean
clean: 
	rm -f $(OBJECTS) 
	rm -f $(MODULES) 
	rm -f $(EXEC)
