##################################################
FC = ifort #gfortran
FFLAGS=-assume byterecl -heap-arrays -nogen-interfaces -mcmodel=large -fPIC -shared-intel -g -traceback -parallel

.SUFFIXES : .o .f90 .F
.f90.o:
	$(FC) -c ${FFLAGS} ${INC} $(MODS) $*.f90 -o $*.o
src: $(OBJECTS)
	ar -rv srcs.a $(OBJECTS)
##################################################
# simulation  
##################################################
TARGET1 = print_grid
TARGET2 = txt_vector
TARGET3 = print_rivvec
TARGET4 = print_fldpth

##################################################
$(TARGET1) : $(TARGET1).o $(LIBS)
	$(FC) $(LFLAGS) $@.o $(LIBS) $(MOD) -o $@ $(LIB)

$(TARGET2) : $(TARGET2).o $(LIBS)
	$(FC) $(LFLAGS) $@.o $(LIBS) $(MOD) -o $@ $(LIB)

$(TARGET3) : $(TARGET3).o $(LIBS)
	$(FC) $(LFLAGS) $@.o $(LIBS) $(MOD) -o $@ $(LIB)

$(TARGET4) : $(TARGET4).o $(LIBS)
	$(FC) $(LFLAGS) $@.o $(LIBS) $(MOD) -o $@ $(LIB)

##################################################
# Commands
##################################################
all : ${TARGET1} ${TARGET2} ${TARGET3} ${TARGET4}
clean: 
	$(RM) -f *.o *.s core *~ *trace *.mod

