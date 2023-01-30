#=========
include ../gosh/Mkinclude
#=========
#*** module
OBJECTS=\
mod_varxf.o\
mod_patch.o\
mod_obser.o\
mod_letkf.o\
common.o

#=========
TARGET=\
data_assim

# Suffix
.SUFFIXES: .o .f90

%.o: %.f90
	$(FC) ${FFLAGS} ${INC} ${LIB} -c $<

$(TARGET): $(OBJECTS) $(TARGET).o
	$(FC) $(FFLAGS) $@.o $(OBJECTS) -o $@ $(LIB)

# Dependecies
mod_varxf.o  : mod_varxf.f90 common.o
mod_patch.o  : mod_patch.f90 common.o
mod_obser.o  : mod_obser.f90 common.o
mod_letkf.o  : mod_letkf.f90 common.o
common.o     : common.f90

# commands
all: $(TARGET)
# clean
.PHONY: clean
clean:
	rm -rf *.o *~ *.mod $(TARGET)
