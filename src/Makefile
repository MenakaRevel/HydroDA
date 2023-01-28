#=========
include ../gosh/Mkinclude
#=========
#*** module
OBJECTS=\
letkf.o\
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
letkf.o  : letkf.f90 common.o
common.o : common.f90

# commands
all: $(TARGET)
# clean
.PHONY: clean
clean:
	rm -rf *.o *~ *.mod $(TARGET)
