#=========
include ../gosh/Mkinclude
#=========
TARGETS=\
data_assim\
make_restart

# Suffix
.SUFFIXES : .f90
.f90:
	$(FC) ${FFLAGS} ${INC} ${LIB} $^ -o $@

# commands
all: $(TARGETS)
# clean
clean:
	rm -rf  $(TARGETS)
