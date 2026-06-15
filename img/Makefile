#=========
include ../gosh/Mkinclude
#=========
TARGETS=calc_stat

# Suffix
.SUFFIXES : .f90
.f90:
	$(FC) ${FFLAGS} $^ -o $@

# commands
all: $(TARGETS)
# clean
clean:
	rm -rf  $(TARGETS)
