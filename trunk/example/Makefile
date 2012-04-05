SRC = ../src/PGBF.chpl ../src/CGBF.chpl ../src/basis.chpl ../src/cints.chpl \
      ../src/constants.chpl ../src/crys.chpl ../src/classical.chpl \
      ../src/Ints.chpl ../src/eigen.chpl ../src/hartree_fock.chpl \
      ../src/io.chpl

CHPL = chpl

TARGETS = exotimo

default: all

all: $(TARGETS)

exotimo: exotimo.chpl $(SRC)
	$(CHPL) -o $@ exotimo.chpl $(SRC)

clean:
	rm -f $(TARGETS)

