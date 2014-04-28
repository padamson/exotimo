SRC = ../src/PGBF.chpl ../src/CGBF.chpl ../src/basis.chpl ../src/cints.chpl \
      ../src/constants.chpl ../src/crys.chpl ../src/classical.chpl \
      ../src/Ints.chpl ../src/eigen.chpl ../src/hartree_fock.chpl \
      ../src/io.chpl

CHPL = chpl

TARGETS = exotimo h2o

default: h2o

all: $(TARGETS)

exotimo: exotimo.chpl $(SRC)
	$(CHPL) -o $@ exotimo.chpl $(SRC)

h2o: ../src/io.chpl h2o.chpl
	$(CHPL) -o ~/Research/bin/h2o ../src/io.chpl h2o.chpl

clean:
	rm -f $(TARGETS)

