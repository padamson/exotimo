EXSRC = ../src/PGBF.chpl ../src/CGBF.chpl ../src/basis.chpl ../src/cints.chpl \
      ../src/constants.chpl ../src/crys.chpl ../src/classical.chpl \
      ../src/Ints.chpl ../src/eigen.chpl ../src/hartree_fock.chpl \
      ../src/io.chpl

H2OSRC =  ../src/io.chpl ../src/classical.chpl

CHPL = chpl

TARGETS = exotimo h2o

default: h2o

all: $(TARGETS)

exotimo: exotimo.chpl $(EXSRC)
	$(CHPL) -o ../bin/exotimo exotimo.chpl

h2o: h2o.chpl $(H2OSRC)
	$(CHPL) -compopts '-M ../src' -o ../bin/h2o h2o.chpl

clean:
	rm -f ../bin/exotimo ../bin/h2o

