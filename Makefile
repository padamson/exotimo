
MAKEFLAGS = --no-print-directory

CHPL = chpl

TARGETS = \
	h2o \
	exotimo \

default: all

all: $(TARGETS)

clean: FORCE
	rm -f $(TARGETS)

h2o: h2o.chpl molecule.chpl
	$(CHPL) -o $@ $<

FORCE:
