MAKEFLAGS = --no-print-directory

CHPL = chpl

TARGETS = \
	testbench \

default: all

all: $(TARGETS)

clean: FORCE
	rm -f $(TARGETS)

testbench: testbench.chpl io.chpl constants.chpl
	$(CHPL) -o ~/Research/bin/testbench testbench.chpl io.chpl constants.chpl

FORCE:
