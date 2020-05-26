################################################################################
# LIO MAKEFILE ################################################################################

all: liosolo liblio g2g tools

again:
	make cuda=2 sm61=1 intel=0 precision=1 libxc=2 libint=1

.PHONY: liosolo
liosolo: liblio
	$(MAKE) -C liosolo


.PHONY: liblio
liblio: g2g
	$(MAKE) -C lioamber


.PHONY: g2g
g2g:
	$(MAKE) -C g2g

.PHONY: tools
tools:
	$(MAKE) -C tools

.PHONY: check
check:
	$(MAKE) check -C test

.PHONY: compile
compile:
	$(MAKE) compile -C test

.PHONY: clean
clean:
	$(MAKE) clean -C liosolo
	$(MAKE) clean -C lioamber
	$(MAKE) clean -C g2g
	$(MAKE) clean -C test
	$(MAKE) clean -C tools

################################################################################
