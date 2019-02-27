################################################################################
# LIO MAKEFILE
################################################################################

all: liosolo liblio g2g
again:
	make cuda=0 intel=0 precision=1 libxc=1 libint=1 -j


.PHONY: liosolo
liosolo: liblio
	$(MAKE) -C liosolo


.PHONY: liblio
liblio: g2g
	$(MAKE) -C lioamber


.PHONY: g2g
g2g:
	$(MAKE) -C g2g


.PHONY: clean
clean:
	$(MAKE) clean -C liosolo
	$(MAKE) clean -C lioamber
	$(MAKE) clean -C g2g

################################################################################
