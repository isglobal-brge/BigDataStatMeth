
R ?= R

.PHONY: all
all:
    $(MAKE) clean
$(MAKE) build
$(MAKE) install
$(MAKE) test
$(MAKE) check

.PHONY: some
some:
    $(MAKE) clean
$(MAKE) build
$(MAKE) install
$(MAKE) test

.PHONY: clean
clean:
    $(RM) BigDataStatMeth_0.99.16.tar.gz
$(RM) src/*.o
$(RM) src/*.so

.PHONY: build
build:
    $(R) CMD build . --no-build-vignettes

.PHONY: install
install:
    $(R) CMD INSTALL BigDataStatMeth_0.99.16.tar.gz

.PHONY: uninstall
uninstall:
    $(R) CMD REMOVE BigDataStatMeth || true

.PHONY: test
test:
    $(R) -e 'require(BigDataStatMeth); test.BigDataStatMeth()'

.PHONY: check
check:
    _R_CHECK_CRAN_INCOMING_REMOTE_=false $(R) CMD check BigDataStatMeth_0.99.16.tar.gz --as-cran --ignore-vignettes --no-stop-on-test-error

.PHONY: revision
revision:
    echo "Revision: $(shell git rev-parse HEAD)" >> DESCRIPTION
