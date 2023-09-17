VERSION := $(shell grep Version: DESCRIPTION | perl -pe 's/.+: //')
PKG     := $(shell basename `pwd`)
DIR     := $(shell pwd)
MAN     := $(/man/*)
build:
	Rscript bin/rman.R R/canalization.R;
	R CMD build .;
	R CMD check $(PKG)_$(VERSION).tar.gz;
check: 
	R CMD check $(PKG)_$(VERSION).tar.gz

man/%.Rd: R/%.R
	Rscript bin/rman.R $<
