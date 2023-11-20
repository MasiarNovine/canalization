VERSION := $(shell grep Version: DESCRIPTION | perl -pe 's/.+: //')
PKG     := $(shell basename `pwd`)
DIR     := $(shell pwd)
MAN     := $(/man/*)
build:
	rm man/*
	Rscript bin/rman.R R/canalization.R;
	R CMD build .;
	R CMD check $(PKG)_$(VERSION).tar.gz;
	R CMD INSTALL $(PKG)_$(VERSION).tar.gz;
