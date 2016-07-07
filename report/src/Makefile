ROOT_DIR:=$(shell dirname $(realpath $(lastword $(MAKEFILE_LIST))))

SRCDIR = $(ROOT_DIR)
OUTDIR = $(realpath $(ROOT_DIR)/..)
AUXDIR = $(OUTDIR)/aux

SOURCE = $(SRCDIR)/main.tex
TARGET = $(OUTDIR)/jones_thesis.pdf

SUBDIRS := $(sort $(dir $(wildcard $(SRCDIR)/*/)))
AUXSUBDIRS := $(SUBDIRS:$(SRCDIR)/%=$(AUXDIR)/%)

DEPEXTS:=png pdf jpg jpeg tiff bmp
BRACKETS:=)|(

DEPEXTENSIONS := ($(subst $(eval) ,$(BRACKETS),$(DEPEXTS)))

DEPFILES := $(shell find $(SRCDIR) -regextype posix-extended -regex "$(SRCDIR)/.*($(DEPEXTENSIONS))")
AUXDEPFILES := $(subst $(SRCDIR),$(AUXDIR),$(DEPFILES))

all: $(TARGET) 

$(TARGET): $(DEPFILES) | $(AUXSUBDIRS)
	@echo "Making target"
	@latexmk -outdir=$(AUXDIR) $(notdir $(SOURCE))
	cp "$(SOURCE:$(SRCDIR)/%.tex=$(AUXDIR)/%.pdf)" "$(TARGET)"

$(AUXSUBDIRS):
	mkdir -p $@

$(DEPFILES): | $(AUXSUBDIRS)
	cp "$@" "$(subst $(SRCDIR),$(AUXDIR),$@)"

clean:
	rm -rf $(AUXDIR) $(TARGET)
