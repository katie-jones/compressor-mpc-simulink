ROOT_DIR:=$(shell dirname $(realpath $(lastword $(MAKEFILE_LIST))))

SRCDIR = $(ROOT_DIR)
OUTDIR = $(realpath $(ROOT_DIR)/..)
AUXDIR = $(OUTDIR)/aux

SOURCE = $(SRCDIR)/main.tex
TARGET = $(OUTDIR)/jones_thesis.pdf

PDF_OUTPUT = $(SOURCE:$(SRCDIR)/%.tex=$(AUXDIR)/%.pdf)

TEXFILES = $(shell find $(SRCDIR) -name "*.tex")
GLOSSARIESFILES = $(wildcard $(AUXDIR)/*glo)

SUBDIRS := $(sort $(dir $(wildcard $(SRCDIR)/*/)))
AUXSUBDIRS := $(SUBDIRS:$(SRCDIR)/%=$(AUXDIR)/%)

DEPEXTS:=png pdf jpg jpeg tiff bmp bib
BRACKETS:=)|(

DEPEXTENSIONS := ($(subst $(eval) ,$(BRACKETS),$(DEPEXTS)))

IMGFILES := $(shell find $(SRCDIR) -regextype posix-extended -regex "$(SRCDIR)/.*($(DEPEXTENSIONS))")
AUXIMGFILES := $(subst $(SRCDIR),$(AUXDIR),$(IMGFILES))

all: $(TARGET) 

$(TARGET): $(PDF_OUTPUT)
	cp "$(PDF_OUTPUT)" "$(TARGET)"

$(PDF_OUTPUT): $(IMGFILES) $(TEXFILES) | $(AUXSUBDIRS)
	@latexmk -bibtex -outdir=$(AUXDIR) $(notdir $(SOURCE))

$(AUXSUBDIRS):
	mkdir -p $@

$(IMGFILES): | $(AUXSUBDIRS)
	cp "$@" "$(subst $(SRCDIR),$(AUXDIR),$@)"

glossary: firstrun | $(AUXSUBDIRS)
	$(eval GLOSSARIESFILES = $(wildcard $(AUXDIR)/*glo))
	$(foreach x,$(GLOSSARIESFILES), xindy -L english -C utf8 -I xindy -M $(basename $(x)) -t $(x:%glo=%glg) -o $(x:%glo=%gls) $(x); )
	@rm -f $(PDF_OUTPUT)

firstrun: | $(AUXSUBDIRS)
	pdflatex --shell-escape  -recorder -output-directory="$(AUXDIR)"  "$(SOURCE)"

clean:
	rm -rf $(PDF_OUTPUT) $(TARGET)

cleaner: clean
	rm -rf $(AUXDIR)

.PHONY: all clean cleaner glossary firstrun
