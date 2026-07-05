# Makefile for EuFindtRNA golden file generation

SRCDIR = ../../../original/src
OBJS = $(SRCDIR)/pavesi.o \
       $(SRCDIR)/revcomp.o \
       $(SRCDIR)/iupac.o \
       $(SRCDIR)/sqerror.o \
       $(SRCDIR)/sre_ctype.o \
       $(SRCDIR)/sre_string.o \
       $(SRCDIR)/sre_math.o

CFLAGS = -I$(SRCDIR)
LDFLAGS = -lm

all: gen_eufind golden_files

gen_eufind: gen_eufind.c $(OBJS)
	gcc $(CFLAGS) -o gen_eufind gen_eufind.c $(OBJS) $(LDFLAGS)

golden_files: gen_eufind
	./gen_eufind

clean:
	rm -f gen_eufind bbox_scores.txt abox_scores.txt trna_detection.txt

.PHONY: all golden_files clean
