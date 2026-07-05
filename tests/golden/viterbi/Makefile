# Makefile for generating Viterbi golden files
#
# This compiles gen_viterbi.c against the original tRNAscan-SE C library

ORIGINAL_SRC = ../../../original/src
ORIGINAL_LIB = ../../../original/lib
MODEL_DIR = ../../../original/lib/models

CC = gcc
CFLAGS = -I$(ORIGINAL_SRC) -g -Wall
LDFLAGS =

# Source files needed from original (core CM/Viterbi functions)
SRCS = $(ORIGINAL_SRC)/align.c \
       $(ORIGINAL_SRC)/debug.c \
       $(ORIGINAL_SRC)/emit.c \
       $(ORIGINAL_SRC)/fastmodelmaker.c \
       $(ORIGINAL_SRC)/fast-dbviterbi.c \
       $(ORIGINAL_SRC)/konings.c \
       $(ORIGINAL_SRC)/lengthdist.c \
       $(ORIGINAL_SRC)/maxmodelmaker.c \
       $(ORIGINAL_SRC)/misc.c \
       $(ORIGINAL_SRC)/model.c \
       $(ORIGINAL_SRC)/modelmaking.c \
       $(ORIGINAL_SRC)/prior.c \
       $(ORIGINAL_SRC)/probify.c \
       $(ORIGINAL_SRC)/save.c \
       $(ORIGINAL_SRC)/scorestack.c \
       $(ORIGINAL_SRC)/smallviterbi.c \
       $(ORIGINAL_SRC)/structs.c \
       $(ORIGINAL_SRC)/trace.c \
       $(ORIGINAL_SRC)/viterbi.c

# SQUID library files
SQUID_SRCS = $(ORIGINAL_SRC)/sqio.c \
             $(ORIGINAL_SRC)/selex.c \
             $(ORIGINAL_SRC)/sqerror.c \
             $(ORIGINAL_SRC)/types.c \
             $(ORIGINAL_SRC)/sre_ctype.c \
             $(ORIGINAL_SRC)/sre_math.c \
             $(ORIGINAL_SRC)/sre_string.c \
             $(ORIGINAL_SRC)/alignio.c \
             $(ORIGINAL_SRC)/interleaved.c \
             $(ORIGINAL_SRC)/msf.c \
             $(ORIGINAL_SRC)/iupac.c \
             $(ORIGINAL_SRC)/revcomp.c \
             $(ORIGINAL_SRC)/stack.c \
             $(ORIGINAL_SRC)/gnuregex.c

ALL_SRCS = $(SRCS) $(SQUID_SRCS)
OBJS = $(ALL_SRCS:.c=.o)

# Default target
all: gen_viterbi

# Build the generator program
gen_viterbi: gen_viterbi.o $(OBJS)
	$(CC) $(CFLAGS) -o $@ $^ $(LDFLAGS) -lm

gen_viterbi.o: gen_viterbi.c
	$(CC) $(CFLAGS) -c $< -o $@

# Generate golden files using TRNA2.cm model (standard eukaryotic model)
generate: gen_viterbi
	@echo "Generating golden files with TRNA2.cm model..."
	./gen_viterbi $(MODEL_DIR)/TRNA2.cm
	@echo "Golden files generated in current directory"

# Clean
clean:
	rm -f gen_viterbi gen_viterbi.o viterbi_scores.txt traceback.txt

.PHONY: all generate clean
