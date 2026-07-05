/* gen_eufind.c - Generate golden files for EuFindtRNA phase
 *
 * Tests the Pavesi algorithm functions for finding eukaryotic tRNA
 * transcriptional control regions (A-box, B-box) and termination signals.
 * Outputs intermediate scoring values to golden files for Rust reimplementation.
 *
 * Compile:
 *   gcc -I../../original/squid -I../../original/src \
 *       gen_eufind.c \
 *       ../../original/src/pavesi.o \
 *       ../../original/squid/libsquid.a \
 *       -o gen_eufind -lm
 */

#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include "squid.h"
#include "eufind_const.h"
#include "pavesi.h"

/* Test sequence: tRNA-Phe from C. elegans Example1.fa around line 254-256 */
const char *test_trna_seq =
    "TCGCTGTTAGTTACCATCGCACGGATGGCCGAGTGGTCTAAGGCGCCAGA"
    "CTCAAGCGAAATGCTTGCCTCATGCTCGAGGTCGACTGGGTGTTCTGGTA"
    "CTCGTATGGGTGCGTGGGTTCGAATCCCACTTCGTGCAG";

/* Another test sequence with B-box pattern */
const char *test_bbox_seq =
    "GGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGG"
    "GCGGTTTTGGGGGTGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGG"
    "GGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGG";

void test_bbox_scanning(FILE *fp_bbox);
void test_abox_detection(FILE *fp_abox);
void test_trna_full_detection(FILE *fp_trna);
void print_trna_info(FILE *fp, TRNA_TYPE *tRNA, const char *label);

int main(int argc, char **argv) {
    FILE *fp_bbox, *fp_abox, *fp_trna;

    printf("Generating EuFindtRNA golden files...\n");

    /* Open output files */
    fp_bbox = fopen("bbox_scores.txt", "w");
    if (!fp_bbox) {
        fprintf(stderr, "Error: Cannot open bbox_scores.txt\n");
        return 1;
    }

    fp_abox = fopen("abox_scores.txt", "w");
    if (!fp_abox) {
        fprintf(stderr, "Error: Cannot open abox_scores.txt\n");
        fclose(fp_bbox);
        return 1;
    }

    fp_trna = fopen("trna_detection.txt", "w");
    if (!fp_trna) {
        fprintf(stderr, "Error: Cannot open trna_detection.txt\n");
        fclose(fp_bbox);
        fclose(fp_abox);
        return 1;
    }

    /* Write headers */
    fprintf(fp_bbox, "# B-box Score Detection Golden File\n");
    fprintf(fp_bbox, "# Format: position score\n");
    fprintf(fp_bbox, "# Cutoff: %.2f\n\n", BBOX_CUTOFF);

    fprintf(fp_abox, "# A-box Detection Golden File\n");
    fprintf(fp_abox, "# Format: position Abox_score ABdist ABdist_score total_score\n\n");

    fprintf(fp_trna, "# Full tRNA Detection Golden File\n");
    fprintf(fp_trna, "# Shows complete TRNA_TYPE structure after detection\n\n");

    /* Run tests */
    test_bbox_scanning(fp_bbox);
    test_abox_detection(fp_abox);
    test_trna_full_detection(fp_trna);

    /* Close files */
    fclose(fp_bbox);
    fclose(fp_abox);
    fclose(fp_trna);

    printf("Golden files generated successfully:\n");
    printf("  - bbox_scores.txt\n");
    printf("  - abox_scores.txt\n");
    printf("  - trna_detection.txt\n");

    return 0;
}

void test_bbox_scanning(FILE *fp_bbox) {
    char intseq[200];
    float score;
    int seqidx, found;
    int seqlen = strlen(test_trna_seq);

    fprintf(fp_bbox, "## Test 1: B-box scanning on tRNA-Phe sequence\n");
    fprintf(fp_bbox, "Sequence length: %d\n", seqlen);
    fprintf(fp_bbox, "Sequence: %s\n\n", test_trna_seq);

    /* Integer encode the sequence */
    IntEncodeSeq(intseq, (char*)test_trna_seq, seqlen);

    /* Scan for B-boxes */
    seqidx = BBOX_START_IDX - 1;
    fprintf(fp_bbox, "Position  Score     Status\n");
    fprintf(fp_bbox, "--------  --------  ------\n");

    while (seqidx < seqlen - BBOX_LEN) {
        int start_idx = seqidx;
        found = GetBbox(&score, &seqidx, intseq, seqlen, 0, 0);

        if (found) {
            fprintf(fp_bbox, "%8d  %8.4f  FOUND (end=%d)\n",
                   seqidx, score, seqidx + BBOX_LEN);
            /* Continue searching after this B-box */
        } else {
            /* No more B-boxes found */
            break;
        }
    }

    fprintf(fp_bbox, "\n## Test 2: B-box scanning on artificial sequence with strong pattern\n");
    seqlen = strlen(test_bbox_seq);
    fprintf(fp_bbox, "Sequence length: %d\n\n", seqlen);

    IntEncodeSeq(intseq, (char*)test_bbox_seq, seqlen);

    seqidx = BBOX_START_IDX - 1;
    fprintf(fp_bbox, "Position  Score     Status\n");
    fprintf(fp_bbox, "--------  --------  ------\n");

    /* Show first few positions */
    for (int i = 0; i < 5 && seqidx < seqlen - BBOX_LEN; i++) {
        found = GetBbox(&score, &seqidx, intseq, seqlen, 0, 0);
        if (found) {
            fprintf(fp_bbox, "%8d  %8.4f  FOUND\n", seqidx, score);
        } else {
            break;
        }
    }

    fprintf(fp_bbox, "\n");
}

void test_abox_detection(FILE *fp_abox) {
    TRNA_TYPE tRNA;
    char intseq[200];
    int seqlen = strlen(test_trna_seq);

    fprintf(fp_abox, "## Test: A-box detection for tRNA-Phe sequence\n");
    fprintf(fp_abox, "Sequence: %s\n\n", test_trna_seq);

    /* Initialize tRNA structure */
    Init_tRNA(&tRNA);

    /* Integer encode sequence */
    IntEncodeSeq(intseq, (char*)test_trna_seq, seqlen);

    /* First find a B-box */
    float score;
    int seqidx = BBOX_START_IDX - 1;
    if (GetBbox(&score, &seqidx, intseq, seqlen, 0, 0)) {
        tRNA.Bbox_st = seqidx;
        tRNA.Bbox_end = seqidx + BBOX_LEN - 1;
        tRNA.BboxSc = score;

        fprintf(fp_abox, "B-box found at position %d (score=%.4f)\n",
               tRNA.Bbox_st, tRNA.BboxSc);
        fprintf(fp_abox, "B-box end: %d\n\n", tRNA.Bbox_end);

        /* Now search for best A-box */
        fprintf(fp_abox, "Searching for A-box (AB_BOX_DIST_RANGE=%d):\n\n",
               AB_BOX_DIST_RANGE);

        GetBestABox(&tRNA, (char*)test_trna_seq, intseq, seqlen, 0, 0,
                   MIN_AB_BOX_DIST + AB_BOX_DIST_RANGE, 0);

        fprintf(fp_abox, "Best A-box results:\n");
        fprintf(fp_abox, "  Position:      %d - %d\n", tRNA.Abox_st, tRNA.Abox_end);
        fprintf(fp_abox, "  A-box score:   %.4f\n", tRNA.AboxSc);
        fprintf(fp_abox, "  AB distance:   %d\n", tRNA.Bbox_st - tRNA.Abox_end - 1);
        fprintf(fp_abox, "  AB dist score: %.4f\n", tRNA.ABdistSc);
        fprintf(fp_abox, "  Combined:      %.4f\n",
               tRNA.AboxSc + tRNA.BboxSc + tRNA.ABdistSc);
    } else {
        fprintf(fp_abox, "No B-box found in test sequence\n");
    }

    fprintf(fp_abox, "\n");
}

void test_trna_full_detection(FILE *fp_trna) {
    TRNA_TYPE tRNA;
    char intseq[200];
    int seqlen = strlen(test_trna_seq);
    float score;
    int seqidx;

    fprintf(fp_trna, "## Full tRNA Detection Test\n");
    fprintf(fp_trna, "Sequence: %s\n", test_trna_seq);
    fprintf(fp_trna, "Length: %d\n\n", seqlen);

    /* Initialize */
    Init_tRNA(&tRNA);
    IntEncodeSeq(intseq, (char*)test_trna_seq, seqlen);

    /* Step 1: Find B-box */
    fprintf(fp_trna, "=== Step 1: B-box Detection ===\n");
    seqidx = BBOX_START_IDX - 1;
    if (GetBbox(&score, &seqidx, intseq, seqlen, 0, 0)) {
        tRNA.Bbox_st = seqidx;
        tRNA.Bbox_end = seqidx + BBOX_LEN - 1;
        tRNA.BboxSc = score;
        fprintf(fp_trna, "B-box found: pos=%d score=%.4f\n\n",
               tRNA.Bbox_st, tRNA.BboxSc);
    } else {
        fprintf(fp_trna, "No B-box found\n\n");
        return;
    }

    /* Step 2: Find best A-box */
    fprintf(fp_trna, "=== Step 2: A-box Detection ===\n");
    GetBestABox(&tRNA, (char*)test_trna_seq, intseq, seqlen, 0, 0,
               MIN_AB_BOX_DIST + AB_BOX_DIST_RANGE, 0);
    fprintf(fp_trna, "A-box found: pos=%d-%d score=%.4f ABdist=%d ABscore=%.4f\n\n",
           tRNA.Abox_st, tRNA.Abox_end, tRNA.AboxSc,
           tRNA.Bbox_st - tRNA.Abox_end - 1, tRNA.ABdistSc);

    /* Step 3: Find termination signal */
    fprintf(fp_trna, "=== Step 3: Termination Signal ===\n");
    int term_found = GetBestTrxTerm(&tRNA, (char*)test_trna_seq, seqlen, 0.0);
    fprintf(fp_trna, "Term signal: found=%d pos=%d score=%.4f\n\n",
           term_found, tRNA.Term_st, tRNA.TermSc);

    /* Step 4: Calculate total score */
    fprintf(fp_trna, "=== Step 4: Total Score ===\n");
    tRNA.totSc = tRNA.AboxSc + tRNA.BboxSc + tRNA.ABdistSc + tRNA.TermSc;
    fprintf(fp_trna, "Total score: %.4f\n", tRNA.totSc);
    fprintf(fp_trna, "  Components: A=%.4f B=%.4f ABdist=%.4f Term=%.4f\n\n",
           tRNA.AboxSc, tRNA.BboxSc, tRNA.ABdistSc, tRNA.TermSc);

    /* Step 5: Determine start/end positions */
    if (tRNA.Abox_st > 0) {
        tRNA.start = tRNA.Abox_st - 5;  /* 5bp upstream of A-box */
    } else {
        tRNA.start = 0;
    }

    if (tRNA.Term_st > 0) {
        tRNA.end = tRNA.Term_st + 4;  /* End of TTTT terminator */
    } else {
        tRNA.end = seqlen - 1;
    }

    /* Print complete structure */
    fprintf(fp_trna, "=== Complete TRNA_TYPE Structure ===\n");
    print_trna_info(fp_trna, &tRNA, "Final");
}

void print_trna_info(FILE *fp, TRNA_TYPE *tRNA, const char *label) {
    fprintf(fp, "%s tRNA info:\n", label);
    fprintf(fp, "  iso_type:  %s\n", tRNA->iso_type);
    fprintf(fp, "  acodon:    %s\n", tRNA->acodon);
    fprintf(fp, "  start:     %d\n", tRNA->start);
    fprintf(fp, "  end:       %d\n", tRNA->end);
    fprintf(fp, "  Abox_st:   %d\n", tRNA->Abox_st);
    fprintf(fp, "  Abox_end:  %d\n", tRNA->Abox_end);
    fprintf(fp, "  Abox_gap:  %d\n", tRNA->Abox_gap);
    fprintf(fp, "  Bbox_st:   %d\n", tRNA->Bbox_st);
    fprintf(fp, "  Bbox_end:  %d\n", tRNA->Bbox_end);
    fprintf(fp, "  Term_st:   %d\n", tRNA->Term_st);
    fprintf(fp, "  totSc:     %.6f\n", tRNA->totSc);
    fprintf(fp, "  AboxSc:    %.6f\n", tRNA->AboxSc);
    fprintf(fp, "  BboxSc:    %.6f\n", tRNA->BboxSc);
    fprintf(fp, "  ABdistSc:  %.6f\n", tRNA->ABdistSc);
    fprintf(fp, "  TermSc:    %.6f\n", tRNA->TermSc);
    fprintf(fp, "  intron:    %d\n", tRNA->intron);
    fprintf(fp, "  acodon_idx:%d\n", tRNA->acodon_idx);
    fprintf(fp, "  idno:      %d\n", tRNA->idno);
    fprintf(fp, "\n");
}
