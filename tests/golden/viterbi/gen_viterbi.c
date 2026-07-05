/* gen_viterbi.c
 * Generate golden files for Viterbi alignment testing
 *
 * This program:
 * 1. Loads a CM model
 * 2. Converts it to integer model with RearrangeCM()
 * 3. Runs ViterbiAlign() on test sequences
 * 4. Outputs scores and traceback trees
 */

#include <stdio.h>
#include <stdlib.h>
#include <string.h>

#include "funcs.h"
#include "structs.h"
#include "squid.h"

/* Test sequences - short tRNA-like sequences */
static char *test_sequences[] = {
    /* Test 1: Simple 73bp tRNA (Phe from Example1) */
    "GCCTCGATAGCTCAGTTGGGAGAGCGTACGACTGAAGATCGTAAGGTCACCAGTTCGATCCTGGTTCGGGGCA",

    /* Test 2: 82bp tRNA (Ser from Example1) */
    "GCAGTCATGTCCGAGTGGTAAGGAGATTGACTAGAAATCAATTGGGCTCTGCCCGCGTAGGTTCGAATCCTGCTGACTGCG",

    /* Test 3: Minimal test case (50bp fragment) */
    "GCCTCGATAGCTCAGTTGGGAGAGCGTACGACTGAAGATCGTAAGGTCAC",

    NULL
};

static char *test_names[] = {
    "Phe-73bp",
    "Ser-82bp",
    "Fragment-50bp",
    NULL
};

/* Helper function to print traceback tree recursively */
void print_traceback_recursive(FILE *fp, struct trace_s *tr, int depth)
{
    int i;

    if (tr == NULL)
        return;

    /* Print indentation */
    for (i = 0; i < depth; i++)
        fprintf(fp, "  ");

    /* Print node info: nodeidx type emitl emitr */
    fprintf(fp, "node=%d type=%d emitl=%d emitr=%d (%s)\n",
            tr->nodeidx, tr->type, tr->emitl, tr->emitr,
            UstatetypeName(tr->type));

    /* Recursively print children */
    if (tr->nxtl != NULL)
        print_traceback_recursive(fp, tr->nxtl, depth + 1);
    if (tr->nxtr != NULL)
        print_traceback_recursive(fp, tr->nxtr, depth + 1);
}

/* Print traceback tree in a structured format */
void print_traceback_tree(FILE *fp, struct trace_s *tr)
{
    fprintf(fp, "Traceback tree (depth-first traversal):\n");
    fprintf(fp, "Format: node=<nodeidx> type=<type> emitl=<emitl> emitr=<emitr> (<statetype_name>)\n");
    fprintf(fp, "---\n");

    if (tr != NULL && tr->nxtl != NULL) {
        print_traceback_recursive(fp, tr->nxtl, 0);
    }

    fprintf(fp, "---\n\n");
}

int main(int argc, char **argv)
{
    struct cm_s      *cm;
    struct istate_s  *icm;
    int               statenum;
    char             *modelfile;
    double            rfreq[ALPHASIZE];
    int               i;
    FILE             *score_fp;
    FILE             *trace_fp;

    /* Check arguments */
    if (argc != 2) {
        fprintf(stderr, "Usage: %s <cmfile>\n", argv[0]);
        fprintf(stderr, "Example: %s ../../original/lib/models/TRNA2.cm\n", argv[0]);
        return 1;
    }

    modelfile = argv[1];

    /* Load CM model */
    printf("Loading CM model from %s...\n", modelfile);
    if (!ReadCM(modelfile, &cm)) {
        fprintf(stderr, "Failed to read CM from %s\n", modelfile);
        return 1;
    }
    printf("Model loaded: %d nodes\n", cm->nodes);

    /* Set uniform background frequencies for RNA */
    for (i = 0; i < ALPHASIZE; i++)
        rfreq[i] = 0.25;

    /* Convert to integer model */
    printf("Converting to integer model...\n");
    if (!RearrangeCM(cm, rfreq, &icm, &statenum)) {
        fprintf(stderr, "Failed to rearrange CM\n");
        FreeCM(cm);
        return 1;
    }
    printf("Integer model created: %d states\n", statenum);

    /* Open output files */
    score_fp = fopen("viterbi_scores.txt", "w");
    trace_fp = fopen("traceback.txt", "w");

    if (!score_fp || !trace_fp) {
        fprintf(stderr, "Failed to open output files\n");
        free(icm);
        FreeCM(cm);
        return 1;
    }

    /* Write headers */
    fprintf(score_fp, "# Viterbi alignment scores\n");
    fprintf(score_fp, "# Format: <test_name> <sequence_length> <score>\n");
    fprintf(score_fp, "# Model: %s\n", modelfile);
    fprintf(score_fp, "# States: %d\n\n", statenum);

    fprintf(trace_fp, "# Viterbi traceback trees\n");
    fprintf(trace_fp, "# Model: %s\n", modelfile);
    fprintf(trace_fp, "# States: %d\n\n", statenum);

    /* Process each test sequence */
    for (i = 0; test_sequences[i] != NULL; i++) {
        char            *seq = test_sequences[i];
        char            *prepseq;
        double           score;
        struct trace_s  *trace;
        int              seqlen;

        seqlen = strlen(seq);
        printf("\nProcessing %s (length=%d)...\n", test_names[i], seqlen);

        /* Prepare sequence (convert to internal format) */
        prepseq = (char *) malloc((seqlen + 1) * sizeof(char));
        strcpy(prepseq, seq);
        PrepareSequence(prepseq);

        /* Run Viterbi alignment */
        if (!ViterbiAlign(icm, statenum, prepseq, &score, &trace)) {
            fprintf(stderr, "ViterbiAlign failed for %s\n", test_names[i]);
            free(prepseq);
            continue;
        }

        printf("  Score: %.4f\n", score);

        /* Write score to file */
        fprintf(score_fp, "%s\t%d\t%.4f\n", test_names[i], seqlen, score);

        /* Write traceback to file */
        fprintf(trace_fp, "=== %s (length=%d, score=%.4f) ===\n",
                test_names[i], seqlen, score);
        print_traceback_tree(trace_fp, trace);

        /* Clean up */
        free(prepseq);
        /* Note: trace is freed by FreeTrace in real code, but we'll skip for simplicity */
    }

    /* Close files */
    fclose(score_fp);
    fclose(trace_fp);

    printf("\nGolden files generated successfully:\n");
    printf("  - viterbi_scores.txt\n");
    printf("  - traceback.txt\n");

    /* Cleanup */
    free(icm);
    FreeCM(cm);

    return 0;
}
