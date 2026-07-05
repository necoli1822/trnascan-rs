/* gen_cm.c - Generate golden files for CM structure parsing
 *
 * This program reads a CM file and dumps its internal structure
 * to verify Rust reimplementation correctness.
 *
 * Generates:
 *   - cm_structure.txt: Full CM structure dump
 *   - node_details.txt: Per-node type, transitions, emissions
 *   - istate_dump.txt: Integer state array after RearrangeCM
 */

#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include "funcs.h"
#include "structs.h"

/* Function prototypes are already in funcs.h */

/* Node type names for clarity */
static const char *node_type_names[] = {
    "BIFURC_NODE",
    "MATP_NODE",
    "MATL_NODE",
    "MATR_NODE",
    "BEGINL_NODE",
    "BEGINR_NODE",
    "ROOT_NODE"
};

/* State type names */
static const char *state_type_names[] = {
    "DEL_ST",
    "MATP_ST",
    "MATL_ST",
    "MATR_ST",
    "INSL_ST",
    "INSR_ST"
};

/* Unique state type names for flags */
static const char *unique_state_names(int statetype) {
    switch (statetype) {
        case uDEL_ST:    return "uDEL_ST";
        case uMATP_ST:   return "uMATP_ST";
        case uMATL_ST:   return "uMATL_ST";
        case uMATR_ST:   return "uMATR_ST";
        case uINSL_ST:   return "uINSL_ST";
        case uINSR_ST:   return "uINSR_ST";
        case uBEGIN_ST:  return "uBEGIN_ST";
        case uEND_ST:    return "uEND_ST";
        case uBIFURC_ST: return "uBIFURC_ST";
        default:         return "UNKNOWN";
    }
}

void dump_cm_structure(struct cm_s *cm, const char *filename) {
    FILE *fp = fopen(filename, "w");
    if (!fp) {
        perror("fopen cm_structure.txt");
        return;
    }

    fprintf(fp, "CM Structure Dump\n");
    fprintf(fp, "=================\n\n");
    fprintf(fp, "Number of nodes: %d\n\n", cm->nodes);

    fprintf(fp, "Constants:\n");
    fprintf(fp, "  STATETYPES = %d\n", STATETYPES);
    fprintf(fp, "  ALPHASIZE = %d\n", ALPHASIZE);
    fprintf(fp, "  NODETYPES = %d\n\n", NODETYPES);

    for (int i = 0; i < cm->nodes; i++) {
        struct node_s *nd = &cm->nd[i];
        fprintf(fp, "Node %d: type=%d (%s), nxt=%d, nxt2=%d\n",
                i, nd->type, node_type_names[nd->type], nd->nxt, nd->nxt2);
    }

    fclose(fp);
    printf("Generated: %s\n", filename);
}

void dump_node_details(struct cm_s *cm, const char *filename) {
    FILE *fp = fopen(filename, "w");
    if (!fp) {
        perror("fopen node_details.txt");
        return;
    }

    fprintf(fp, "Node Details Dump\n");
    fprintf(fp, "=================\n\n");

    for (int i = 0; i < cm->nodes; i++) {
        struct node_s *nd = &cm->nd[i];

        fprintf(fp, "Node %d: %s\n", i, node_type_names[nd->type]);
        fprintf(fp, "  Connections: nxt=%d, nxt2=%d\n", nd->nxt, nd->nxt2);

        /* Transition matrix */
        fprintf(fp, "  Transition matrix [%dx%d]:\n", STATETYPES, STATETYPES);
        for (int row = 0; row < STATETYPES; row++) {
            fprintf(fp, "    %s -> ", state_type_names[row]);
            for (int col = 0; col < STATETYPES; col++) {
                fprintf(fp, "%10.6f ", nd->tmx[row][col]);
            }
            fprintf(fp, "\n");
        }

        /* Emission probabilities */
        fprintf(fp, "  MATP emissions [%dx%d]:\n", ALPHASIZE, ALPHASIZE);
        for (int row = 0; row < ALPHASIZE; row++) {
            fprintf(fp, "    ");
            for (int col = 0; col < ALPHASIZE; col++) {
                fprintf(fp, "%10.6f ", nd->mp_emit[row][col]);
            }
            fprintf(fp, "\n");
        }

        fprintf(fp, "  INSL emissions: ");
        for (int j = 0; j < ALPHASIZE; j++) {
            fprintf(fp, "%10.6f ", nd->il_emit[j]);
        }
        fprintf(fp, "\n");

        fprintf(fp, "  INSR emissions: ");
        for (int j = 0; j < ALPHASIZE; j++) {
            fprintf(fp, "%10.6f ", nd->ir_emit[j]);
        }
        fprintf(fp, "\n");

        fprintf(fp, "  MATL emissions: ");
        for (int j = 0; j < ALPHASIZE; j++) {
            fprintf(fp, "%10.6f ", nd->ml_emit[j]);
        }
        fprintf(fp, "\n");

        fprintf(fp, "  MATR emissions: ");
        for (int j = 0; j < ALPHASIZE; j++) {
            fprintf(fp, "%10.6f ", nd->mr_emit[j]);
        }
        fprintf(fp, "\n\n");
    }

    fclose(fp);
    printf("Generated: %s\n", filename);
}

void dump_istate_array(struct istate_s *istate, int numstates, const char *filename) {
    FILE *fp = fopen(filename, "w");
    if (!fp) {
        perror("fopen istate_dump.txt");
        return;
    }

    fprintf(fp, "Integer State Array Dump (after RearrangeCM)\n");
    fprintf(fp, "============================================\n\n");
    fprintf(fp, "Total states: %d\n\n", numstates);

    fprintf(fp, "Field verification:\n");
    fprintf(fp, "  istate_s has nodeidx: YES\n");
    fprintf(fp, "  istate_s has bifr: NO (only in pstate_s)\n\n");

    for (int i = 0; i < numstates; i++) {
        struct istate_s *st = &istate[i];

        fprintf(fp, "State %d:\n", i);
        fprintf(fp, "  nodeidx: %d\n", st->nodeidx);
        fprintf(fp, "  statetype: 0x%x (%s)\n", st->statetype, unique_state_names(st->statetype));
        fprintf(fp, "  offset: %d\n", st->offset);
        fprintf(fp, "  connectnum: %d\n", st->connectnum);

        fprintf(fp, "  tmx[%d]: ", STATETYPES);
        for (int j = 0; j < STATETYPES; j++) {
            fprintf(fp, "%d ", st->tmx[j]);
        }
        fprintf(fp, "\n");

        /* Determine emit array size based on state type */
        int emit_size = 0;
        if (st->statetype == uMATP_ST) {
            emit_size = ALPHASIZE * ALPHASIZE;
        } else if (st->statetype == uMATL_ST || st->statetype == uMATR_ST ||
                   st->statetype == uINSL_ST || st->statetype == uINSR_ST) {
            emit_size = ALPHASIZE;
        }

        if (emit_size > 0) {
            fprintf(fp, "  emit[%d]: ", emit_size);
            for (int j = 0; j < emit_size; j++) {
                fprintf(fp, "%d ", st->emit[j]);
            }
            fprintf(fp, "\n");
        } else {
            fprintf(fp, "  emit: (none)\n");
        }

        fprintf(fp, "\n");
    }

    fclose(fp);
    printf("Generated: %s\n", filename);
}

void verify_magic_constants() {
    printf("\nMagic constant verification:\n");
    printf("  v20magic = 0xe3edb2b0 (expected binary format magic)\n");
    printf("  STATETYPES = %d (expected: 6)\n", STATETYPES);
    printf("  ALPHASIZE = %d (expected: 4)\n", ALPHASIZE);
    printf("  NODETYPES = %d (expected: 7)\n", NODETYPES);

    if (STATETYPES != 6 || ALPHASIZE != 4 || NODETYPES != 7) {
        fprintf(stderr, "ERROR: Constants do not match expected values!\n");
        exit(1);
    }
}

int main(int argc, char **argv) {
    if (argc != 2) {
        fprintf(stderr, "Usage: %s <cmfile>\n", argv[0]);
        fprintf(stderr, "  Generates golden files in current directory:\n");
        fprintf(stderr, "    cm_structure.txt\n");
        fprintf(stderr, "    node_details.txt\n");
        fprintf(stderr, "    istate_dump.txt\n");
        return 1;
    }

    const char *cmfile = argv[1];

    printf("Reading CM file: %s\n", cmfile);
    verify_magic_constants();

    /* Read the CM */
    struct cm_s *cm = NULL;
    if (!ReadCM((char *)cmfile, &cm)) {
        fprintf(stderr, "ERROR: Failed to read CM file\n");
        return 1;
    }
    if (!cm) {
        fprintf(stderr, "ERROR: CM is NULL after ReadCM\n");
        return 1;
    }

    printf("Successfully read CM with %d nodes\n", cm->nodes);

    /* Generate golden files */
    dump_cm_structure(cm, "cm_structure.txt");
    dump_node_details(cm, "node_details.txt");

    /* Rearrange CM to integer state array */
    int numstates = 0;
    struct istate_s *istate = NULL;
    double rfreq[ALPHASIZE] = {0.25, 0.25, 0.25, 0.25}; /* uniform background */

    if (!RearrangeCM(cm, rfreq, &istate, &numstates)) {
        fprintf(stderr, "ERROR: RearrangeCM failed\n");
        FreeCM(cm);
        return 1;
    }
    if (!istate) {
        fprintf(stderr, "ERROR: istate is NULL after RearrangeCM\n");
        FreeCM(cm);
        return 1;
    }

    printf("RearrangeCM produced %d states\n", numstates);
    dump_istate_array(istate, numstates, "istate_dump.txt");

    /* Cleanup */
    free(istate);
    FreeCM(cm);

    printf("\nAll golden files generated successfully!\n");
    return 0;
}
