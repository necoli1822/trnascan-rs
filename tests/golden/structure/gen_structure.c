/* gen_structure.c
 *
 * Generate golden files for Phase 7 (tRNA secondary structure) testing
 * Tests the following konings.c functions:
 * - Trace2KHS() - Convert traceback to secondary structure
 * - KHS2ct() - Convert KHS format to CT format
 * - IsRNAComplement() - Check if bases can pair
 */

#include <stdio.h>
#include <stdlib.h>
#include <string.h>

/* Simple trace node structure for testing */
struct trace_s {
  int emitl;
  int emitr;
  int nodeidx;
  int type;
  struct trace_s *nxtl;
  struct trace_s *nxtr;
  struct trace_s *prv;
};

/* State types from structs.h */
#define uMATP_ST    (1<<1)
#define uMATL_ST    (1<<2)
#define uMATR_ST    (1<<3)
#define uDEL_ST     (1<<0)
#define uINSL_ST    (1<<4)
#define uINSR_ST    (1<<5)

#define TRUE  1
#define FALSE 0

/* Minimal trace stack implementation */
struct tracestack_s {
  int next;
  int num;
  struct trace_s **list;
};

struct tracestack_s *InitTracestack(void) {
  struct tracestack_s *stack = (struct tracestack_s *) malloc(sizeof(struct tracestack_s));
  stack->next = 0;
  stack->num = 64;
  stack->list = (struct trace_s **) malloc(sizeof(struct trace_s *) * stack->num);
  return stack;
}

void PushTracestack(struct tracestack_s *stack, struct trace_s *node) {
  if (stack->next >= stack->num) {
    stack->num *= 2;
    stack->list = (struct trace_s **) realloc(stack->list, sizeof(struct trace_s *) * stack->num);
  }
  stack->list[stack->next++] = node;
}

struct trace_s *PopTracestack(struct tracestack_s *stack) {
  if (stack->next == 0) return NULL;
  return stack->list[--stack->next];
}

void FreeTracestack(struct tracestack_s *stack) {
  free(stack->list);
  free(stack);
}

/* Minimal int stack implementation */
struct intstack_s {
  int next;
  int num;
  int *list;
};

struct intstack_s *InitIntStack(void) {
  struct intstack_s *stack = (struct intstack_s *) malloc(sizeof(struct intstack_s));
  stack->next = 0;
  stack->num = 64;
  stack->list = (int *) malloc(sizeof(int) * stack->num);
  return stack;
}

void PushIntStack(struct intstack_s *stack, int data) {
  if (stack->next >= stack->num) {
    stack->num *= 2;
    stack->list = (int *) realloc(stack->list, sizeof(int) * stack->num);
  }
  stack->list[stack->next++] = data;
}

int PopIntStack(struct intstack_s *stack, int *ret_data) {
  if (stack->next == 0) return FALSE;
  *ret_data = stack->list[--stack->next];
  return TRUE;
}

int FreeIntStack(struct intstack_s *stack) {
  int remaining = stack->next;
  free(stack->list);
  free(stack);
  return remaining;
}

/* Function: IsRNAComplement()
 * Returns TRUE if sym1, sym2 are Watson-Crick complementary.
 * If allow_gu is TRUE, GU pairs also return TRUE.
 */
int IsRNAComplement(char sym1, char sym2, int allow_gu) {
  if (sym1 >= 'a' && sym1 <= 'z') sym1 = sym1 - 'a' + 'A';
  if (sym2 >= 'a' && sym2 <= 'z') sym2 = sym2 - 'a' + 'A';
  if (sym1 == 'T') sym1 = 'U';
  if (sym2 == 'T') sym2 = 'U';

  if ((sym1 == 'A' && sym2 == 'U') ||
      (sym1 == 'C' && sym2 == 'G') ||
      (sym1 == 'G' && sym2 == 'C') ||
      (sym1 == 'U' && sym2 == 'A') ||
      (allow_gu && sym1 == 'G' && sym2 == 'U') ||
      (allow_gu && sym1 == 'U' && sym2 == 'G'))
    return TRUE;
  else
    return FALSE;
}

/* Function: Trace2KHS()
 * Convert a traceback tree to a secondary structure string.
 */
void Trace2KHS(struct trace_s *tr, char *seq, int rlen, int watsoncrick, char **ret_ss) {
  struct tracestack_s *dolist;
  struct trace_s *curr;
  char *ss;

  ss = (char *) malloc(sizeof(char) * (rlen + 1));
  memset(ss, '.', rlen);
  ss[rlen] = '\0';

  dolist = InitTracestack();
  if (tr->nxtl) PushTracestack(dolist, tr->nxtl);

  while ((curr = PopTracestack(dolist)) != NULL) {
    if (curr->type == uMATP_ST) {
      if (!watsoncrick || IsRNAComplement(seq[curr->emitl], seq[curr->emitr], TRUE)) {
        ss[curr->emitl] = '>';
        ss[curr->emitr] = '<';
      }
    }

    if (curr->nxtr) PushTracestack(dolist, curr->nxtr);
    if (curr->nxtl) PushTracestack(dolist, curr->nxtl);
  }

  FreeTracestack(dolist);
  *ret_ss = ss;
}

/* Function: KHS2ct()
 * Convert a secondary structure string to a CT array.
 */
int KHS2ct(char *ss, int len, int allow_pseudoknots, int **ret_ct) {
  struct intstack_s *dolist[27];
  int *ct;
  int i;
  int pos, pair;
  int status = TRUE;

  for (i = 0; i < 27; i++)
    dolist[i] = InitIntStack();

  ct = (int *) malloc(len * sizeof(int));
  for (pos = 0; pos < len; pos++)
    ct[pos] = -1;

  for (pos = 0; ss[pos] != '\0'; pos++) {
    if (ss[pos] == '>') {
      PushIntStack(dolist[0], pos);
    }
    else if (ss[pos] == '<') {
      if (!PopIntStack(dolist[0], &pair)) {
        status = FALSE;
      } else {
        ct[pos] = pair;
        ct[pair] = pos;
      }
    }
    else if (allow_pseudoknots && ss[pos] >= 'A' && ss[pos] <= 'Z') {
      PushIntStack(dolist[ss[pos] - 'A' + 1], pos);
    }
    else if (allow_pseudoknots && ss[pos] >= 'a' && ss[pos] <= 'z') {
      if (!PopIntStack(dolist[ss[pos] - 'a' + 1], &pair)) {
        status = FALSE;
      } else {
        ct[pos] = pair;
        ct[pair] = pos;
      }
    }
  }

  for (i = 0; i < 27; i++)
    if (FreeIntStack(dolist[i]) > 0)
      status = FALSE;

  *ret_ct = ct;
  return status;
}

/* Helper to create a trace node */
struct trace_s *create_trace_node(int emitl, int emitr, int nodeidx, int type) {
  struct trace_s *node = (struct trace_s *) malloc(sizeof(struct trace_s));
  node->emitl = emitl;
  node->emitr = emitr;
  node->nodeidx = nodeidx;
  node->type = type;
  node->nxtl = NULL;
  node->nxtr = NULL;
  node->prv = NULL;
  return node;
}

/* Test 1: Trace2KHS with a simple tRNA cloverleaf structure */
void test_trace2khs(FILE *fp) {
  fprintf(fp, "=== Test 1: Trace2KHS() - Convert traceback to KHS structure ===\n\n");

  /* Example tRNA sequence (76 nt) */
  char *seq = "GCGGAUUUAGCUCAGUUGGGAGAGCGCCAGACUGAAGAUCUGGAGGUCCUGUGUUCGAUCCACAGAAUUCGCACCA";
  int rlen = strlen(seq);

  fprintf(fp, "Test sequence (%d nt):\n%s\n\n", rlen, seq);

  /* Build a simplified traceback tree for tRNA structure
   * Positions are 0-indexed
   * Typical tRNA has these base pairs:
   * Acceptor stem: 0-6 pairs with 69-75
   * D stem: 10-12 pairs with 22-24
   * Anticodon stem: 27-32 pairs with 38-43
   * T stem: 49-52 pairs with 61-64
   */

  struct trace_s *root = create_trace_node(-1, -1, 0, uDEL_ST);

  /* Create pairing nodes (simplified) */
  struct trace_s *pair1 = create_trace_node(0, 75, 1, uMATP_ST);
  struct trace_s *pair2 = create_trace_node(1, 74, 2, uMATP_ST);
  struct trace_s *pair3 = create_trace_node(2, 73, 3, uMATP_ST);
  struct trace_s *pair4 = create_trace_node(10, 24, 4, uMATP_ST);
  struct trace_s *pair5 = create_trace_node(11, 23, 5, uMATP_ST);
  struct trace_s *pair6 = create_trace_node(27, 43, 6, uMATP_ST);
  struct trace_s *pair7 = create_trace_node(28, 42, 7, uMATP_ST);
  struct trace_s *pair8 = create_trace_node(49, 64, 8, uMATP_ST);
  struct trace_s *pair9 = create_trace_node(50, 63, 9, uMATP_ST);

  /* Link nodes in tree */
  root->nxtl = pair1;
  pair1->nxtl = pair2;
  pair2->nxtl = pair3;
  pair3->nxtl = pair4;
  pair4->nxtl = pair5;
  pair5->nxtl = pair6;
  pair6->nxtl = pair7;
  pair7->nxtl = pair8;
  pair8->nxtl = pair9;

  /* Test with Watson-Crick enforcement */
  char *ss_wc;
  Trace2KHS(root, seq, rlen, TRUE, &ss_wc);
  fprintf(fp, "Watson-Crick only (watsoncrick=TRUE):\n%s\n\n", ss_wc);

  /* Test without Watson-Crick enforcement */
  char *ss_all;
  Trace2KHS(root, seq, rlen, FALSE, &ss_all);
  fprintf(fp, "All pairs (watsoncrick=FALSE):\n%s\n\n", ss_all);

  /* Clean up */
  free(ss_wc);
  free(ss_all);
  free(pair1);
  free(pair2);
  free(pair3);
  free(pair4);
  free(pair5);
  free(pair6);
  free(pair7);
  free(pair8);
  free(pair9);
  free(root);
}

/* Test 2: KHS2ct conversion */
void test_khs2ct(FILE *fp) {
  fprintf(fp, "=== Test 2: KHS2ct() - Convert KHS to CT format ===\n\n");

  /* Test case 1: Simple hairpin */
  char *ss1 = ">>>....<<<";
  int len1 = strlen(ss1);
  int *ct1;
  int status1 = KHS2ct(ss1, len1, FALSE, &ct1);

  fprintf(fp, "Test case 1: Simple hairpin\n");
  fprintf(fp, "Structure: %s\n", ss1);
  fprintf(fp, "Status: %s\n", status1 ? "SUCCESS" : "FAILED");
  fprintf(fp, "CT array:\n");
  for (int i = 0; i < len1; i++) {
    fprintf(fp, "  Position %d: paired to %d\n", i, ct1[i]);
  }
  fprintf(fp, "\n");
  free(ct1);

  /* Test case 2: tRNA-like structure */
  char *ss2 = ">>>>>>>......<<<.....>>>....<<<........>>>>>>>......<<<<<<<.......<<<<<<";
  int len2 = strlen(ss2);
  int *ct2;
  int status2 = KHS2ct(ss2, len2, FALSE, &ct2);

  fprintf(fp, "Test case 2: tRNA-like structure\n");
  fprintf(fp, "Structure: %s\n", ss2);
  fprintf(fp, "Status: %s\n", status2 ? "SUCCESS" : "FAILED");
  fprintf(fp, "CT array (first 20 positions):\n");
  for (int i = 0; i < 20 && i < len2; i++) {
    fprintf(fp, "  Position %d: paired to %d\n", i, ct2[i]);
  }
  fprintf(fp, "\n");
  free(ct2);

  /* Test case 3: With pseudoknots */
  char *ss3 = ">>>AAA<<<aaa";
  int len3 = strlen(ss3);
  int *ct3;
  int status3 = KHS2ct(ss3, len3, TRUE, &ct3);

  fprintf(fp, "Test case 3: Pseudoknot structure\n");
  fprintf(fp, "Structure: %s\n", ss3);
  fprintf(fp, "Status: %s\n", status3 ? "SUCCESS" : "FAILED");
  fprintf(fp, "CT array:\n");
  for (int i = 0; i < len3; i++) {
    fprintf(fp, "  Position %d: paired to %d\n", i, ct3[i]);
  }
  fprintf(fp, "\n");
  free(ct3);
}

/* Test 3: IsRNAComplement */
void test_rna_complement(FILE *fp) {
  fprintf(fp, "=== Test 3: IsRNAComplement() - Base pairing rules ===\n\n");

  /* Watson-Crick pairs */
  fprintf(fp, "Watson-Crick pairs (allow_gu=FALSE):\n");
  fprintf(fp, "  A-U: %s\n", IsRNAComplement('A', 'U', FALSE) ? "TRUE" : "FALSE");
  fprintf(fp, "  U-A: %s\n", IsRNAComplement('U', 'A', FALSE) ? "TRUE" : "FALSE");
  fprintf(fp, "  G-C: %s\n", IsRNAComplement('G', 'C', FALSE) ? "TRUE" : "FALSE");
  fprintf(fp, "  C-G: %s\n", IsRNAComplement('C', 'G', FALSE) ? "TRUE" : "FALSE");
  fprintf(fp, "  G-U: %s (should be FALSE)\n", IsRNAComplement('G', 'U', FALSE) ? "TRUE" : "FALSE");
  fprintf(fp, "\n");

  /* With GU wobble pairs */
  fprintf(fp, "With GU wobble pairs (allow_gu=TRUE):\n");
  fprintf(fp, "  A-U: %s\n", IsRNAComplement('A', 'U', TRUE) ? "TRUE" : "FALSE");
  fprintf(fp, "  G-C: %s\n", IsRNAComplement('G', 'C', TRUE) ? "TRUE" : "FALSE");
  fprintf(fp, "  G-U: %s (should be TRUE)\n", IsRNAComplement('G', 'U', TRUE) ? "TRUE" : "FALSE");
  fprintf(fp, "  U-G: %s (should be TRUE)\n", IsRNAComplement('U', 'G', TRUE) ? "TRUE" : "FALSE");
  fprintf(fp, "\n");

  /* Non-complementary pairs */
  fprintf(fp, "Non-complementary pairs:\n");
  fprintf(fp, "  A-A: %s\n", IsRNAComplement('A', 'A', TRUE) ? "TRUE" : "FALSE");
  fprintf(fp, "  A-G: %s\n", IsRNAComplement('A', 'G', TRUE) ? "TRUE" : "FALSE");
  fprintf(fp, "  C-C: %s\n", IsRNAComplement('C', 'C', TRUE) ? "TRUE" : "FALSE");
  fprintf(fp, "\n");

  /* Case insensitive testing */
  fprintf(fp, "Case insensitive:\n");
  fprintf(fp, "  a-u: %s\n", IsRNAComplement('a', 'u', TRUE) ? "TRUE" : "FALSE");
  fprintf(fp, "  G-c: %s\n", IsRNAComplement('G', 'c', TRUE) ? "TRUE" : "FALSE");
  fprintf(fp, "\n");

  /* DNA (T instead of U) */
  fprintf(fp, "DNA compatibility (T converted to U):\n");
  fprintf(fp, "  A-T: %s\n", IsRNAComplement('A', 'T', TRUE) ? "TRUE" : "FALSE");
  fprintf(fp, "  T-A: %s\n", IsRNAComplement('T', 'A', TRUE) ? "TRUE" : "FALSE");
  fprintf(fp, "  G-T: %s\n", IsRNAComplement('G', 'T', TRUE) ? "TRUE" : "FALSE");
  fprintf(fp, "\n");
}

int main(int argc, char **argv) {
  FILE *fp_khs, *fp_ct, *fp_comp;

  /* Generate khs_output.txt */
  fp_khs = fopen("tests/golden/structure/khs_output.txt", "w");
  if (!fp_khs) {
    fprintf(stderr, "Error: Cannot open khs_output.txt\n");
    return 1;
  }
  fprintf(fp_khs, "Phase 7: tRNA Secondary Structure - Trace2KHS Test Output\n");
  fprintf(fp_khs, "Generated by gen_structure.c\n\n");
  test_trace2khs(fp_khs);
  fclose(fp_khs);
  printf("Generated: tests/golden/structure/khs_output.txt\n");

  /* Generate ct_output.txt */
  fp_ct = fopen("tests/golden/structure/ct_output.txt", "w");
  if (!fp_ct) {
    fprintf(stderr, "Error: Cannot open ct_output.txt\n");
    return 1;
  }
  fprintf(fp_ct, "Phase 7: tRNA Secondary Structure - KHS2ct Test Output\n");
  fprintf(fp_ct, "Generated by gen_structure.c\n\n");
  test_khs2ct(fp_ct);
  fclose(fp_ct);
  printf("Generated: tests/golden/structure/ct_output.txt\n");

  /* Generate rna_complement.txt */
  fp_comp = fopen("tests/golden/structure/rna_complement.txt", "w");
  if (!fp_comp) {
    fprintf(stderr, "Error: Cannot open rna_complement.txt\n");
    return 1;
  }
  fprintf(fp_comp, "Phase 7: tRNA Secondary Structure - IsRNAComplement Test Output\n");
  fprintf(fp_comp, "Generated by gen_structure.c\n\n");
  test_rna_complement(fp_comp);
  fclose(fp_comp);
  printf("Generated: tests/golden/structure/rna_complement.txt\n");

  printf("\nAll golden files generated successfully!\n");
  return 0;
}
