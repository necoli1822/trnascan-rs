/* gen_squid.c
 * Golden file generator for SQUID library functions
 * Tests sequence format detection and SQINFO structure parsing
 */

#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include "squid.h"

void print_format_name(int format, FILE *fp) {
    switch(format) {
        case kUnknown:   fprintf(fp, "kUnknown"); break;
        case kIG:        fprintf(fp, "kIG"); break;
        case kGenBank:   fprintf(fp, "kGenBank"); break;
        case kNBRF:      fprintf(fp, "kNBRF"); break;
        case kEMBL:      fprintf(fp, "kEMBL"); break;
        case kGCG:       fprintf(fp, "kGCG"); break;
        case kStrider:   fprintf(fp, "kStrider"); break;
        case kPearson:   fprintf(fp, "kPearson"); break;
        case kZuker:     fprintf(fp, "kZuker"); break;
        case kIdraw:     fprintf(fp, "kIdraw"); break;
        case kSelex:     fprintf(fp, "kSelex"); break;
        case kMSF:       fprintf(fp, "kMSF"); break;
        case kPIR:       fprintf(fp, "kPIR"); break;
        case kRaw:       fprintf(fp, "kRaw"); break;
        case kSquid:     fprintf(fp, "kSquid"); break;
        case kXPearson:  fprintf(fp, "kXPearson"); break;
        case kGCGdata:   fprintf(fp, "kGCGdata"); break;
        case kClustal:   fprintf(fp, "kClustal"); break;
        default:         fprintf(fp, "UNKNOWN_FORMAT"); break;
    }
}

void test_format_detection(const char *filepath, FILE *out) {
    int format;
    int ret = SeqfileFormat((char*)filepath, &format, NULL);

    fprintf(out, "File: %s\n", filepath);
    fprintf(out, "  Return: %d\n", ret);
    fprintf(out, "  Format code: %d (", format);
    print_format_name(format, out);
    fprintf(out, ")\n");
    fprintf(out, "\n");
}

void test_read_seq(const char *filepath, int expected_format, FILE *out) {
    SQFILE *sqfp;
    char *seq;
    SQINFO sqinfo;
    int count = 0;

    fprintf(out, "=== Reading file: %s ===\n", filepath);
    fprintf(out, "Expected format: %d (", expected_format);
    print_format_name(expected_format, out);
    fprintf(out, ")\n\n");

    sqfp = SeqfileOpen((char*)filepath, expected_format, NULL);
    if (sqfp == NULL) {
        fprintf(out, "ERROR: Could not open file\n\n");
        return;
    }

    while (ReadSeq(sqfp, expected_format, &seq, &sqinfo)) {
        count++;
        fprintf(out, "Sequence #%d:\n", count);

        fprintf(out, "  flags: 0x%x\n", sqinfo.flags);

        if (sqinfo.flags & SQINFO_NAME) {
            fprintf(out, "  name: \"%s\"\n", sqinfo.name);
        }

        if (sqinfo.flags & SQINFO_ID) {
            fprintf(out, "  id: \"%s\"\n", sqinfo.id);
        }

        if (sqinfo.flags & SQINFO_ACC) {
            fprintf(out, "  acc: \"%s\"\n", sqinfo.acc);
        }

        if (sqinfo.flags & SQINFO_DESC) {
            fprintf(out, "  desc: \"%s\"\n", sqinfo.desc);
        }

        if (sqinfo.flags & SQINFO_LEN) {
            fprintf(out, "  len: %d\n", sqinfo.len);
        }

        if (sqinfo.flags & SQINFO_START) {
            fprintf(out, "  start: %d\n", sqinfo.start);
        }

        if (sqinfo.flags & SQINFO_STOP) {
            fprintf(out, "  stop: %d\n", sqinfo.stop);
        }

        if (sqinfo.flags & SQINFO_OLEN) {
            fprintf(out, "  olen: %d\n", sqinfo.olen);
        }

        if (sqinfo.flags & SQINFO_TYPE) {
            fprintf(out, "  type: %d", sqinfo.type);
            switch(sqinfo.type) {
                case kDNA:   fprintf(out, " (kDNA)"); break;
                case kRNA:   fprintf(out, " (kRNA)"); break;
                case kAmino: fprintf(out, " (kAmino)"); break;
                case kOtherSeq: fprintf(out, " (kOtherSeq)"); break;
                default:     fprintf(out, " (UNKNOWN)"); break;
            }
            fprintf(out, "\n");
        }

        if (sqinfo.flags & SQINFO_WGT) {
            fprintf(out, "  weight: %.3f\n", sqinfo.weight);
        }

        fprintf(out, "  sequence: %s\n", seq);
        fprintf(out, "\n");

        FreeSequence(seq, &sqinfo);
    }

    fprintf(out, "Total sequences read: %d\n", count);
    fprintf(out, "\n");

    SeqfileClose(sqfp);
}

int main(int argc, char **argv) {
    FILE *format_out;
    FILE *sqinfo_out;

    printf("SQUID Library Golden File Generator\n");
    printf("====================================\n\n");

    /* Test format detection */
    format_out = fopen("tests/golden/squid/format_detection.txt", "w");
    if (!format_out) {
        fprintf(stderr, "Error: Cannot create format_detection.txt\n");
        return 1;
    }

    fprintf(format_out, "SQUID Format Detection Test Results\n");
    fprintf(format_out, "====================================\n\n");

    test_format_detection("tests/golden/squid/inputs/test.fasta", format_out);
    test_format_detection("tests/golden/squid/inputs/test.gb", format_out);
    test_format_detection("tests/golden/squid/inputs/test.embl", format_out);
    test_format_detection("tests/golden/squid/inputs/test.raw", format_out);

    fclose(format_out);
    printf("Created: tests/golden/squid/format_detection.txt\n");

    /* Test SQINFO parsing */
    sqinfo_out = fopen("tests/golden/squid/sqinfo_fields.txt", "w");
    if (!sqinfo_out) {
        fprintf(stderr, "Error: Cannot create sqinfo_fields.txt\n");
        return 1;
    }

    fprintf(sqinfo_out, "SQUID ReadSeq and SQINFO Test Results\n");
    fprintf(sqinfo_out, "======================================\n\n");

    test_read_seq("tests/golden/squid/inputs/test.fasta", kPearson, sqinfo_out);
    test_read_seq("tests/golden/squid/inputs/test.gb", kGenBank, sqinfo_out);
    test_read_seq("tests/golden/squid/inputs/test.embl", kEMBL, sqinfo_out);

    fclose(sqinfo_out);
    printf("Created: tests/golden/squid/sqinfo_fields.txt\n");

    printf("\nGolden file generation complete!\n");
    return 0;
}
