#include <stdio.h>

#define ABOX_LEN 21
#define BBOX_LEN 11
#define ABDIST_MAT_SIZE 7
#define BTERM_MAT_SIZE 9

int main(void) {
    /* A-Box weight matrix from pavesi.c (lines 41-62) */
    float Abox_Mat[6][ABOX_LEN] = {
        {-1.268,-3.651,-0.899,-4.749,-5.442,-2.351,-3.363,-0.009,-1.977,-3.497,-5.442,
         -5.442,-5.442,-2.498,-4.749,-5.442,-0.031,-1.417,-1.180,-1.048,-4.344},

        {-3.651,-5.442,-4.056,-2.958,-0.480,-1.073,-0.857,-5.442,-5.442,-1.887,-2.498,
         -5.442,-5.442,-2.958,-2.224,-5.442,-5.442,-3.363,-1.417,-3.651,-0.393},

        {-0.779,-5.442,-0.598,-0.076,-3.651,-1.435,-1.614,-4.749,-0.154,-2.803,-5.442,
         0.000,0.000,-3.363,-3.651,-5.442,-3.497,-0.672,-1.012,-0.473,-3.651},

        {-1.453,-0.026,-3.651,-4.344,-1.036,-1.125,-1.073,-5.442,-5.442,-0.278,-1.399,
         -5.442,-5.442,-0.185,-0.827,-2.041,-5.442,-1.551,-2.447,-5.442,-1.253},

        {-5.442,-5.442,-5.442,-5.442,-5.442,-5.442,-5.442,-5.442,-5.442,-5.442,-0.412,
         -5.442,-5.442,-5.442,-0.868,-0.144,-5.442,-5.442,-5.442,-5.442,-5.442},

        {-0.779,-0.026,-0.598,-0.076,-0.480,-1.073,-0.857,-0.009,-0.154,-0.278,-1.399,
         0.000,0.000,-0.185,-0.827,-2.041,-0.031,-0.672,-1.012,-0.473,-0.393}
    };

    /* B-Box weight matrix from pavesi.c (lines 67-77) */
    float Bbox_Mat[6][BBOX_LEN] = {
        {-2.351,-5.442,-2.670,-5.442,-5.442,-1.472,0.000,-0.798,-2.498,-5.442,-3.497},
        {-3.245,-5.442,-5.442,-5.442,-0.004,-5.442,-5.442,-2.498,-1.435,-0.009,-0.190},
        {-0.175,-0.004,-5.442,-5.442,-5.442,-0.272,-5.442,-2.147,-5.442,-5.442,-3.651},
        {-3.651,-5.442,-0.072,0.000,-5.442,-4.749,-5.442,-1.048,-0.393,-5.442,-2.147},
        {-5.442,-5.442,-5.442,-5.442,-5.442,-5.442,-5.442,-5.442,-5.442,-5.442,-5.442},
        {-0.175,-0.004,-0.072,0.000,-0.004,-0.272,0.000,-0.798,-0.393,-0.009,-0.190}
    };

    /* ABDistIdx_Mat from pavesi.c (line 81) */
    int ABDistIdx_Mat[ABDIST_MAT_SIZE] = {30,36,42,48,54,60,66};

    /* ABDistSc_Mat from pavesi.c (lines 82-83) */
    float ABDistSc_Mat[ABDIST_MAT_SIZE] = {-0.46,-1.83,-2.35,-3.24,
                                           -4.06,-3.83,-4.75};

    /* BTermDistIdx_Mat from pavesi.c (line 86) */
    int BTermDistIdx_Mat[BTERM_MAT_SIZE] = {17,23,29,35,41,47,53,59,100};

    /* BTermDistSc_Mat from pavesi.c (lines 88-89) */
    float BTermDistSc_Mat[BTERM_MAT_SIZE] = {-0.54,-1.40,-2.80,-3.36,
                                             -3.24,-5.44,-5.44,-4.06,-5.44};

    int i, j;

    /* Print A-Box matrix */
    for (i = 0; i < 6; i++) {
        for (j = 0; j < ABOX_LEN; j++) {
            printf("Abox_Mat[%d][%d] = %g\n", i, j, Abox_Mat[i][j]);
        }
    }

    /* Print B-Box matrix */
    for (i = 0; i < 6; i++) {
        for (j = 0; j < BBOX_LEN; j++) {
            printf("Bbox_Mat[%d][%d] = %g\n", i, j, Bbox_Mat[i][j]);
        }
    }

    /* Print ABDistIdx_Mat */
    for (i = 0; i < ABDIST_MAT_SIZE; i++) {
        printf("ABDistIdx_Mat[%d] = %d\n", i, ABDistIdx_Mat[i]);
    }

    /* Print ABDistSc_Mat */
    for (i = 0; i < ABDIST_MAT_SIZE; i++) {
        printf("ABDistSc_Mat[%d] = %g\n", i, ABDistSc_Mat[i]);
    }

    /* Print BTermDistIdx_Mat */
    for (i = 0; i < BTERM_MAT_SIZE; i++) {
        printf("BTermDistIdx_Mat[%d] = %d\n", i, BTermDistIdx_Mat[i]);
    }

    /* Print BTermDistSc_Mat */
    for (i = 0; i < BTERM_MAT_SIZE; i++) {
        printf("BTermDistSc_Mat[%d] = %g\n", i, BTermDistSc_Mat[i]);
    }

    return 0;
}
