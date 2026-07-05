#include <stdio.h>

#define RANGE 268435456
#define DIV   16384
#define MULT  72530821

static int sre_reseed = 0;
static int sre_randseed = 666;

float sre_random(void) {
    static long rnd;
    static int firsttime = 1;
    long high1, low1, high2, low2;

    if (sre_reseed || firsttime) {
        sre_reseed = firsttime = 0;
        if (sre_randseed <= 0) sre_randseed = 666;
        high1 = sre_randseed / DIV;
        low1 = sre_randseed % DIV;
        high2 = MULT / DIV;
        low2 = MULT % DIV;
        rnd = (((high2*low1 + high1*low2) % DIV)*DIV + low1*low2) % RANGE;
    }
    high1 = rnd / DIV;
    low1 = rnd % DIV;
    high2 = MULT / DIV;
    low2 = MULT % DIV;
    rnd = (((high2*low1 + high1*low2) % DIV)*DIV + low1*low2) % RANGE;
    return ((float) rnd / (float) RANGE);
}

void sre_srandom(int seed) {
    if (seed < 0) seed = -1 * seed;
    sre_reseed = 1;
    sre_randseed = seed;
}

int main() {
    sre_srandom(666);
    for (int i = 0; i < 1000; i++) {
        printf("random[%d] = %.10f\n", i, sre_random());
    }
    return 0;
}
