#include <stdio.h>
#include <math.h>

#define INTPRECISION 1000.0
#define NEGINFINITY -999999
#define ILOG2(a) (((a) > 0.0) ? (int)(log(a) / 0.69314718 * INTPRECISION) : NEGINFINITY)

int main() {
    double test_values[] = {
        0.0,
        -1.0,
        0.5,
        0.25,
        0.125,
        1.0,
        2.0,
        4.0,
        8.0,
        0.001,
        0.999,
        1e-10,
        0.333333
    };

    int num_tests = sizeof(test_values) / sizeof(test_values[0]);

    for (int i = 0; i < num_tests; i++) {
        double value = test_values[i];
        int result = ILOG2(value);
        printf("ILOG2(%.15g) = %d\n", value, result);
    }

    return 0;
}
