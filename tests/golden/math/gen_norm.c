#include <stdio.h>
#include <stdlib.h>
#include <math.h>

/* DNorm - normalize a double vector to sum to 1.0 */
int DNorm(double *vec, int n)
{
  int x;
  double sum;

  sum = 0.0;
  for (x = 0; x < n; x++) sum += vec[x];
  if (sum != 0.0)
    for (x = 0; x < n; x++) vec[x] /= sum;
  else
    { return 0; }  /* return 0 if all values were zero */
  return 1;
}

/* FNorm - normalize a float vector to sum to 1.0 */
int FNorm(float *vec, int n)
{
  int x;
  float sum;

  sum = 0.0;
  for (x = 0; x < n; x++) sum += vec[x];
  if (sum != 0.0)
    for (x = 0; x < n; x++) vec[x] /= sum;
  else
    { return 0; }  /* return 0 if all values were zero */
  return 1;
}

/* Helper function to print double vector */
void print_double_vector(const char *label, double *vec, int n)
{
  printf("%s: ", label);
  for (int i = 0; i < n; i++) {
    printf("%.15f", vec[i]);
    if (i < n - 1) printf(", ");
  }
  printf("\n");
}

/* Helper function to print float vector */
void print_float_vector(const char *label, float *vec, int n)
{
  printf("%s: ", label);
  for (int i = 0; i < n; i++) {
    printf("%.8f", vec[i]);
    if (i < n - 1) printf(", ");
  }
  printf("\n");
}

/* Helper function to sum vector */
double sum_double_vector(double *vec, int n)
{
  double sum = 0.0;
  for (int i = 0; i < n; i++) sum += vec[i];
  return sum;
}

float sum_float_vector(float *vec, int n)
{
  float sum = 0.0;
  for (int i = 0; i < n; i++) sum += vec[i];
  return sum;
}

int main(void)
{
  int result;
  double dsum;
  float fsum;

  printf("==================================================\n");
  printf("DNorm and FNorm Vector Normalization Test Program\n");
  printf("==================================================\n\n");

  /* Test 1: {1.0, 2.0, 3.0, 4.0} */
  printf("TEST 1: {1.0, 2.0, 3.0, 4.0}\n");
  printf("---------\n");

  double d_vec1[] = {1.0, 2.0, 3.0, 4.0};
  float f_vec1[] = {1.0f, 2.0f, 3.0f, 4.0f};

  print_double_vector("Double Before", d_vec1, 4);
  dsum = sum_double_vector(d_vec1, 4);
  printf("Sum before: %.15f\n", dsum);

  result = DNorm(d_vec1, 4);
  printf("DNorm result: %d\n", result);
  print_double_vector("Double After ", d_vec1, 4);
  dsum = sum_double_vector(d_vec1, 4);
  printf("Sum after: %.15f\n\n", dsum);

  print_float_vector("Float Before ", f_vec1, 4);
  fsum = sum_float_vector(f_vec1, 4);
  printf("Sum before: %.8f\n", fsum);

  result = FNorm(f_vec1, 4);
  printf("FNorm result: %d\n", result);
  print_float_vector("Float After  ", f_vec1, 4);
  fsum = sum_float_vector(f_vec1, 4);
  printf("Sum after: %.8f\n\n", fsum);

  /* Test 2: {0.1, 0.2, 0.3, 0.4} */
  printf("TEST 2: {0.1, 0.2, 0.3, 0.4}\n");
  printf("---------\n");

  double d_vec2[] = {0.1, 0.2, 0.3, 0.4};
  float f_vec2[] = {0.1f, 0.2f, 0.3f, 0.4f};

  print_double_vector("Double Before", d_vec2, 4);
  dsum = sum_double_vector(d_vec2, 4);
  printf("Sum before: %.15f\n", dsum);

  result = DNorm(d_vec2, 4);
  printf("DNorm result: %d\n", result);
  print_double_vector("Double After ", d_vec2, 4);
  dsum = sum_double_vector(d_vec2, 4);
  printf("Sum after: %.15f\n\n", dsum);

  print_float_vector("Float Before ", f_vec2, 4);
  fsum = sum_float_vector(f_vec2, 4);
  printf("Sum before: %.8f\n", fsum);

  result = FNorm(f_vec2, 4);
  printf("FNorm result: %d\n", result);
  print_float_vector("Float After  ", f_vec2, 4);
  fsum = sum_float_vector(f_vec2, 4);
  printf("Sum after: %.8f\n\n", fsum);

  /* Test 3: {10.0, 20.0, 30.0, 40.0} */
  printf("TEST 3: {10.0, 20.0, 30.0, 40.0}\n");
  printf("---------\n");

  double d_vec3[] = {10.0, 20.0, 30.0, 40.0};
  float f_vec3[] = {10.0f, 20.0f, 30.0f, 40.0f};

  print_double_vector("Double Before", d_vec3, 4);
  dsum = sum_double_vector(d_vec3, 4);
  printf("Sum before: %.15f\n", dsum);

  result = DNorm(d_vec3, 4);
  printf("DNorm result: %d\n", result);
  print_double_vector("Double After ", d_vec3, 4);
  dsum = sum_double_vector(d_vec3, 4);
  printf("Sum after: %.15f\n\n", dsum);

  print_float_vector("Float Before ", f_vec3, 4);
  fsum = sum_float_vector(f_vec3, 4);
  printf("Sum before: %.8f\n", fsum);

  result = FNorm(f_vec3, 4);
  printf("FNorm result: %d\n", result);
  print_float_vector("Float After  ", f_vec3, 4);
  fsum = sum_float_vector(f_vec3, 4);
  printf("Sum after: %.8f\n\n", fsum);

  /* Test 4: {0.0, 0.0, 0.0, 0.0} - edge case */
  printf("TEST 4: {0.0, 0.0, 0.0, 0.0} - Edge Case (All Zeros)\n");
  printf("---------\n");

  double d_vec4[] = {0.0, 0.0, 0.0, 0.0};
  float f_vec4[] = {0.0f, 0.0f, 0.0f, 0.0f};

  print_double_vector("Double Before", d_vec4, 4);
  dsum = sum_double_vector(d_vec4, 4);
  printf("Sum before: %.15f\n", dsum);

  result = DNorm(d_vec4, 4);
  printf("DNorm result: %d (0 = division by zero error)\n", result);
  print_double_vector("Double After ", d_vec4, 4);
  dsum = sum_double_vector(d_vec4, 4);
  printf("Sum after: %.15f\n\n", dsum);

  print_float_vector("Float Before ", f_vec4, 4);
  fsum = sum_float_vector(f_vec4, 4);
  printf("Sum before: %.8f\n", fsum);

  result = FNorm(f_vec4, 4);
  printf("FNorm result: %d (0 = division by zero error)\n", result);
  print_float_vector("Float After  ", f_vec4, 4);
  fsum = sum_float_vector(f_vec4, 4);
  printf("Sum after: %.8f\n\n", fsum);

  printf("==================================================\n");
  printf("All tests completed.\n");
  printf("==================================================\n");

  return 0;
}
