#ifndef GASTONLD
#define GASTONLD
double LD_colxx(matrix4 & A, double mu1, double mu2, double v, size_t x1, size_t x2);
void LD_col(matrix4 & A, bar & mu, bar & sd, size_t c1, size_t c2, lou & LD);
void LD_chunk(matrix4 & A, bar & mu, bar & sd, int c1, int c2, int d1, int d2, lou & LD);
#endif
