/* Example of calling a Fortran routine from C */

const char* dgemm_desc = "Fortran dgemm.";

int dgemm_(const char*, const char*, const int*, const int*, const int*,
           const double*, const double*, const int*,
           const double*, const int*,
           const double*, double*, const int*);

void square_dgemm(const int M, double *A, double *B, double *C)
{
    double one = 1.0;
    double zero = 0.0;
    dgemm_("N", "N", &M, &M, &M, &one, A, &M, B, &M, &zero, C, &M);
}
