struct mymatrix3 findQmatrix(double sigma,double beta,double rho);
struct mymatrix3 matrix_product(struct mymatrix3 A,struct mymatrix3 B);
struct mymatrix3 transpose(struct mymatrix3 A);
struct myvector matrix_vector(struct mymatrix3 A,struct myvector v);
struct mymatrix3 matrix_sum(struct mymatrix3 A,struct mymatrix3 B);
void print_matrix(struct mymatrix3 U,char *name );