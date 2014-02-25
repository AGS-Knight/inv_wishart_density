double matrix_determ(gsl_matrix *X);
gsl_matrix fill_matrix(gsl_matrix *X);
double matrix_trace(gsl_matrix *X);
gsl_matrix inv_matrix(gsl_matrix *X, gsl_matrix *inv);
double mv_gamma(double a, double d);
void print_matrix(gsl_matrix *X);
double iwishpdf(gsl_matrix *X, gsl_matrix *Scale, gsl_matrix *inv, double dof);