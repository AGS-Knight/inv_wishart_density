inv_wishart_density
===================

A C function utilizing GSL which calculates the pdf of the inverse Wishart.

For a toy example, just make the program.  

    ---> ./invwishpdf 3 3 

    Inverse Wishart matrix
    00.18 -0.09 -0.01 
    -0.09 00.28 00.02 
    -0.01 00.02 00.14 

    Scale matrix
    01.00 00.00 00.00 
    00.00 01.00 00.00 
    00.00 00.00 01.00 

    Wishart matrix
    06.74 02.20 00.33 
    02.20 04.28 -0.32 
    00.33 -0.32 07.00 

    Density
    296.446395

Included are a matrix trace function, a matrix determinant function, a matrix inversion function, and a multivariate gamma function.
