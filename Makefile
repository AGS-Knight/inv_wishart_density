CC = gcc
CFLAGS = -Wall -lm -lgsl -lgslcblas -llapack
	
InvWishPDF: invwishpdf.c
	$(CC) invwishpdf.c  -O2 -o invwishpdf $(CFLAGS)