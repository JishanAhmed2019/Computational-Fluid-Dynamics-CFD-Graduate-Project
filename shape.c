/* ********************************************************************
 *
 * ***************************************************************** */
# include "dg1D.h"

/* ****************************************************************** */
double shapeF (double x, Cell *cell, int nshape) {
/* ****************************************************************** */
	double f;
	
	x = 2.0 * (x - cell->x) / cell->h;
	
	f = sqrt(2.0 * nshape + 1);
	
	return (f * legendreF (x, nshape));
}
/* ****************************************************************** */
double shapeDF (double x, Cell *cell, int nshape) {
/* ****************************************************************** */

	double f;

	x = (x - cell->x) / cell->h;
	
	f = 2.0 * sqrt(2.0 * nshape + 1) / cell->h;
	
	return (f * legendreDF (x, nshape));
}
/* ****************************************************************** */
double legendreF (double x, int n) {
/* ****************************************************************** */
	switch (n) {
		case 0:
			return 1.0;
		case 1:
			return x;
		case 2:
			return 0.5 * (3.0 * x * x - 1.0);
		default:
			printf ("Legendre: Unknown Legendre function no = %d\n", n);
			exit(0);
	}
}
/* ****************************************************************** */
double legendreDF (double x, int n) {
/* ****************************************************************** */
	switch (n) {
		case 0:
			return 0.0;
		case 1:
			return 1.0;
		case 2:
			return 3.0 * x;
		default:
			printf("LegendreDeriv: Unknown Legendre function no = %d\n", n);
			exit(0);
	}
}
/* ****************************************************************** */




