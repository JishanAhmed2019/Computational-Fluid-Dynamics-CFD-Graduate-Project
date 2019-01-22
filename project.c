/* ********************************************************************
 *
 * ***************************************************************** */
# include "dg1D.h"

//double minmod(double, double, double);
//double minmod2(double, double, double, double, double);

/* ****************************************************************** */
void project (Mesh *pMesh) {
/* ****************************************************************** */

	Cell *cell = pMesh->cell;
   int i, j, k;
   double u[NVAR], ux[NVAR], uxb[NVAR], dul[NVAR], dur[NVAR], R[NVAR][NVAR],
      Ri[NVAR][NVAR], fact;

   fact = sqrt(3.0);
   
	int Mfact = 10.;

   for (i = 1; i < nCells - 1; i++) {
      for (j = 0; j < NVAR; j++) {

         dul[j] = cell[i].Un[j][0] - cell[i - 1].Un[j][0];
         dur[j] = cell[i + 1].Un[j][0] - cell[i].Un[j][0];
         u[j] = cell[i].Un[j][0];
         ux[j] = fact * cell[i].Un[j][1];
      }

      eigenMat (u, R, Ri);
      mxv (Ri, ux);
      mxv (Ri, dul);
      mxv (Ri, dur);

      for (j = 0; j < NVAR; j++) {

         if(fabs(ux[j]) <= Mfact * cell[i].h * cell[i].h)
            uxb[j] = ux[j];
         else
            uxb[j] = minmod(ux[j], dul[j], dur[j]);
      }

      mxv (R, uxb);

      for (j = 0; j < NVAR; j++) {
         uxb[j] = uxb[j] / fact;
         if (cell[i].Un[j][1] != uxb[j]) {
            cell[i].Un[j][1] = uxb[j];

            for (k = 2; k < cell[i].p; k++)
               cell[i].Un[j][k] = 0.0;
         }

      }
   }
}
/* ****************************************************************** */
double minmod (double a, double b, double c) {
/* ****************************************************************** */
   double sgn, m;

   if(a * b <= 0.0 || b * c <= 0.0)
      return 0.0;

   sgn = (a > 0.0) ? 1.0 : -1.0;
   a = fabs(a);
   b = fabs(b);
   c = fabs(c);
   m = (a < b) ? a : b;
   m = (c < m) ? c : m;
   return sgn * m;

}
/* ****************************************************************** */
void eigenMat (double *U, double R[][3], double Ri[][3]) {
/* ****************************************************************** */
  double d, v, p, c, h, g1, g2;
	g1 = GAMMA - 1.0;
	g2 = g1 / 2.0;

	d = U[0];
	v = U[1] / d;
	p = (GAMMA - 1.0) * (U[2] - 0.5 * d * v * v);
	c = sqrt(GAMMA * p / d);
	h = c * c / g1 + 0.5 * v * v;
	double alpha = g2 /c /c;

	Ri[0][0] = 0.5*(alpha*v*v+v/c);
	Ri[1][0] = 1.-alpha*v*v;
	Ri[2][0] = 0.5*(alpha*v*v-v/c);

	Ri[0][1] = -alpha*v-0.5/c;
	Ri[1][1] = 2.*alpha*v;
	Ri[2][1] = -alpha*v+0.5/c;

	Ri[0][2] = alpha;
	Ri[1][2] = -2.*alpha;
	Ri[2][2] = alpha;

	R[0][0] = 1.;
	R[1][0] = v-c;
	R[2][0] = h-v*c;

	R[0][1] = 1.;
	R[1][1] = v;
	R[2][1] = 0.5*v*v;

	R[0][2] = 1.;
	R[1][2] = v+c;
	R[2][2] = h+v*c;

	return;
}
/* ****************************************************************** */
void mxv (double R[][3], double *U) {
/* ****************************************************************** */
   int i, j;
   double Ut[NVAR];

   for(i = 0; i < NVAR; i++)
      Ut[i] = U[i];

   for(i = 0; i < NVAR; i++) {
      U[i] = 0.0;
      for(j = 0; j < NVAR; j++)
         U[i] += R[i][j] * Ut[j];
   }
}
/* ****************************************************************** */
