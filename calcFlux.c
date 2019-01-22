/* ********************************************************************
 *
 * ***************************************************************** */
# include "dg1D.h"

double MaxEigVal (double *, double *);

/* ****************************************************************** */
void Flux (Mesh *pMesh) {
/* ****************************************************************** */
	
	int i, j, k, l, cl, cr;
	Cell *cell = pMesh->cell;
	
	double UL[NVAR], UR[NVAR], fl[nCells + 1][NVAR], UG[NVAR], flg[NVAR], vl, vr, vx;

	for (i = 0; i < nCells; i++)
      for (j = 0; j < NVAR; j++)
         for (k = 0; k < cell[i].p; k++)
            cell[i].Re[j][k] = 0.0;

   /* Loop over cell faces and find flux */
	for (i = 0; i <= nCells; i++) {

      cl = (i == 0) ? i : (i - 1);
      cr = (i == nCells) ? (i - 1) : i;

      Uphi (&cell[cl], cell[cl].xr, UL);

      Uphi (&cell[cr], cell[cr].xl, UR);

			LFFlux (UL, UR, fl[i]);
	}

   /* Add interface flux to the cells */
	for (i = 0; i < nCells; i++)

      for (j = 0; j < NVAR; j++)

         for (k = 0; k < cell[i].p; k++) {

            vl = shapeF (cell[i].xl, &cell[i], k);

            vr = shapeF (cell[i].xr, &cell[i], k);

            cell[i].Re[j][k] += fl[i + 1][j] * vr - fl[i][j] * vl;
         }

   /* Flux quadrature */
	for (i = 0; i < nCells; i++)
      for (j = 0; j < cell[i].ng; j++) {

         Uphi (&cell[i], cell[i].xg[j], UG);

         EulerFlux (UG, flg);

         for (k = 0; k < cell[i].p; k++) {

            vx = shapeDF (cell[i].xg[j], &cell[i], k);

            for (l = 0; l < NVAR; l++)
            
               cell[i].Re[l][k] -=
                  0.5 * cell[i].h * flg[l] * vx * wg[cell[i].ng - 1][j];
         }
      }
	return;
}
/* ****************************************************************** */
void EulerFlux (double *U, double *flux) {
/* ****************************************************************** */
	double p;

	p = (GAMMA - 1.0) * (U[2] - 0.5 * U[1] * U[1] / U[0]);

	flux[0] = U[1];
	flux[1] = p + U[1] * U[1] / U[0];
	flux[2] = (U[2] + p) * U[1] / U[0];

	return;
}
/* ****************************************************************** */
void LFFlux (double *Ul, double *Ur, double *flux) {
/* ****************************************************************** */

	int i;
	double Fl[NVAR], Fr[NVAR], lam;

	EulerFlux (Ul, Fl);
	EulerFlux (Ur, Fr);

	lam = MaxEigVal (Ul, Ur);

	for (i = 0; i < NVAR; i++)
		flux[i] = 0.5 * (Fl[i] + Fr[i] - lam * (Ur[i] - Ul[i]));

	return;
}
/* ****************************************************************** */
double MaxEigVal (double *Ul, double *Ur) {
/* ****************************************************************** */
	double dl, ul, pl, al, ll, dr, ur, pr, ar, lr;

	dl = Ul[0];
	ul = Ul[1] / dl;
	pl = (GAMMA - 1.0) * (Ul[2] - 0.5 * dl * ul * ul);
	al = sqrt(GAMMA * pl / dl);
	ll = fabs(ul) + al;

	dr = Ur[0];
	ur = Ur[1] / dr;
	pr = (GAMMA - 1.0) * (Ur[2] - 0.5 * dr * ur * ur);
	ar = sqrt(GAMMA * pr / dr);
	lr = fabs(ur) + ar;

	if (ll > lr)
		return ll;
	else
		return lr;
}
/* ****************************************************************** */
void Uphi (Cell *cell, double x, double *U) {
/* ****************************************************************** */
	int iv, ip;

	for (iv = 0; iv < NVAR; iv++) {
		U[iv] = 0.0;
		
		for(ip = 0; ip < cell->p; ip++)
			U[iv] += cell->Un[iv][ip] * shapeF (x, cell, ip);	
	}
	return;
}
/* ****************************************************************** */













