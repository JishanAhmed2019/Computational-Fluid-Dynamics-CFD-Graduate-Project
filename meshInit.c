/* ********************************************************************
 *
 * ***************************************************************** */
# include "dg1D.h"

/* ****************************************************************** */
void defineMesh (Mesh *pMesh) {
/* ****************************************************************** */

	readData (pMesh);

	pMesh->nNodes = nCells + 1;
	pMesh->dx = (pMesh->maxX - pMesh->minX) / nCells;
	
	int i, j, k, l;

	int pOrd = pMesh->pOrd;
	
	int itgPts = (2*pOrd-1);

	double U[NVAR];
	double phi;

	pMesh->cell = (Cell *) calloc (nCells, sizeof (Cell));	
		if (pMesh->cell == NULL) {
			printf("Init: Could not allocate cell\n");
			exit(0);
	}

	Cell *cell;
	cell = pMesh->cell;

   for (i = 0; i < nCells; i++) {
      cell[i].xl = pMesh->minX + i * pMesh->dx;
      cell[i].xr = cell[i].xl + pMesh->dx;
      cell[i].x = 0.5 * (cell[i].xl + cell[i].xr);
      cell[i].h = cell[i].xr - cell[i].xl;

      cell[i].p = pOrd;

      cell[i].ng = itgPts;
      cell[i].xg = (double *) calloc (cell[i].ng, sizeof (double));
      
      integrationPoints (&cell[i]);

      cell[i].Un = (double **) calloc (NVAR, sizeof (double *));
      cell[i].Uo = (double **) calloc (NVAR, sizeof (double *));
      cell[i].Re = (double **) calloc (NVAR, sizeof (double *));
      for(j = 0; j < NVAR; j++) {   
         cell[i].Un[j] = (double *) calloc (cell[i].p, sizeof (double));
         cell[i].Uo[j] = (double *) calloc (cell[i].p, sizeof (double));
         cell[i].Re[j] = (double *) calloc (cell[i].p, sizeof (double));
      }
   }
	printf ("polynomial order=%d\n", pMesh->pOrd);
	printf ("Domain: [%g,%g] #Cells=%d #Nodes=%d Cell length=%g\n",
					pMesh->minX, pMesh->maxX, nCells, pMesh->nNodes, pMesh->dx);

   /* Set initial condition by L2 projection */
	for (i = 0; i < nCells; i++) {

		for (j = 0; j < NVAR; j++)    
			for (k = 0; k < cell[i].p; k++)
				cell[i].Un[j][k] = 0.0;
            
            
		for (j = 0; j < cell[i].p; j++)
		
			for (k = 0; k < cell[i].ng; k++) {

				initCondEuler (cell[i].xg[k], U);

				phi = shapeF (cell[i].xg[k], &cell[i], j);

				for (l = 0; l < NVAR; l++)
					cell[i].Un[l][j] += 0.5 * U[l] * phi * wg[cell[i].ng - 1][k];
			}
	}
	return;
}
/* ****************************************************************** */
void initIntegration () {
/* ****************************************************************** */
   int n = 3, i, j;

   printf("Calculating Gauss integration points and weights ...\n");

   for (i = 0; i < n; i++)
      for(j = 0; j < n; j++) {
         xg[i][j] = 0.0;
         wg[i][j] = 0.0;
      }

   /* Gauss integration points */
   xg[0][0] = 0.0;

   xg[1][0] = -1.0 / sqrt(3.0);
   xg[1][1] = 1.0 / sqrt(3.0);

   xg[2][0] = -sqrt(15.0) / 5.0;
   xg[2][1] = 0.0;
   xg[2][2] = sqrt(15.0) / 5.0;

   /* Gaussian weights */
   wg[0][0] = 2.0;

   wg[1][0] = 1.0;
   wg[1][1] = 1.0;

   wg[2][0] = 5.0 / 9.0;
   wg[2][1] = 8.0 / 9.0;
   wg[2][2] = 5.0 / 9.0;

	return;
}
/* ****************************************************************** */
void integrationPoints (Cell *cell) {
/* ****************************************************************** */
	double xl, xr;
	int i;
	
	xl = cell->xl;
	xr = cell->xr;

	for (i = 0; i < cell->ng; i++)
		cell->xg[i] = 0.5 * (xl * (1.0 - xg[cell->ng - 1][i]) +
                         xr * (1.0 + xg[cell->ng - 1][i]));
	return;
}
/* ****************************************************************** */
void initCondEuler (double x, double *U) {
/* ****************************************************************** */
	double d, u, p;

	if (x < diaph) {
		d = d_left;
		u = u_left;
		p = p_left;
	}
	else {
		d = d_right;
		u = u_right;
		p = p_right;
	}

	U[0] = d;
	U[1] = d * u;
	U[2] = p / (GAMMA - 1.0) + 0.5 * d * u * u;
   
	return;
}
/* ****************************************************************** */




