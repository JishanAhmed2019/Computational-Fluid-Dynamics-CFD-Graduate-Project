/* ********************************************************************
 *
 * ***************************************************************** */
# include "dg1D.h"

/* ****************************************************************** */
void readData (Mesh *pMesh) {
/* ****************************************************************** */
	FILE *fp;
	char dummy[100];
	
	fp = fopen("data.i", "r");
	if(fp == NULL) {
		printf("Error: Could not open data.i\n");
		exit(0);
	}
   
	int pOrd;
	double xmin, xmax;

	fscanf(fp, "%s%lf", dummy, &cfl);
	fscanf(fp, "%s%lf", dummy, &timeOut);
	fscanf(fp, "%s%d", dummy, &nCells);
	fscanf(fp, "%s%d", dummy, &pOrd);
	fscanf(fp, "%s%lf%lf", dummy, &xmin, &xmax);
	fscanf(fp, "%s%lf", dummy, &diaph);
	fscanf(fp, "%s%lf%lf%lf", dummy, &d_left, &u_left, &p_left);
	fscanf(fp, "%s%lf%lf%lf", dummy, &d_right, &u_right, &p_right);
	fclose(fp);
   
	pMesh->minX = xmin;
	pMesh->maxX = xmax;
	pMesh->pOrd = pOrd;
   
  return; 
}
/* ****************************************************************** */
void writeSolution (Mesh *pMesh) {
/* ****************************************************************** */

	int i, j;
	double dx, x, U[NVAR], d, u, p, c, m;
	Cell *cell = pMesh->cell;
	   
	char name[9] = "solution";
	FILE *fileOut;
	
	if (!(fileOut = fopen (name,"w+"))) {
		printf ("Error! Unable to open file.\n");
	}

	for (i = 0; i < nCells; i++) {

		dx = 0.5 * cell[i].h;
		
		x = cell[i].xl + dx;

		Uphi (&cell[i], x, U);

		d = U[0];
		u = U[1] / d;
		p = (GAMMA - 1.0) * (U[2] - 0.5 * d * u * u);
		c = sqrt (GAMMA * p / d);
//		m = u / c;

//		fprintf (fileOut, "%f %f %f %f %f\n", x, d, u, p, m);
		fprintf (fileOut, "%f %f %f %f \n", x, d, u, p);
	}
	fclose(fileOut);

	return;
}
/* ****************************************************************** */
