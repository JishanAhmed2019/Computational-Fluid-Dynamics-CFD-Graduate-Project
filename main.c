/* *********************************************************************
 * 
 * 
 ******************************************************************** */

# include <stdio.h>
# include <string.h>
# include <stdlib.h>
# include <math.h>

# include "dg1D.h"

/* ****************************************************************** */
int main (int argc, char *argv[]) {
/* ****************************************************************** */

	Mesh mesh;

	initIntegration ();	
	defineMesh (&mesh);

	cfl = cfl / mesh.pOrd;
	double time = 0.0;
	int iter = 0;
	
	int rk, RK = 3;

	printf ("cfl=%g timeOut=%g\n",cfl,timeOut);

	printf("Beginning of iterations ...\n");
	while (time < timeOut) {

		getUo (&mesh);
		timeStep (&mesh);

		if (time + mesh.dt > timeOut)
			mesh.dt = timeOut - time;

		for (rk = 0; rk < RK; rk++) {
			Flux (&mesh);
			update (rk, &mesh);
			project (&mesh);
		}

		time += mesh.dt;
		++iter;
      printf("%8d \t%g\t%g\n", iter, mesh.dt, time);      
	}
	writeSolution (&mesh);
		
	return EXIT_SUCCESS;
}
/* ****************************************************************** */
void getUo (Mesh *pMesh) {
/* ****************************************************************** */
	int i, j, k;
	Cell *cell = pMesh->cell;
	
	for (i = 0; i < nCells; i++)
		for (j = 0; j < NVAR; j++)
			for (k = 0; k < cell[i].p; k++)
				cell[i].Uo[j][k] = cell[i].Un[j][k];
	return;
}
/* ****************************************************************** */
void timeStep (Mesh *pMesh) {
/* ****************************************************************** */
	int i;
	double d, u, p, c, t;
	Cell *cell = pMesh->cell;

	double dt = 1.0e20;

	for (i = 0; i < nCells; i++) {
		d = cell[i].Un[0][0];
		u = cell[i].Un[1][0] / d;
		p = (GAMMA - 1.0) * (cell[i].Un[2][0] - 0.5 * d * u * u);
		c = sqrt(GAMMA * p / d);
		t = cell[i].h / (fabs(u) + c);
		dt = (t < dt) ? t : dt;
	}
	pMesh->dt = cfl * dt;
	
	return;
}
/* ****************************************************************** */
void update (int rk, Mesh *pMesh) {
/* ****************************************************************** */
	int i, j, k;
	double dtdx;
	double dt = pMesh->dt;
	
	Cell *cell = pMesh->cell;
	
	double ark[3],brk[3];
	
	ark[0] = 0.0;
	ark[1] = 3.0 / 4.0;
	ark[2] = 1.0 / 3.0;

	brk[0] = 1.0;
	brk[1] = 1.0 / 4.0;
	brk[2] = 2.0 / 3.0;
	
	for (i = 0; i < nCells; i++) {
	
		dtdx = dt / cell[i].h;
		for (j = 0; j < NVAR; j++)

			for (k = 0; k < cell[i].p; k++)
			
				cell[i].Un[j][k] =
               ark[rk] * cell[i].Uo[j][k] + brk[rk] * (cell[i].Un[j][k] -
                                                       dtdx * cell[i].Re[j][k]);
	}
	return;
}
/* ****************************************************************** */



