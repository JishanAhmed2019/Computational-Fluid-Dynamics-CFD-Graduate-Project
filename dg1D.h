# ifndef DG1D_H
# define DG1D_H

/* ************************************************** */

/* ************************************************** */

# include <math.h>
# include <stdlib.h>
# include <string.h>
# include <stdio.h>
# include <assert.h>

# define GAMMA 1.4
# define CFL 0.8
# define NVAR 3

typedef struct sCell {
	double x, xl, xr, h, *xg;
	int p, ng;
	double **Un, **Uo, **Re;
}Cell;

typedef struct sMesh{
	double 		minX;
	double 		maxX;
//  int 			nx;
  int				nNodes;
	double 		dx;
	double 		time, dt;
	int				pOrd;
	Cell			*cell;
}Mesh;

/* xg = Gauss integration points in [-1,+1]
 * wg = corresponding weights
 */
double xg[3][3], wg[3][3];

double d_left, u_left, p_left;
double d_right, u_right, p_right;

int nCells;
double cfl, timeOut, diaph;

void getUo (Mesh *pMesh);
void timeStep (Mesh *pMesh);
void update (int rk, Mesh *pMesh);

void integrationPoints (Cell *cell);
void initIntegration ();
void readData (Mesh *pMesh);
void initCondEuler (double x, double *U);

void Flux (Mesh *pMesh);
void EulerFlux (double *U, double *flux);
void LFFlux (double *Ul, double *Ur, double *flux);
void Uphi (Cell *cell, double x, double *U);

double shapeF (double x, Cell *cell, int nshape);
double shapeDF (double x, Cell *cell, int nshape);

double legendreF (double, int);
double legendreDF (double, int);

void defineMesh (Mesh *pMesh);
void initVars (Mesh *pMesh);
void writeSolution (Mesh *pMesh);

double phi (int i);
double dphi (int i);

void calcFlux (double* u, double* flux);
void timeStep (Mesh *pMesh);
void calcResidual(Mesh *pMesh);

void initProblem (Mesh *pMesh);
void initCond (double *u, double x);

void project (Mesh *pMesh);
double minmod (double a, double b, double c);
void eigenMat (double *U, double R[][3], double Ri[][3]);
void mxv (double R[][3], double *U);


# endif
