
#include <mkl.h>

#include "QuasiStaticAnalysis.h"
#include "FiniteCellManager.h"
#include "DirichletBoundaryCondition.h"
#include "NeumannBoundaryCondition.h"


void QuasiStaticAnalysis::execute()
{
	//------------------------------------------------------------------------------ ==============��������==============
	const int &NEQ = fcManager->mesh->totalDofs;											// �ܷ�����
	printf("NEQ =  %d\n", NEQ);
	int *MHT = new int[NEQ];																// �п����飬�ܸ���ÿһ�еĿ��
	int* Kdiag = new int[NEQ + 1];															// �Խ�Ԫ��һά�ܸ��е�λ��

	//------------------------------------------------------------------------------ ==============�����������==============
	std::vector<int> columns;																// ����Ԫ�����ܸ��е����±�
	fcManager->calcMHTandColumns(MHT, columns, fcManager->mesh->polynomialDegree_);			// �γ� �п����� �����±�����
	fcManager->calcKdiag(Kdiag, MHT, NEQ);													// ����Խ�Ԫλ������
	delete[]MHT;  MHT = nullptr;															// �ͷſռ�

	//------------------------------------------------------------------------------ ==============�����ܸ�==============
	double* Kv = new double[fcManager->NWK]();												// �ܸ�����
	fcManager->calcKsparse(Kv, Kdiag, columns);

	//------------------------------------------------------------------------------ ==============�����غ�����==============
	double* Fv = new double[NEQ]();															// ���غ�����
	for (size_t nBC = 0; nBC < NeumannBC.size(); nBC++)
	{
		NeumannBC[nBC]->calcLoadVector(fcManager->mesh, Fv);
	}

	//------------------------------------------------------------------------------ ==============ʩ��λ�Ʊ߽�����==============
	for (size_t dBC = 0; dBC < DirichletBC.size(); dBC++)
	{
		DirichletBC[dBC]->modifyLinearSystem(fcManager->mesh, Kv, Kdiag, Fv);
	}

	//------------------------------------------------------------------------------ ==============�������==============
	double* Uv = new double[NEQ]();															// ��λ������

	pardisoSolver(Kv, Kdiag, &columns[0], Fv, Uv, NEQ);
	
	//------------------------------------------------------------------------------ ==============���λ��==============
	printf("\nThe solution vector is: \n");
	for (int i = 0; i < NEQ; i++)
	{
		printf("U[%d] = %f\n", i + 1, Uv[i]);
	}




	//------------------------------------------------------------------------------ ==============��λ�Ʒֲ�����Ԫ��==============
	



	//------------------------------------------------------------------------------ ==============�ռ��ͷ�==============
	delete[]Kdiag;	Kdiag = nullptr;
	delete[]Kv;		Kv = nullptr;
	delete[]Fv;		Fv = nullptr;
	delete[]Uv;		Uv = nullptr;
}





void QuasiStaticAnalysis::pardisoSolver(double *Kv, int *Kdiag, int *Columns, double *Lv, double *Uv, const int &Neq)
{
	int n = Neq;
	//int mtype = 2;		/* Real and symmetric positive definite */
	int mtype = -2;			/* Real symmetric matrix */
	int nrhs = 1;			/* Number of right hand sides. */
	void *pt[64];			/* Internal solver memory pointer pt, */
	/* Pardiso control parameters. */
	MKL_INT iparm[64];
	MKL_INT maxfct, mnum, phase, error, msglvl;
	/* Auxiliary variables. */
	MKL_INT i;
	double ddum;          /* Double dummy */
	MKL_INT idum;         /* Integer dummy. */
/* -------------------------------------*/
/* .. Setup Pardiso control parameters. */
/* -------------------------------------*/
	for (i = 0; i < 64; i++)
	{
		iparm[i] = 0;
	}
	iparm[0] = 1;         /* No solver default */
	iparm[1] = 2;         /* Fill-in reordering from METIS */
	iparm[3] = 0;         /* No iterative-direct algorithm */
	iparm[4] = 0;         /* No user fill-in reducing permutation */
	iparm[5] = 0;         /* Write solution into x */
	iparm[7] = 0;         /* Max numbers of iterative refinement steps */
	iparm[9] = 7;         /* Perturb the pivot elements with 1E-8 */ //��Ԫ�Ŷ��������ڶԳƲ�������
	//iparm[10] = 1;        /* Use nonsymmetric permutation and scaling MPS */
	iparm[10] = 0;        /* Default for symmetric indefinite matrices */
	iparm[12] = 0;        /* Maximum weighted matching algorithm is switched-off (default for symmetric). Try iparm[12] = 1 in case of inappropriate accuracy */
	iparm[13] = 0;        /* Output: Number of perturbed pivots */
	iparm[17] = -1;       /* Output: Number of nonzeros in the factor LU */
	iparm[18] = -1;       /* Output: Mflops for LU factorization */
	iparm[19] = 0;        /* Output: Numbers of CG Iterations */
	iparm[26] = 1;
	iparm[34] = 1;        /* PARDISO use C-style indexing for ia and ja arrays */
	maxfct = 1;           /* Maximum number of numerical factorizations. */
	mnum = 1;			  /* Which factorization to use. */
	msglvl = 0;           /* Not Print statistical information in file */
	error = 0;            /* Initialize error flag */
/* ----------------------------------------------------------------*/
/* .. Initialize the internal solver memory pointer. This is only  */
/*   necessary for the FIRST call of the PARDISO solver.           */
/* ----------------------------------------------------------------*/
	for (i = 0; i < 64; i++)
	{
		pt[i] = 0;
	}
	/* --------------------------------------------------------------------*/
	/* .. Reordering and Symbolic Factorization. This step also allocates  */
	/*    all memory that is necessary for the factorization.              */
	/* --------------------------------------------------------------------*/
	phase = 11;
	PARDISO(pt, &maxfct, &mnum, &mtype, &phase,
		&n, Kv, Kdiag, Columns, &idum, &nrhs, iparm, &msglvl, &ddum, &ddum, &error);
	if (error != 0)
	{
		printf("\nERROR during symbolic factorization: %d", error);
		exit(1);
	}
	printf("\nReordering completed ... ");
	printf("\nNumber of nonzeros in factors = %d", iparm[17]);
	printf("\nNumber of factorization MFLOPS = %d", iparm[18]);
	/* ----------------------------*/
	/* .. Numerical factorization. */
	/* ----------------------------*/
	phase = 22;
	PARDISO(pt, &maxfct, &mnum, &mtype, &phase,
		&n, Kv, Kdiag, Columns, &idum, &nrhs, iparm, &msglvl, &ddum, &ddum, &error);
	if (error != 0)
	{
		printf("\nERROR during numerical factorization: %d", error);
		exit(2);
	}
	printf("\nFactorization completed ... ");
	/* -----------------------------------------------*/
	/* .. Back substitution and iterative refinement. */
	/* -----------------------------------------------*/
	phase = 33;
	iparm[7] = 2;         /* Max numbers of iterative refinement steps. */
	PARDISO(pt, &maxfct, &mnum, &mtype, &phase,
		&n, Kv, Kdiag, Columns, &idum, &nrhs, iparm, &msglvl, Lv, Uv, &error);
	if (error != 0)
	{
		printf("\nERROR during solution: %d", error);
		exit(3);
	}
	printf("\nSolve completed ... \n");

	/* --------------------------------------*/
	/* .. Termination and release of memory. */
	/* --------------------------------------*/
	phase = 0;           /* Release internal memory. */

	PARDISO(pt, &maxfct, &mnum, &mtype, &phase,
		&n, &ddum, Kdiag, Columns, &idum, &nrhs,
		iparm, &msglvl, &ddum, &ddum, &error);

}













