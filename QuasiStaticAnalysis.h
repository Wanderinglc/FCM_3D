#ifndef QuasiStaticAnalysis_H
#define QuasiStaticAnalysis_H

#include "FiniteCellManager.h"
#include "DirichletBoundaryCondition.h"
#include "NeumannBoundaryCondition.h"

class QuasiStaticAnalysis
{
public:
	QuasiStaticAnalysis(): fcManager(nullptr), Dim(3) { }
	~QuasiStaticAnalysis()	{ }

	void execute();
	void pardisoSolver(double *Kv, int *Kdiag, int *Columns, double *Lv, double *Uv, const int &Neq);


public:
	int Dim;
	FiniteCellManager* fcManager;
	std::vector<DirichletBoundaryCondition *> DirichletBC;
	std::vector<NeumannBoundaryCondition *> NeumannBC;
};











#endif
