#ifndef FiniteCellManager_H
#define FiniteCellManager_H



#include "Mesh3D.h"
#include "EmbeddedDomain.h"
#include "Material.h"
#include "GaussIntegrator.h"
#include "NeumannBoundaryCondition.h"

class FiniteCellManager
{
public:
	FiniteCellManager() :mesh(nullptr), domain(nullptr), material(nullptr), integrator(nullptr)
	{

	}

	virtual ~FiniteCellManager()
	{
		delete[]mesh;
		delete[]domain;
		delete[]material;
		delete[]integrator;
	}


	void initialize();

	void generateMesh();

	void assemble_Kc_To_Kv(double *Kc, double *Kv, int *Kdiag, std::vector<int> &Columns, int *LM, const int NCDof);

	void calcMHTandColumns(int *MHT, std::vector<int>&columnVec, const int& PD);

	void getShareCellLM(int* shareLM, const std::vector<int> &cellId, const int cellDofs);

	void calcKdiag(int *Kdiag, int *MHT, const int neq);

	void calcKsparse(double* Kv, int* Kdiag, std::vector<int>& columns);

	void printStiffness(const double* Kv, const int* Kdiag, const int NEQ) const;


public:
	int NWK;								// Number of wide of K×Ü¸ÕµÄ¿í¶È
	int spaceTreeDepth;
	//int polynomialDegree;
	Mesh3D* mesh;
	EmbeddedDomain* domain;
	AbsMaterial* material;
	GaussIntegrator* integrator;


};










#endif