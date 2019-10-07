#ifndef DirichletBoundaryCondition_H
#define DirichletBoundaryCondition_H

#include "Mesh3D.h"

/*------------------------------------------------------------------------------------------------------------------
	��������
------------------------------------------------------------------------------------------------------------------*/
class StrongPenaltyAlgorithm
{
public:
	StrongPenaltyAlgorithm(double Penalty) :penaltyValue(Penalty) {	}

	~StrongPenaltyAlgorithm() { }

	void modifyLinearSystem(int dofId, double prescribedValues, double* K, int* Kdiag, double* F)
	{
		K[Kdiag[dofId]] += penaltyValue;
		F[dofId] += penaltyValue * prescribedValues;
	}


public:
	double penaltyValue;
};

/*------------------------------------------------------------------------------------------------------------------
	����߽��������� 
------------------------------------------------------------------------------------------------------------------*/
class DirichletBoundaryCondition
{
public:
	DirichletBoundaryCondition(double* pre, int* dir, double Penalty):PenaltyAlgorithm(Penalty)
	{
		prescribedValue[0] = pre[0];
		prescribedValue[1] = pre[1];
		prescribedValue[2] = pre[2];

		direction[0] = dir[0];
		direction[1] = dir[1];
		direction[2] = dir[2];
	}

	virtual ~DirichletBoundaryCondition()	{ }	

	virtual void modifyLinearSystem(Mesh3D* mesh, double* K, int* Kdiag, double* F) = 0;

	virtual void constrainFaces(Mesh3D* mesh, std::vector<int>& faces, double* K, int* Kdiag, double* F, const int& Dim);

	virtual void constrainEdges(Mesh3D* mesh, std::vector<int>& edges, double* K, int* Kdiag, double* F, const int& Dim);

	virtual void constrainNodes(Mesh3D* mesh, std::vector<int>& nodes, double* K, int* Kdiag, double* F, const int& Dim);

public:
	double prescribedValue[3];
	int direction[3];
	StrongPenaltyAlgorithm PenaltyAlgorithm;
};







/*------------------------------------------------------------------------------------------------------------------
	�� λ�Ʊ߽�������
------------------------------------------------------------------------------------------------------------------*/

class FaceDirichletBoundaryCondition :public DirichletBoundaryCondition
{
public:
	FaceDirichletBoundaryCondition(double*face_start, double*face_end, double* pre, int* dir, double Penalty)
		: DirichletBoundaryCondition(pre, dir, Penalty)
	{
		faceStart[0] = face_start[0];
		faceStart[1] = face_start[1];
		faceStart[2] = face_start[2];

		faceEnd[0] = face_end[0];
		faceEnd[1] = face_end[1];
		faceEnd[2] = face_end[2];
	}

	~FaceDirichletBoundaryCondition()	{ }


	virtual void modifyLinearSystem(Mesh3D* mesh, double* K, int* Kdiag, double* F);

public:
	double faceStart[3];
	double faceEnd[3];
};






/*------------------------------------------------------------------------------------------------------------------
	�� λ�Ʊ߽�������
------------------------------------------------------------------------------------------------------------------*/
class EdgeDirichletBoundaryCondition :public DirichletBoundaryCondition
{
public:
	EdgeDirichletBoundaryCondition(double*line_start, double*line_end, double* pre, int* dir, double Penalty)
		: DirichletBoundaryCondition(pre, dir, Penalty)
	{
		lineStart[0] = line_start[0];
		lineStart[1] = line_start[1];
		lineStart[2] = line_start[2];

		lineEnd[0] = line_end[0];
		lineEnd[1] = line_end[1];
		lineEnd[2] = line_end[2];
	}

	~EdgeDirichletBoundaryCondition()	{ }


	virtual void modifyLinearSystem(Mesh3D* mesh, double* K, int* Kdiag, double* F);

public:
	double lineStart[3];
	double lineEnd[3];
};




/*------------------------------------------------------------------------------------------------------------------
	�� λ�Ʊ߽�������
------------------------------------------------------------------------------------------------------------------*/
class NodeDirichletBoundaryCondition :public DirichletBoundaryCondition
{
public:
	NodeDirichletBoundaryCondition(double*point, double* pre, int* dir, double Penalty)
		: DirichletBoundaryCondition(pre, dir, Penalty)
	{
		position[0] = point[0];
		position[1] = point[1];
		position[2] = point[2];
	}

	~NodeDirichletBoundaryCondition()	{ }


	virtual void modifyLinearSystem(Mesh3D* mesh, double* K, int* Kdiag, double* F);

public:
	double position[3];

};









#endif