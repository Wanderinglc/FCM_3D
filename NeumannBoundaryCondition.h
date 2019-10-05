#ifndef NeumannBoundaryCondition_H
#define NeumannBoundaryCondition_H

#include "GaussIntegrator.h"
#include "Mesh3D.h"

class NeumannBoundaryCondition
{
public:
	NeumannBoundaryCondition()
	{

	}

	virtual ~NeumannBoundaryCondition()
	{

	}


	virtual void calcLoadVector(Mesh3D* mesh, double* F) = 0;





public:
	GaussIntegrator* pIntegrator;
};





/*---------------------------------------------------------
	面 载荷边界条件类
---------------------------------------------------------*/
class FaceNeumannBoundaryCondition:public NeumannBoundaryCondition
{
public:
	FaceNeumannBoundaryCondition(double*face_start, double*face_end):NeumannBoundaryCondition()
	{
		faceStart[0] = face_start[0];
		faceStart[1] = face_start[1];
		faceStart[2] = face_start[2];

		faceEnd[0] = face_end[0];
		faceEnd[1] = face_end[1];
		faceEnd[2] = face_end[2];
	}

	~FaceNeumannBoundaryCondition()
	{

	}

	virtual void calcLoadVector(Mesh3D* mesh, double* F);


public:
	double faceStart[3];
	double faceEnd[3];
};







#endif