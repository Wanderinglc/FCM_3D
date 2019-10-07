#ifndef NeumannBoundaryCondition_H
#define NeumannBoundaryCondition_H

#include "GaussIntegrator.h"
#include "Mesh3D.h"

class LoadFunctionObject
{
public:
	LoadFunctionObject(double* Value)
	{
		value[0] = Value[0];
		value[1] = Value[1];
		value[2] = Value[2];
	}

	~LoadFunctionObject()	{ }

	double* operator() (const double* coords)
	{
		double* loadValue = new double[3];
		loadValue[0] = value[0];
		loadValue[1] = value[1];
		loadValue[2] = value[2];
		return loadValue;
	}

public:
	double value[3];
};



class NeumannBoundaryCondition
{
public:

	NeumannBoundaryCondition(LoadFunctionObject* loadFuc, GaussIntegrator* integrator):
		pIntegrator(integrator), pLoadFuction(loadFuc)
	{

	}

	virtual ~NeumannBoundaryCondition()	 { }

	virtual void calcLoadVector(Mesh3D* mesh, double* F) = 0;


public:
	LoadFunctionObject* pLoadFuction;
	GaussIntegrator* pIntegrator;
};





/*---------------------------------------------------------
	面 载荷边界条件类
---------------------------------------------------------*/
class FaceNeumannBoundaryCondition:public NeumannBoundaryCondition
{
public:
	FaceNeumannBoundaryCondition(double*face_start, double*face_end, LoadFunctionObject* loadFuc, GaussIntegrator* integrator):
		NeumannBoundaryCondition(loadFuc, integrator)
	{
		faceStart[0] = face_start[0];
		faceStart[1] = face_start[1];
		faceStart[2] = face_start[2];

		faceEnd[0] = face_end[0];
		faceEnd[1] = face_end[1];
		faceEnd[2] = face_end[2];
	}

	~FaceNeumannBoundaryCondition()	 { }

	virtual void calcLoadVector(Mesh3D* mesh, double* F);


public:
	double faceStart[3];
	double faceEnd[3];
};





#endif