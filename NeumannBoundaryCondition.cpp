

#include <mkl.h>
#include "NeumannBoundaryCondition.h"


void FaceNeumannBoundaryCondition::calcLoadVector(Mesh3D* mesh, double* F)
{
	std::vector<int> Faces;
	mesh->findFaces(Faces, faceStart, faceEnd);
	
	const int &nip = pIntegrator->NIP;							// ��������ֵ����
	double integratePoint[2];

	const int Dim = 3;
	double MiddlePoint[Dim];									// ����е�

	int node0, node1;
	int ElementId = 0;
	double* FVC = nullptr;										// Face Vertex Coordinates

	double* globalPoint;										// ���ϵĻ��ֵ�ӳ��Ϊȫ������
	double* localPoint;											// ���ϵĻ��ֵ�ӳ��Ϊ�ֲ�����
	double loadValue[Dim];										// �û��ֵ��ϵ��غ�ֵ
	double* ShapeN;												// �������ڵ�Ԫ���κ�������
	double detJ;												// ����ſ˱Ⱦ�������ʽ
	double* loadVector = new double[mesh->NCDof];				// �õ�Ԫ���غ�����
	int* LM;													// �õ�Ԫ����ϵ����

	double* GaussPoints = pIntegrator->getGaussCoordinates();	// ���ָ�벻��ɾ,����������������ͷ�
	double* GaussWeights = pIntegrator->getGaussWeights();		// ���ָ�벻��ɾ,����������������ͷ�

	for (size_t iFace = 0; iFace < Faces.size(); iFace++)
	{
		memset(loadVector, 0, mesh->NCDof * sizeof(double));								// ����Ԫ�غ����� ����Ϊ0

		detJ = mesh->calcFaceDetJ(Faces[iFace]);											// ���������ſ˱Ⱦ�������ʽ
		FVC = mesh->getFaceVertexCoords(Faces[iFace]);										// �õ��غ�����ĸ����������

		MiddlePoint[0] = (FVC[0] + FVC[6]) / 2.0;											// ���������������е�
		MiddlePoint[1] = (FVC[1] + FVC[7]) / 2.0;
		MiddlePoint[2] = (FVC[2] + FVC[8]) / 2.0;		

		ElementId = mesh->findElementByPoint(MiddlePoint);									// �����е��ҵ������ڵĵ�Ԫ
		Hexahedral* hex = &mesh->elementVector_[ElementId];									//
		LM = hex->getLocationMatrix();														// �õ��õ�Ԫ����ϵ����

		double alpha, beta = 1.0;
		int m = mesh->NCDof, n = 1, k = Dim;
		int lda = mesh->NCDof, ldb = 1, ldc = 1;

		for (int jj = 0; jj < nip; jj++)													//
		{
			for (int ii = 0; ii < nip; ii++)
			{
				integratePoint[0] = GaussPoints[ii];										//
				integratePoint[1] = GaussPoints[jj];

				globalPoint = mesh->mapLocalToGlobalForFace(integratePoint, FVC);			// ��2ά��˹��ת��Ϊȫ�����ꣻ
				LoadFunction(globalPoint, loadValue);										// ����õ���غ�ֵ

				localPoint = hex->mapGlobalToLocal(globalPoint);							//��ȫ������ת��Ϊ��Ԫ�ڵľֲ�����
				ShapeN = hex->calcShapeN(localPoint);										//���ݾֲ��������κ��������ֵ
				alpha = detJ * GaussWeights[ii] * GaussWeights[jj];
				cblas_dgemm(CblasRowMajor, CblasTrans, CblasNoTrans, m, n, k, alpha, ShapeN, lda, loadValue, ldb, beta, loadVector, ldc);


				delete[]globalPoint;
				delete[]localPoint;
			}
		}

		mesh->assembleToFv(LoadVector, F, LM);
		delete[]LM;

	}
	delete[]loadVector;
}







