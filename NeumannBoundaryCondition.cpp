

#include <mkl.h>
#include "NeumannBoundaryCondition.h"


void FaceNeumannBoundaryCondition::calcLoadVector(Mesh3D* mesh, double* F)
{
	std::vector<int> Faces;
	mesh->findFaces(Faces, faceStart, faceEnd);
	
	const int &nip = pIntegrator->NIP;							// 单方向积分点个数
	double integratePoint[2];

	const int Dim = 3;
	double MiddlePoint[Dim];									// 面的中点

	int node0, node1;
	int ElementId = 0;
	double* FVC = nullptr;										// Face Vertex Coordinates

	double* globalPoint;										// 面上的积分点映射为全局坐标
	double* localPoint;											// 面上的积分点映射为局部坐标
	double loadValue[Dim];										// 该积分点上的载荷值
	double* ShapeN;												// 该面所在单元的形函数矩阵
	double detJ;												// 面的雅克比矩阵行列式
	double* loadVector = new double[mesh->NCDof];				// 该单元的载荷向量
	int* LM;													// 该单元的联系数组

	double* GaussPoints = pIntegrator->getGaussCoordinates();	// 这个指针不用删,留给类的析构函数释放
	double* GaussWeights = pIntegrator->getGaussWeights();		// 这个指针不用删,留给类的析构函数释放

	for (size_t iFace = 0; iFace < Faces.size(); iFace++)
	{
		memset(loadVector, 0, mesh->NCDof * sizeof(double));								// 将单元载荷向量 重置为0

		detJ = mesh->calcFaceDetJ(Faces[iFace]);											// 计算该面的雅克比矩阵行列式
		FVC = mesh->getFaceVertexCoords(Faces[iFace]);										// 得到载荷面的四个顶点的坐标

		MiddlePoint[0] = (FVC[0] + FVC[6]) / 2.0;											// 根据坐标计算面的中点
		MiddlePoint[1] = (FVC[1] + FVC[7]) / 2.0;
		MiddlePoint[2] = (FVC[2] + FVC[8]) / 2.0;		

		ElementId = mesh->findElementByPoint(MiddlePoint);									// 根据中点找到其所在的单元
		Hexahedral* hex = &mesh->elementVector_[ElementId];									//
		LM = hex->getLocationMatrix();														// 得到该单元的联系数组

		double alpha, beta = 1.0;
		int m = mesh->NCDof, n = 1, k = Dim;
		int lda = mesh->NCDof, ldb = 1, ldc = 1;

		for (int jj = 0; jj < nip; jj++)													//
		{
			for (int ii = 0; ii < nip; ii++)
			{
				integratePoint[0] = GaussPoints[ii];										//
				integratePoint[1] = GaussPoints[jj];

				globalPoint = mesh->mapLocalToGlobalForFace(integratePoint, FVC);			// 把2维高斯点转换为全局坐标；
				LoadFunction(globalPoint, loadValue);										// 计算该点的载荷值

				localPoint = hex->mapGlobalToLocal(globalPoint);							//将全局坐标转换为单元内的局部坐标
				ShapeN = hex->calcShapeN(localPoint);										//根据局部坐标获得形函数矩阵的值
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







