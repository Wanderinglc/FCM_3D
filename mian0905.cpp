#include <iostream>
#include <vector>
#include <algorithm>

#include "Mesh3D.h"
#include "Edge.h"
#include "Node.h"
#include "EmbeddedDomain.h"
#include "Material.h"
#include "GaussIntegrator.h"
#include "FiniteCellManager.h"
#include "NeumannBoundaryCondition.h"
#include "DirichletBoundaryCondition.h"
#include "QuasiStaticAnalysis.h"


int main()
{
// 初始参数设置
	double MeshOrigin[3] = { 0.0, 0.0, 0.0 };
	double MeshLengths[3] = { 8.0, 8.0, 2.0 };
	int Division[3] = { 2, 2, 2 };	
	int DofDimension = 3;
	int PolynomialDegree = 3;

	Mesh3D mesh(MeshOrigin, MeshLengths, Division, PolynomialDegree, DofDimension);	


// 测试联系数组是够正确
	//mesh.elementVector_[0].getLocationMatrix();

// 测试 嵌入域 函数
	
	double center[3] = { 8.0,0.0 };
	double holeRadius = 1;
	double penaltyAlpha = 1E-10;
	ThickPlateWithAHole plate(MeshOrigin, MeshLengths, center, holeRadius, penaltyAlpha);

	//double coords1[3] = { 7.0,0.0,1.0 };	int index1 = plate.getDomainIndex(coords1);
	//printf("\nThe index of point 1 is: %d  (1: Fictitious, 2: Physics)\n", index1);
	

// 测试 材料参数 矩阵	
	double YoungModulus = 200;
	double PoissonsRatio = 0.2;
	double Density = 1;
	Hooke3D mat3D(YoungModulus, PoissonsRatio, Density);
	//mat3D.printMaterialMatrix();	


// 积分器	
	GaussIntegrator integrator(PolynomialDegree);	
	

// 测试 Element 的成员函数
	
	//Hexahedral* ele = &mesh.elementVector_[1];
	//// 顶点坐标输出
	//ele->printCellVertexCoords();
	//double* vc = ele->getCellVertexCoords();

	//// 判断单元是否被分割
	//int state = ele->isIntersected(plate, vc);
	//printf("\nElement position: %d \n", state);

	//// 将局部坐标映射到全局坐标
	//double Local[3] = { 0.0,0.0,0.0 };
	//double* point = ele->mapLocalToGlobal(Local, vc);
	//printf("\nPoint Coordinate in Element center: [%-6.2f, %-6.2f, %-6.2f] \n", point[0], point[1], point[2]);

	//// Jacobi矩阵测试
	//double* Jac = ele->calcJacobiAtCenter();

	//printf("\nJacobi Matrix J: ");
	//for (int i = 0; i < 9; i++)
	//{
	//	if (i % 3 == 0)  printf("| \n |  ");
	//	printf("%-7.2f  ", Jac[i]);
	//}
	//printf("| \n");

	//// 测试形函数
	//double coords2[3] = { 0.2,0.2,0.2 };
	////ele->printShape_dN(coords2);
	////ele->printShapeN(coords2);
	////ele->printMatrixB(coords2);


	//int spaceTreeDepth = 2;

	//double* Kc = ele->calcStiffnessMatrix(plate, mat3D, integrator, spaceTreeDepth);
	//// ele->printStiffnessMatrix(Kc);

// Finite Cell Manager
	FiniteCellManager fcm;
	fcm.mesh = &mesh;
	fcm.domain = &plate;
	fcm.material = &mat3D;
	fcm.integrator = &integrator;

	fcm.spaceTreeDepth = 3;


	fcm.generateMesh();

	//double* K = fcm.calcKsparse();


	//double point11[3] = { 4.0,4.0,0.0 };
	//int isNode = fcm.mesh->checkPointIsNode(point11);
	//printf("\nPoint 11 is Node: %d  (0: not a node;  1: is a node)\n", isNode);

	//double maxValue = *std::max_element(point11, point11 + 3);
	//printf("\nmax is %f\n", maxValue);


	//int* index = fcm.mesh->convertCoordinatesIntoIndices(point11);
	//printf("\nIndex are : [%d, %d, %d]\n", index[0], index[1], index[2]);


	//double faceStart[3] = { 0.0, 0.0, 2.0 };
	//double faceEnd[3] = { 8.0, 8.0, 2.0 };
	//std::vector<int> facesId;

	//fcm.mesh->findFaces(facesId, faceStart, faceEnd);
	//for (size_t i = 0; i < facesId.size(); i++)
	//{
	//	printf("face to constrain: %d\n", facesId[i]);
	//}

	//------------------------------------------------------------------------------ ==============计算载荷向量==============
	double tractionValue[3] = { 0.0, 100.0, 0.0 };
	LoadFunctionObject* pLoadFunc = new  LoadFunctionObject(tractionValue);

	double loadFaceStart[3] = { 0.0, 8.0, 0.0 };
	double loadFaceEnd[3] = { 8.0, 8.0, 2.0 };
	GaussIntegrator* pIntegrator = &integrator;

	FaceNeumannBoundaryCondition NeumannBC1(loadFaceStart, loadFaceEnd, pLoadFunc, pIntegrator);
	//double* Fv = new double[mesh.totalDofs]();
	//NeumannBC.calcLoadVector(&mesh, Fv);
	//printf("\nThe load vector is:\n");
	//for (int i = 0; i < mesh.totalDofs; i++)
	//{
	//	printf("F[%d] = %f\n", i + 1, Fv[i]);
	//}

	double prescribedValue[3] = { 0.0,0.0,0.0 };
	double PenaltyValue = 10E10;
	// Symmetry Boundary Conditions
	// XY
	double constraintFaceStart1[3] = { 0.0, 0.0, 0.0 };
	double constraintFaceEnd1[3] =   { 8.0, 8.0, 0.0 };
	int fixedDirection1[3] = { 0, 0, 1 };
	FaceDirichletBoundaryCondition DirichletBC_XY(constraintFaceStart1, constraintFaceEnd1, prescribedValue, fixedDirection1, PenaltyValue);
	// YZ
	double constraintFaceStart2[3] = { 8.0, 0.0, 0.0 };
	double constraintFaceEnd2[3] =   { 8.0, 8.0, 2.0 };
	int fixedDirection2[3] = { 1, 0, 0 };
	FaceDirichletBoundaryCondition DirichletBC_YZ(constraintFaceStart2, constraintFaceEnd2, prescribedValue, fixedDirection2, PenaltyValue);
	// XZ
	double constraintFaceStart3[3] = { 0.0, 0.0, 0.0 };
	double constraintFaceEnd3[3] =   { 8.0, 0.0, 2.0 };
	int fixedDirection3[3] = { 0, 1, 0 };
	FaceDirichletBoundaryCondition DirichletBC_XZ(constraintFaceStart3, constraintFaceEnd3, prescribedValue, fixedDirection3, PenaltyValue);

	QuasiStaticAnalysis myAnalysis;

	myAnalysis.fcManager = &fcm;
	myAnalysis.NeumannBC.push_back(&NeumannBC1);
	myAnalysis.DirichletBC.push_back(&DirichletBC_XY);
	myAnalysis.DirichletBC.push_back(&DirichletBC_YZ);
	myAnalysis.DirichletBC.push_back(&DirichletBC_XZ);

	
	myAnalysis.execute();







	getchar();
	return 0;
}


