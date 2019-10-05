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

int main()
{
// ��ʼ��������
	double MeshOrigin[3] = { 0.0, 0.0, 0.0 };
	double MeshLengths[3] = { 8.0, 8.0, 2.0 };
	int Division[3] = { 2, 2, 1 };
	int PolynomialDegree = 3;
	int DofDimension = 3;

	Mesh3D mesh(MeshOrigin, MeshLengths, Division, PolynomialDegree, DofDimension);	


// ������ϵ�����ǹ���ȷ
	//mesh.elementVector_[0].getLocationMatrix();

// ���� Ƕ���� ����
	
	double center[3] = { 8.0,0.0 };
	double holeRadius = 1;
	double penaltyAlpha = 1E-10;
	ThickPlateWithAHole plate(MeshOrigin, MeshLengths, center, holeRadius, penaltyAlpha);

	//double coords1[3] = { 7.0,0.0,1.0 };	int index1 = plate.getDomainIndex(coords1);
	//printf("\nThe index of point 1 is: %d  (1: Fictitious, 2: Physics)\n", index1);
	

// ���� ���ϲ��� ����	
	double YoungModulus = 2000;
	double PoissonsRatio = 0.2;
	double Density = 1;
	Hooke3D mat3D(YoungModulus, PoissonsRatio, Density);
	//mat3D.printMaterialMatrix();	


// ������	
	GaussIntegrator integrator(PolynomialDegree);	
	

// ���� Element �ĳ�Ա����
	
	//Hexahedral* ele = &mesh.elementVector_[1];
	//// �����������
	//ele->printCellVertexCoords();
	//double* vc = ele->getCellVertexCoords();

	//// �жϵ�Ԫ�Ƿ񱻷ָ�
	//int state = ele->isIntersected(plate, vc);
	//printf("\nElement position: %d \n", state);

	//// ���ֲ�����ӳ�䵽ȫ������
	//double Local[3] = { 0.0,0.0,0.0 };
	//double* point = ele->mapLocalToGlobal(Local, vc);
	//printf("\nPoint Coordinate in Element center: [%-6.2f, %-6.2f, %-6.2f] \n", point[0], point[1], point[2]);

	//// Jacobi�������
	//double* Jac = ele->calcJacobiAtCenter();

	//printf("\nJacobi Matrix J: ");
	//for (int i = 0; i < 9; i++)
	//{
	//	if (i % 3 == 0)  printf("| \n |  ");
	//	printf("%-7.2f  ", Jac[i]);
	//}
	//printf("| \n");

	//// �����κ���
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


	double point11[3] = { 4.0,4.0,0.0 };
	int isNode = fcm.mesh->checkPointIsNode(point11);
	printf("\nPoint 11 is Node: %d  (0: not a node;  1: is a node)\n", isNode);

	double maxValue = *std::max_element(point11, point11 + 3);
	printf("\nmax is %f\n", maxValue);


	int* index = fcm.mesh->convertCoordinatesIntoIndices(point11);
	printf("\nIndex are : [%d, %d, %d]\n", index[0], index[1], index[2]);


	double faceStart[3] = { 0.0, 0.0, 2.0 };
	double faceEnd[3] = { 8.0, 8.0, 2.0 };
	std::vector<int> facesId;

	fcm.mesh->findFaces(facesId, faceStart, faceEnd);
	for (size_t i = 0; i < facesId.size(); i++)
	{
		printf("face to constrain: %d\n", facesId[i]);
	}





	getchar();
	return 0;
}


