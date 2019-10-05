#ifndef Mesh3D_H
#define Mesh3D_H

#include <vector>
#include "Dof.h"
#include "Node.h"
#include "Edge.h"
#include "Face.h"
#include "Bulk.h"
#include "Element.h"

class Mesh3D
{
public:
	// Constructor
	Mesh3D() {}

	Mesh3D(double* start, double* length, int* divide, int PD, int dofDim)
		:polynomialDegree_(PD), dofDimension_(dofDim)
	{
		NCNode = (PD + 1)*(PD + 1)*(PD + 1);
		NCDof = NCNode * dofDim;

		meshOrigin[0] = start[0];
		meshOrigin[1] = start[1];
		meshOrigin[2] = start[2];

		lateralSize[0] = length[0];
		lateralSize[1] = length[1];
		lateralSize[2] = length[2];

		numCells[0] = divide[0];
		numCells[1] = divide[1];
		numCells[2] = divide[2];
	}

	// Destructor
	~Mesh3D();

	//void discretize_extended_domian(double original_point, int NCX, int NCY);
	// ���ɽ��
	void creatNode();						// ���ɽ��
	void printNode() const;
	// ���ɱ�
	void creatEdge();
	void printEdge() const;
	void creatEdgeOnDirection1();
	void creatEdgeOnDirection2();
	void creatEdgeOnDirection3();

	// ������
	void creatFace();
	void creatFaceXY();
	void creatFaceYZ();
	void creatFaceXZ();

	// ������
	void creatBulk();
	void printBulk();

	// ���ɵ�Ԫ
	void creatElement();

	// ���ɶȱ��
	void assignDofs();
	int assignNodesToVertices(int &TotalDofCounter);
	int assignNodesToEdges(int CurrentPD, int &TotalDofCounter);
	int assignNodesToFaces(int CurrentPD, int &TotalDofCounter);
	int assignNodesToBulks(int CurrentPD, int &TotalDofCounter);

	void generateStructuredMesh();

	void findFaces(std::vector<int> &faces, double* faceStart, double* faceEnd);
	bool checkPointIsNode(double* Coordinate);
	int *convertCoordinatesIntoIndices(double *Coordinate);
	int findElementByPoint(double *Coordinate);

	double* getFaceVertexCoords(int faceId);
	double* getEdgeVertexCoords(int edgeId);
	double* mapLocalToGlobalForFace(double *Local, double *VC);

	double calcFaceDetJ(int faceId);
	double calcEdgeDetJ(int edgeId);







	// Data member
public:
	static const int Dimention = 3;
	double meshOrigin[3];					// the start point of mesh
	double lateralSize[3];					// length of mesh domain in three dimensions
	int numCells[3];						// number of cells in three dimensions

	
	int totalNodes_;
	int totalEdges_;
	int totalFaces_;
	int totalBulks_;
	int totalCells_;
	int totalDofs;
	int NCDof;								// ��Ԫ�ڵ������ɶ���
	int NCNode;								// ��Ԫ�ڵ��ܽڵ���

	int polynomialDegree_;
	int dofDimension_;

	std::vector<Node> nodeVector_;          // global nodes list
	std::vector<Edge> edgeVector_;
	std::vector<Face> faceVector_;
	std::vector<Bulk> bulkVector_;
	std::vector<Hexahedral> elementVector_;   // global element list

};


#endif
