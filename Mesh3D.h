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
	Mesh3D();

	// Destructor
	virtual ~Mesh3D();

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





	// Data member
public:
	double MeshOrigin[3];					// the start point of mesh
	double lateralSize[3];					// length of mesh domain in three dimensions
	int numCells[3];						// number of cells in three dimensions

	
	int totalNodes_;
	int totalEdges_;
	int totalFaces_;
	int totalBulks_;
	int totalCells_;

	int polynomialDegree_;
	int dofDimension_;

	std::vector<Node> nodeVector_;          // global nodes list
	std::vector<Edge> edgeVector_;
	std::vector<Face> faceVector_;
	std::vector<Bulk> bulkVector_;
	std::vector<Hexahedral> elementVector_;   // global element list

};


#endif
