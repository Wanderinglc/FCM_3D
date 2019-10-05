#ifndef Edge_H
#define Edge_H

/*=======================================================================
						 Edge
	Define the attributes of Edge for finite element method.
 =======================================================================*/
#include <vector>
#include "Dof.h"
#include "Node.h"

class Edge
{
	// ��Ա����
public:
	// constructor
	Edge() :available(true), shareCells(0) {}

	Edge(const int polynomialDegree, const int dofDimension):
		numberOfDofsPerFieldComponent(polynomialDegree - 1)
	{
		dofs.resize(numberOfDofsPerFieldComponent*dofDimension);
	}

	// Destructor
	~Edge()
	{
		//delete[] nodes;
	}

	// ��Ա����
public:
	int id;									// �ߵı��
	int numberOfDofsPerFieldComponent;		// �����������ɶ���
	int nodesId[2];
	//Node* nodes[2];
	std::vector<Dof> dofs;

	bool available;						// �����жϱ��Ƿ����
	int shareCells;						// �߹���Cell�ĸ���
	std::vector<int> shareCellsId;		// �߹����Cell���
};



#endif
