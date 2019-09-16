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
	Edge() {}

	Edge(const int polynomialDegree, const int dofDimension):
		numberOfDofsPerFieldComponent(polynomialDegree - 1)
	{
		dofs.resize(numberOfDofsPerFieldComponent*dofDimension);
	}

	// Destructor
	virtual ~Edge()
	{
		//delete[] nodes;
	}

	// ��Ա����
public:
	int id;									// �ߵı��
	int numberOfDofsPerFieldComponent;		// �����������ɶ���
	Node* nodes[2];
	std::vector<Dof> dofs;
};



#endif
