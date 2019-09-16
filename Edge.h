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
	// 成员函数
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

	// 成员变量
public:
	int id;									// 边的编号
	int numberOfDofsPerFieldComponent;		// 场变量的自由度数
	Node* nodes[2];
	std::vector<Dof> dofs;
};



#endif
