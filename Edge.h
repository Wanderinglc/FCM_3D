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

	// 成员变量
public:
	int id;									// 边的编号
	int numberOfDofsPerFieldComponent;		// 场变量的自由度数
	int nodesId[2];
	//Node* nodes[2];
	std::vector<Dof> dofs;

	bool available;						// 用来判断边是否存在
	int shareCells;						// 边共享Cell的个数
	std::vector<int> shareCellsId;		// 边共享的Cell编号
};



#endif
