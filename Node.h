#ifndef NODE_H
#define NODE_H

/*=======================================================================
						 Node
	Define the attributes of node for finite element method.
 =======================================================================*/
#include <vector>
#include "Dof.h"

class Node
{
	// 成员函数
public:
	// constructor
	Node()	{}

	Node(const int dofDimension)
	{
		dofs.resize(dofDimension);
	}

	// Destructor
	virtual ~Node()
	{
		// nothing for now 
	}

	// 成员变量
public:
	int id;
	double coords_[3];
	std::vector<Dof> dofs;
};


#endif
