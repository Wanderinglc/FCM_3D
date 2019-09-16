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
	// ��Ա����
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

	// ��Ա����
public:
	int id;
	double coords_[3];
	std::vector<Dof> dofs;
};


#endif
