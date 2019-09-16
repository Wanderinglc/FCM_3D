#ifndef Face_H
#define Face_H

/*=======================================================================
						 Face
	Define the attributes of Face for finite element method.
 =======================================================================*/
#include <vector>
#include "Dof.h"
#include "Node.h"
#include "Edge.h"


class Face
{
	// 成员函数
public:
	// constructor
	Face(){}
	Face(const int polynomialDegree, const int dofDimension)
	{
		numberOfDofsPerFieldComponent = (polynomialDegree - 1)*(polynomialDegree - 1);
		dofs.resize(numberOfDofsPerFieldComponent*dofDimension);
	}

	// Destructor
	virtual ~Face()
	{
		//delete[] nodes;
		//delete[] edges;
	}

	// 成员变量
public:
	int id;									// 面的编号
	int numberOfDofsPerFieldComponent;		// 场变量的自由度数					
	
	Node* nodes[4];							// 每个面有四个结点
	Edge* edges[4];							// 每个面有四条边
	std::vector<Dof> dofs;					// 面上的自由度	

};




#endif
