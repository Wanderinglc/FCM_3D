#ifndef Bulk_H
#define Bulk_H

/*=======================================================================
						 Bulk
	Define the attributes of Bulk for finite element method.
 =======================================================================*/
#include <vector>
#include "Dof.h"
#include "Node.h"
#include "Edge.h"
#include "Face.h"

class Bulk
{
	// 成员函数
public:
	// constructor
	Bulk() :available(true) {}
	Bulk(const int polynomialDegree, const int dofDimension)
	{
		numberOfDofsPerFieldComponent = (polynomialDegree - 1)
			*(polynomialDegree - 1)*(polynomialDegree - 1);
		dofs.resize(numberOfDofsPerFieldComponent*dofDimension);
	}

	// Destructor
	~Bulk()
	{
		//delete[] nodes;
		//delete[] edges;
		//delete[] faces;
	}

	// 成员变量
public:
	int id;									// 实体的编号
	int numberOfDofsPerFieldComponent;		// 场变量的自由度数

	int nodesId[8];
	int edgesId[12];
	int facesId[6];

	//Node* nodes[8];						// 每个体有8个结点
	//Edge* edges[12];						// 每个体有12条边
	//Face* faces[6];						// 每个体有6条边
	std::vector<Dof> dofs;					// 实体上的自由度

	
	bool available;							// 用来判断体是否存在
	static const int shareCells = 1;		// 体共享Cell的个数
	int shareCellsId;						// 面共享的胞元编号

};



#endif
