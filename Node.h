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
	Node()	:available(true), shareCells(0)
	{
	}

	Node(const int dofDimension)
	{
		dofs.resize(dofDimension);
	}

	// Destructor
	~Node()
	{
		// nothing for now 
	}

	// 成员变量
public:
	int id;								// 结点编号
	double coords_[3];					// 结点的坐标
	std::vector<Dof> dofs;				// 结点上的自由度

	bool available;						// 用来判断结点是否存在
	int shareCells;						// 结点共享Cell的个数
	std::vector<int> shareCellsId;		// 结点共享的Cell编号

};


#endif
