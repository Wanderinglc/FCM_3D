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

	// ��Ա����
public:
	int id;								// �����
	double coords_[3];					// ��������
	std::vector<Dof> dofs;				// ����ϵ����ɶ�

	bool available;						// �����жϽ���Ƿ����
	int shareCells;						// ��㹲��Cell�ĸ���
	std::vector<int> shareCellsId;		// ��㹲���Cell���

};


#endif
