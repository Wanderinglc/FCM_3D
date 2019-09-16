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
	// ��Ա����
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

	// ��Ա����
public:
	int id;									// ��ı��
	int numberOfDofsPerFieldComponent;		// �����������ɶ���					
	
	Node* nodes[4];							// ÿ�������ĸ����
	Edge* edges[4];							// ÿ������������
	std::vector<Dof> dofs;					// ���ϵ����ɶ�	

};




#endif
