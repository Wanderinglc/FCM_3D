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
	// ��Ա����
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

	// ��Ա����
public:
	int id;									// ʵ��ı��
	int numberOfDofsPerFieldComponent;		// �����������ɶ���

	int nodesId[8];
	int edgesId[12];
	int facesId[6];

	//Node* nodes[8];						// ÿ������8�����
	//Edge* edges[12];						// ÿ������12����
	//Face* faces[6];						// ÿ������6����
	std::vector<Dof> dofs;					// ʵ���ϵ����ɶ�

	
	bool available;							// �����ж����Ƿ����
	static const int shareCells = 1;		// �干��Cell�ĸ���
	int shareCellsId;						// �湲��İ�Ԫ���

};



#endif
