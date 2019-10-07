
#include <iostream>
#include <vector>
#include <algorithm>
#include <numeric>

#include "Mesh3D.h"
#include "Dof.h"
#include "Node.h"
#include "Edge.h"
#include "Face.h"
#include "Bulk.h"
#include "Element.h"



//Mesh3D::Mesh3D()
//{
//
//}

Mesh3D::~Mesh3D()
{

}

/*---------------------------------------
			���ɽ��
  ---------------------------------------*/
void Mesh3D::creatNode()
{
	//int dCell[3];
	double xIncrement = (lateralSize[0]) / numCells[0];		// x ��������
	double yIncrement = (lateralSize[1]) / numCells[1];		// y ��������
	double zIncrement = (lateralSize[2]) / numCells[2];		// y ��������

	int numNodes[3];											// x,y,z ���������ϵĽ����
	for (int i = 0; i < 3; i++)
		numNodes[i] = numCells[i] + 1;

	totalNodes_ = numNodes[0] * numNodes[1] * numNodes[2];		//�ܽ����

	// Initialize nodeVector_
	for (int k = 0; k < numNodes[2]; k++)
	{
		for (int j = 0; j < numNodes[1]; j++)
		{
			for (int i = 0; i < numNodes[0]; i++)
			{
				Node node(dofDimension_);
				node.coords_[0] = meshOrigin[0] + i * xIncrement;
				node.coords_[1] = meshOrigin[1] + j * yIncrement;
				node.coords_[2] = meshOrigin[2] + k * zIncrement;
				node.id = i + j * numNodes[0] + k * numNodes[0]* numNodes[1];
				nodeVector_.push_back(node);
			}
		}
	}

}
//		��������Ϣ
void Mesh3D::printNode() const
{
	std::cout << "----------Print Nodes----------\n";
	int numNodes[3];											// x,y,z ���������ϵĽ����
	for (int i = 0; i < 3; i++)		
		numNodes[i] = numCells[i] + 1;

	int index = 0;
	for (int k = 0; k < numNodes[2]; k++)
	{
		std::cout << "Node num in " << k << " floor\n";
		for (int j = 0; j < numNodes[1]; j++)
		{
			for (int i = 0; i < numNodes[0]; i++)
			{
				index = i + j * numNodes[0] + k * numNodes[0] * numNodes[1];
				Node node = nodeVector_[index];
				printf("%-4d", node.id);
			}
			printf("\n");
		}
	}

	printf("\n Node Coordinate\n");
	for (int k = 0; k < totalNodes_; k++)
	{
		Node node = nodeVector_[k];
		printf("Node %4d : ( %5.1f, %5.1f, %5.1f )\n", node.id, node.coords_[0], node.coords_[1], node.coords_[2]);
	}

	printf("\n Node share cells:\n");
	for (int k = 0; k < totalNodes_; k++)
	{
		Node node = nodeVector_[k];
		printf("Node %4d :   share %d cells:  ", node.id, node.shareCells);
		for (int m = 0; m < node.shareCells; m++)
		{
			printf("%4d", node.shareCellsId[m]);
		}
		printf("\n");
	}
}



/*---------------------------------------
			���� ��
  ---------------------------------------*/
void Mesh3D::creatEdge()
{
	totalEdges_ = numCells[0] * (numCells[1] + 1)*(numCells[2] + 1)
		+ numCells[1] * (numCells[0] + 1)*(numCells[2] + 1)
		+ numCells[2] * (numCells[0] + 1)*(numCells[1] + 1);
	creatEdgeOnDirection1();
	creatEdgeOnDirection2();
	creatEdgeOnDirection3();
}
void Mesh3D::creatEdgeOnDirection1()
{
	int numNodesInXY = (numCells[0] + 1)*(numCells[1] + 1);		//xy ƽ���ڵĽ����
	int count = 0;

	int leftNode, rightNode;
	for (int k = 0; k < numCells[2] + 1; k++)
	{
		for (int j = 0; j < numCells[1] + 1; j++)
		{
			for (int i = 0; i < numCells[0]; i++)
			{
				Edge edg(polynomialDegree_, dofDimension_);
				// �ߵ��������������
				edg.id = count; 
				count++;
				leftNode = k * (numNodesInXY)+j * (numCells[0] + 1) + i;
				rightNode = leftNode + 1;

				edg.nodesId[0] = leftNode;
				edg.nodesId[1] = rightNode;
				edgeVector_.push_back(edg);
			}
		}
	}
}
void Mesh3D::creatEdgeOnDirection2()
{
	int numNodesInXY = (numCells[0] + 1)*(numCells[1] + 1);
	int count = numCells[0] * (numCells[1] + 1)*(numCells[2] + 1);

	int leftNode, rightNode;
	for (int k = 0; k < numCells[2] + 1; k++)
	{
		for (int i = 0; i < numCells[0] + 1; i++)
		{
			for (int j = 0; j < numCells[1]; j++)
			{
				Edge edg(polynomialDegree_, dofDimension_);
				// �ߵ��������������
				edg.id = count;
				count++;
				leftNode = k * (numNodesInXY)+j * (numCells[0] + 1) + i;
				rightNode = leftNode + numCells[0] + 1;

				edg.nodesId[0] = leftNode;
				edg.nodesId[1] = rightNode;
				edgeVector_.push_back(edg);
			}
		}
	}

}
void Mesh3D::creatEdgeOnDirection3()
{
	int numNodesInXY = (numCells[0] + 1)*(numCells[1] + 1);
	int count = numCells[0] * (numCells[1] + 1)*(numCells[2] + 1)
		+ numCells[1] * (numCells[0] + 1)*(numCells[2] + 1);

	int leftNode, rightNode;
	for (int j = 0; j < numCells[1] + 1; j++)
	{
		for (int i = 0; i < numCells[0] + 1; i++)
		{
			for (int k = 0; k < numCells[2]; k++)
			{
				Edge edg(polynomialDegree_, dofDimension_);
				// �ߵ��������������
				edg.id = count;
				count++;
				leftNode = k * (numNodesInXY)+j * (numCells[0] + 1) + i;
				rightNode = leftNode + numNodesInXY;

				edg.nodesId[0] = leftNode;
				edg.nodesId[1] = rightNode;
				edgeVector_.push_back(edg);
			}
		}
	}

}
//        ��� �� ��Ϣ
void Mesh3D::printEdge() const
{
	std::cout << "----------Print Edge----------\n";
	int totalEdges = numCells[0] * (numCells[1] + 1)*(numCells[2] + 1)
		+ numCells[1] * (numCells[0] + 1)*(numCells[2] + 1)
		+ numCells[2] * (numCells[0] + 1)*(numCells[1] + 1);
	//const Node *n1, *n2;
	printf("\n LeftNode  RightNode\n");
	for (int k = 0; k < totalEdges; k++)
	{
		const Edge* edg = &edgeVector_[k];		
		printf("%4d\t%4d\n",edg->nodesId[0], edg->nodesId[1]);
	}
	printf("\n");


	printf("\n Edge share cells:\n");
	for (int k = 0; k < totalEdges; k++)
	{
		const Edge* edg = &edgeVector_[k];
		printf("Edge %4d :   share %d cells:  ", edg->id, edg->shareCells);
		for (int m = 0; m < edg->shareCells; m++)
		{
			printf("%4d", edg->shareCellsId[m]);
		}
		printf("\n");
	}
	printf("\n");

}

/*---------------------------------------
			������
  ---------------------------------------*/
void Mesh3D::creatFace()
{
	totalFaces_ = numCells[0] * numCells[1] * (numCells[2] + 1)
		+ numCells[1] * numCells[2] * (numCells[0] + 1)
		+ numCells[0] * numCells[2] *(numCells[1] + 1);
	creatFaceXY();
	creatFaceYZ();
	creatFaceXZ();
}
void Mesh3D::creatFaceXY()
{
	int numOfXEdges = numCells[0] * (numCells[1] + 1)*(numCells[2] + 1);
	
	int numXEdgesInXY = numCells[0] * (numCells[1] + 1);
	int numYEdgesInXY = numCells[1] * (numCells[0] + 1);
	//printf("\n");
	int count = 0;
	// loop in z direction
	for (int k = 0; k < numCells[2] + 1; k++)
	{
		for (int j = 0; j < numCells[1]; j++)
		{
			for (int i = 0; i < numCells[0]; i++)
			{
				Face fa(polynomialDegree_, dofDimension_);
				fa.id = count;	count++;

				int edge1 = i + j * numCells[0] + k * numXEdgesInXY;
				int edge4 = j + i * numCells[1] + k * numYEdgesInXY + numOfXEdges;
				int edge2 = edge4 + numCells[1];
				int edge3 = edge1 + numCells[0];
				//printf("Face %d :\n", fa.id);
				//printf("%d  %d  %d  %d\n", edge1, edge2, edge3, edge4);

				fa.edgesId[0] = edge1;
				fa.edgesId[1] = edge2;
				fa.edgesId[2] = edge3;
				fa.edgesId[3] = edge4;

				fa.nodesId[0] = edgeVector_[edge1].nodesId[0];
				fa.nodesId[1] = edgeVector_[edge2].nodesId[0];
				fa.nodesId[2] = edgeVector_[edge3].nodesId[1];
				fa.nodesId[3] = edgeVector_[edge4].nodesId[1];

				//for (int n = 0; n < 4; n++)
				//{
				//	printf("-----node id %d: %4.1f, %4.1f, %4.1f\n", fa.nodes[n]->id, fa.nodes[n]->coords_[0], fa.nodes[n]->coords_[1],
				//		fa.nodes[n]->coords_[2]);
				//}
				faceVector_.push_back(fa);

			}
		}
	}
}
void Mesh3D::creatFaceYZ()
{
	int numOfXEdges = numCells[0] * (numCells[1] + 1)*(numCells[2] + 1);
	int numOfYEdges = numCells[1] * (numCells[0] + 1)*(numCells[2] + 1);

	int numYEdgesInXY = numCells[1] * (numCells[0] + 1);
	int numZEdgesInXZ = numCells[2] * (numCells[0] + 1);
	//printf("\n");
	int count = numCells[0] * numCells[1] * (numCells[2] + 1);

	// loop in x direction
	for (int i = 0; i < numCells[0] + 1; i++)
	{
		for (int k = 0; k < numCells[2]; k++)
		{
			for (int j = 0; j < numCells[1]; j++)
			{
				Face fa(polynomialDegree_, dofDimension_);
				fa.id = count;	count++;
				
				int edge1 = j + i * numCells[1] + k * numYEdgesInXY + numOfXEdges;
				int edge4 = k + j * numZEdgesInXZ + i * numCells[2] + numOfXEdges + numOfYEdges;
				int edge3 = edge1 + numYEdgesInXY;
				int edge2 = edge4 + numZEdgesInXZ;

				//printf("Face %d :  ", fa.id);
				//printf("%d  %d  %d  %d\n", edge1, edge2, edge3, edge4);
				fa.edgesId[0] = edge1;
				fa.edgesId[1] = edge2;
				fa.edgesId[2] = edge3;
				fa.edgesId[3] = edge4;

				fa.nodesId[0] = edgeVector_[edge1].nodesId[0];
				fa.nodesId[1] = edgeVector_[edge2].nodesId[0];
				fa.nodesId[2] = edgeVector_[edge3].nodesId[1];
				fa.nodesId[3] = edgeVector_[edge4].nodesId[1];

				//for (int n = 0; n < 4; n++)
				//	printf("-----node id %d: %4.1f, %4.1f, %4.1f\n", fa.nodes[n]->id, fa.nodes[n]->coords_[0], fa.nodes[n]->coords_[1],
				//		fa.nodes[n]->coords_[2]);;
				faceVector_.push_back(fa);

			}
		}
	}
}
void Mesh3D::creatFaceXZ()
{
	int numOfXEdges = numCells[0] * (numCells[1] + 1)*(numCells[2] + 1);
	int numOfYEdges = numCells[1] * (numCells[0] + 1)*(numCells[2] + 1);

	int numXEdgesInXY = numCells[0] * (numCells[1] + 1);
	int numZEdgesInXZ = numCells[2] * (numCells[0] + 1);
	int numZEdgesInYZ = numCells[2] * (numCells[1] + 1);
	//printf("\n");
	int count = numCells[0] * numCells[1] * (numCells[2] + 1)
		+ numCells[1] * numCells[2] * (numCells[0] + 1);

	// loop in y direction
	for (int j = 0; j < numCells[1] + 1; j++)
	{
		for (int k = 0; k < numCells[2]; k++)
		{
			for (int i = 0; i < numCells[0]; i++)
			{
				Face fa(polynomialDegree_, dofDimension_);
				fa.id = count;	count++;

				int edge1 = k + j * numZEdgesInXZ + i * numCells[2] + numOfXEdges + numOfYEdges;
				int edge4 = i + j * numCells[0] + k * numXEdgesInXY;
				int edge2 = edge4 + numXEdgesInXY;
				int edge3 = edge1 + numCells[2];

				//printf("Face %d :  ", fa.id);
				//printf("%d  %d  %d  %d\n", edge1, edge2, edge3, edge4);
				fa.edgesId[0] = edge1;
				fa.edgesId[1] = edge2;
				fa.edgesId[2] = edge3;
				fa.edgesId[3] = edge4;

				fa.nodesId[0] = edgeVector_[edge1].nodesId[0];
				fa.nodesId[1] = edgeVector_[edge2].nodesId[0];
				fa.nodesId[2] = edgeVector_[edge3].nodesId[1];
				fa.nodesId[3] = edgeVector_[edge4].nodesId[1];

				//for (int n = 0; n < 4; n++)
				//	printf("-----node id %d: %4.1f, %4.1f, %4.1f\n", fa.nodes[n]->id, fa.nodes[n]->coords_[0], fa.nodes[n]->coords_[1],
				//		fa.nodes[n]->coords_[2]);;
				faceVector_.push_back(fa);

			}
		}
	}
}

/*---------------------------------------
			������
  ---------------------------------------*/
void Mesh3D::creatBulk()
{
	totalBulks_ = numCells[0] * numCells[1] * numCells[2];

	int numNodesInXY = (numCells[0] + 1)*(numCells[1] + 1);

	int NumOfXEdges = numCells[0] * (numCells[1] + 1)*(numCells[2] + 1);
	int NumOfYEdges = numCells[1] * (numCells[0] + 1)*(numCells[2] + 1);
	
	int NumOfXYFaces = numCells[0] * numCells[1] * (numCells[2] + 1);
	int NumOfYZFaces = numCells[1] * numCells[2] * (numCells[0] + 1);
	
	int count = 0;
	for (int k = 0; k < numCells[2]; k++)
	{
		for (int j = 0; j < numCells[1]; j++)
		{
			for (int i = 0; i < numCells[0]; i++)
			{
				Bulk cube(polynomialDegree_, dofDimension_);
				cube.id = count;	count++;

				// ������ ��� ���
				cube.nodesId[0] = i + j * (numCells[0] + 1) + k * numNodesInXY;
				cube.nodesId[1] = cube.nodesId[0] + 1;
				cube.nodesId[3] = cube.nodesId[0] + numCells[0] + 1;
				cube.nodesId[2] = cube.nodesId[3] + 1;
				cube.nodesId[4] = cube.nodesId[0] + numNodesInXY;
				cube.nodesId[5] = cube.nodesId[4] + 1;
				cube.nodesId[7] = cube.nodesId[4] + numCells[0] + 1;
				cube.nodesId[6] = cube.nodesId[7] + 1;

				// ������ �� ���
				cube.edgesId[0] = i + j * numCells[0] + k * numCells[0] * (numCells[1] + 1);
				cube.edgesId[2] = cube.edgesId[0] + numCells[0];
				cube.edgesId[3] = NumOfXEdges + j + i * numCells[1] + k * numCells[1] * (numCells[0] + 1);
				cube.edgesId[1] = cube.edgesId[3] + numCells[1];

				cube.edgesId[8] = cube.edgesId[0] + numCells[0] * (numCells[1] + 1);
				cube.edgesId[10] = cube.edgesId[8] + numCells[0];
				cube.edgesId[11] = cube.edgesId[3] + numCells[1] * (numCells[0] + 1);
				cube.edgesId[9] = cube.edgesId[11] + numCells[1];

				cube.edgesId[4] = NumOfXEdges + NumOfYEdges + k + i * numCells[2] + j * numCells[2] * (numCells[0] + 1);
				cube.edgesId[5] = cube.edgesId[4] + numCells[2];
				cube.edgesId[7] = cube.edgesId[4] + numCells[2] * (numCells[0] + 1);
				cube.edgesId[6] = cube.edgesId[7] + numCells[2];

				// ������ �� ���
				cube.facesId[0] = i + j * numCells[0] + k * numCells[0] * numCells[1];
				cube.facesId[5] = cube.facesId[0] + numCells[0] * numCells[1];
				cube.facesId[4] = NumOfXYFaces + j + k * numCells[1] + i * numCells[1] * numCells[2];
				cube.facesId[2] = cube.facesId[4] + numCells[1] * numCells[2];
				cube.facesId[1] = NumOfXYFaces + NumOfYZFaces + i + k * numCells[0] + j*numCells[0] * numCells[2];
				cube.facesId[3] = cube.facesId[1] + numCells[0] * numCells[2];
				
				//printf(" Face Num: ");
				//printf("\n");
				bulkVector_.push_back(cube);

			}
		}
	}

}
void Mesh3D::printBulk()
{
	std::cout << "----------Print Bulk----------\n";

	//for (int i = 0; i < bulkVector_.size(); i++) 
	//{
	//	const Bulk* cube = &bulkVector_[i];
	//	printf("\nCube %d: \n", cube->id);
	//	for (int n = 0; n < 8; n++)
	//	{
	//		printf("   node %d: %4.1f, %4.1f, %4.1f\n", cube->nodes[n]->id,
	//			cube->nodes[n]->coords_[0], cube->nodes[n]->coords_[1],
	//			cube->nodes[n]->coords_[2]);
	//	}

	//	printf("   \n");
	//	for (int e = 0; e < 12; e++)
	//	{
	//		printf("   middle point of edge %d: %4.1f, %4.1f, %4.1f\n", cube->edges[e]->id,
	//			(cube->edges[e]->nodes[0]->coords_[0]+ cube->edges[e]->nodes[1]->coords_[0])/2,
	//			(cube->edges[e]->nodes[0]->coords_[1] + cube->edges[e]->nodes[1]->coords_[1]) / 2, 
	//			(cube->edges[e]->nodes[0]->coords_[2] + cube->edges[e]->nodes[1]->coords_[2]) / 2 );
	//	}

	//	printf("   \n");
	//	for (int f = 0; f < 6; f++)
	//	{
	//		printf("   center point of face %d: %4.1f, %4.1f, %4.1f\n", cube->faces[f]->id,
	//			(cube->faces[f]->nodes[0]->coords_[0] + cube->faces[f]->nodes[2]->coords_[0]) / 2,
	//			(cube->faces[f]->nodes[0]->coords_[1] + cube->faces[f]->nodes[2]->coords_[1]) / 2,
	//			(cube->faces[f]->nodes[0]->coords_[2] + cube->faces[f]->nodes[2]->coords_[2]) / 2);
	//	}

	//}


}

/*---------------------------------------
			���ɵ�Ԫ
  ---------------------------------------*/
void Mesh3D::creatElement()
{
	totalCells_ = numCells[0] * numCells[1] * numCells[2];
	//elementVector_.resize(numCells[0] * numCells[1] * numCells[2]);

	for (int i = 0; i < bulkVector_.size(); i++)
	{
		Hexahedral hex(polynomialDegree_);
		hex.id = i;
		//printf("Element %d \n",hex.id);

		// ��Ԫ����
		hex.bulks[0] = &bulkVector_[i];										//�� �� ��ַ
		//hex.bulksId[0] = i;												//�� �� ���
		bulkVector_[i].shareCellsId = i;


		// ��Ԫ�Ľ��
		//printf("\n    Element Node: ");
		int index = 0;
		for (int iNode = 0; iNode < 8; iNode++)
		{
			index = bulkVector_[i].nodesId[iNode];
			hex.nodes[iNode] = &nodeVector_[index];							//�� ��� ��ַ
			//hex.nodesId[iNode] = bulkVector_[i].nodesId[iNode];			//�� ��� ���

			nodeVector_[index].shareCells++;								//��㹲���Cell���� +1
			nodeVector_[index].shareCellsId.push_back(i);					//��㹲���Cell���			
		}
			

		// ��Ԫ�ı�
		//printf("\n    Element Edge: ");
		for (int iEdge = 0; iEdge < 12; iEdge++)
		{
			index = bulkVector_[i].edgesId[iEdge];							//�߱��
			hex.edges[iEdge] = &edgeVector_[index];							//�� �� ��ַ
			//hex.edgesId[iEdge] = bulkVector_[i].edgesId[iEdge];			//�� �� ���

			edgeVector_[index].shareCells++;								//�߹����Cell���� +1
			edgeVector_[index].shareCellsId.push_back(i);					//�߹����Cell���			
		}
			

		// ��Ԫ����
		//printf("\n    Element face: ");
		for (int iFace = 0; iFace < 6; iFace++)
		{
			index = bulkVector_[i].facesId[iFace];
			hex.faces[iFace] = &faceVector_[index];							// �� �� ��ַ
			//hex.facesId[iFace] = bulkVector_[i].facesId[iFace];			// �� �� ���

			faceVector_[index].shareCells++;								//�湲���Cell���� +1
			faceVector_[index].shareCellsId.push_back(i);					//�湲���Cell���
		}
		//printf("\n");		
		elementVector_.push_back(hex);
	}
}


/*---------------------------------------
			���ɶȱ��
  ---------------------------------------*/
void Mesh3D::assignDofs()
{
	int totalDofCounter = 0;
	totalDofCounter = assignNodesToVertices(totalDofCounter);
	//printf("total Dof Counter %d\n", totalDofCounter);
	if (polynomialDegree_ > 1)
	{
		for (int p = 2; p <= polynomialDegree_; p++)
		{
			totalDofCounter = assignNodesToEdges(p, totalDofCounter);	// �����Ͻ����
			//printf("total Dof Counter %d\n", totalDofCounter);
			totalDofCounter = assignNodesToFaces(p, totalDofCounter);	// �����ڽ����
			//printf("total Dof Counter %d\n", totalDofCounter);
			totalDofCounter = assignNodesToBulks(p, totalDofCounter);
			//printf("total Dof Counter %d\n", totalDofCounter);
		}
	}
	totalDofs = totalDofCounter;

}

int Mesh3D::assignNodesToVertices(int &TotalDofCounter)
{
	for (int i = 0; i < nodeVector_.size(); i++)
	{
		for (int dim = 0; dim < dofDimension_; dim++)
		{
			nodeVector_[i].dofs[dim].id = TotalDofCounter;
			TotalDofCounter++;
		}
	}
	return TotalDofCounter;
}

int Mesh3D::assignNodesToEdges(int CurrentPD, int &TotalDofCounter)
{
	int m = (CurrentPD - 2)*dofDimension_;
	for (int i = 0; i < edgeVector_.size(); i++)
	{
		for (int dim = 0; dim < dofDimension_; dim++)
		{
			//printf("m+dim = %d\n", m + dim);
			edgeVector_[i].dofs[m + dim].id = TotalDofCounter;
			TotalDofCounter++;
		}
	}
	return TotalDofCounter;

}

int Mesh3D::assignNodesToFaces(int CurrentPD, int &TotalDofCounter)
{
	int start = (CurrentPD - 2)*(CurrentPD - 2);
	int end = (CurrentPD - 1)*(CurrentPD - 1);
	for (int i = 0; i < faceVector_.size(); i++)
	{
		for (int j = start; j < end; j++)
		{
			for (int dim = 0; dim < dofDimension_; dim++)
			{
				faceVector_[i].dofs[j*dofDimension_ + dim] = TotalDofCounter;
				TotalDofCounter++;
			}

		}
	}
	return TotalDofCounter;
}

int Mesh3D::assignNodesToBulks(int CurrentPD, int &TotalDofCounter)
{
	int start = (CurrentPD - 2)*(CurrentPD - 2)*(CurrentPD - 2);
	int end = (CurrentPD - 1)*(CurrentPD - 1)*(CurrentPD - 1);

	for (int i = 0; i < bulkVector_.size(); i++)
	{
		for (int j = start; j < end; j++)
		{
			for (int dim = 0; dim < dofDimension_; dim++)
			{
				bulkVector_[i].dofs[j*dofDimension_ + dim] = TotalDofCounter;
				TotalDofCounter++;
			}
		}
	}
	return TotalDofCounter;
}



void Mesh3D::generateStructuredMesh()
{
	creatNode();
	//mesh.printNode();
	creatEdge();
	//mesh.printEdge();
	creatFace();
	creatBulk();
	//mesh.printBulk();
	creatElement();

	assignDofs();
}





/*---------------------------------------------------------
	�˺����ҵ� startPoint �� endPoint ֮�������������
	���裺startPoint �� endPoint ������ʵ�Ľ�㣬λ��ʵ�ʵ����ϣ�
		����startPoint ������С�� endPoint ������
---------------------------------------------------------*/
void Mesh3D::findFaces(std::vector<int> &faces, double* faceStart, double* faceEnd)
{
	int i;
	double CoordinatesDifference[Dimention];
	for (i = 0; i < Dimention; i++)  
		CoordinatesDifference[i] = faceEnd[i] - faceStart[i];

	// ���ݼ��裬CoordinatesDifference ֻ��һ���� ��ֵ
	// FaceType = 0 for yz, 1 for zx, 2 for xy
	double StartEndDistance = abs(CoordinatesDifference[0]);
	int FaceType = 0;

	// �ҵ����Ĳ�ֵ�ͷ���
	double absolute = 0;
	for (i = 0; i < Dimention; i++)
	{
		absolute = abs(CoordinatesDifference[i]);
		if (StartEndDistance > absolute)
		{
			StartEndDistance = absolute;
			FaceType = i;
		}
	}

	// ������ʼ����յ��Ƿ���������
	if (!checkPointIsNode(faceStart))
		printf("Error in face load definition: StartPoint has a problem.\n");
	if (!checkPointIsNode(faceEnd))
		printf("Error in face load definition: EndPoint has a problem.\n");

	// ������ʼ����յ��Ƿ��غ�
	if (StartEndDistance != 0)
		printf("Error in face load definition: StartPoint and EndPoint should be on the same plane (xy, yz or zx)\n");

	// ֻ��һ�������ֵΪ0��������������ͬһ��ֱ����
	double sumCoord = std::accumulate(CoordinatesDifference, CoordinatesDifference + Dimention, 0.0);
	double maxCoord = *std:: max_element(CoordinatesDifference, CoordinatesDifference + Dimention);
	if (sumCoord == maxCoord)
		printf("Error in face load definition: the line between StartPoint and EndPoint should not be parallel to x-, y- or z-axis.\n");


	int *IndicesStart, *IndicesEnd;
	IndicesStart = convertCoordinatesIntoIndices(faceStart);
	IndicesEnd = convertCoordinatesIntoIndices(faceEnd);

	int numOfFacesAlongX = 0;
	int numOfFacesAlongY = 0;
	int numOfFacesAlongZ = 0;
	int faceNumStart = 0;
	int numOfFaces = 0;

	// Define number of xy and yz faces
	int NumberOfXYFaces = numCells[0] * numCells[1] * (numCells[2] + 1);
	int NumberOfYZFaces = numCells[1] * numCells[2] * (numCells[0] + 1);

	// FaceType = 0 for yz, 1 for zx, 2 for xy
	switch (FaceType)
	{
	case 0:
		// 0 for yz
		numOfFacesAlongY = IndicesEnd[1] - IndicesStart[1];
		numOfFacesAlongZ = IndicesEnd[2] - IndicesStart[2];
		numOfFaces = numOfFacesAlongY * numOfFacesAlongZ;

		faceNumStart = NumberOfXYFaces + IndicesStart[1] + numCells[1] * IndicesStart[2] + numCells[1] * numCells[2] * IndicesStart[0];

		faces.resize(numOfFaces);

		for (i = 0; i < numOfFaces; i++)
		{
			// faces.push_back(faceNumStart + i);
			faces[i] = faceNumStart + i;
		}

		break;
	case 1:
		// 1 for xz
		numOfFacesAlongX = IndicesEnd[0] - IndicesStart[0];
		numOfFacesAlongZ = IndicesEnd[2] - IndicesStart[2];
		numOfFaces = numOfFacesAlongX * numOfFacesAlongZ;

		faceNumStart = NumberOfXYFaces + NumberOfYZFaces + IndicesStart[0] + numCells[0] * IndicesStart[2] + numCells[0] * numCells[2] * IndicesStart[1];

		faces.resize(numOfFaces);

		for (i = 0; i < numOfFaces; i++)
		{
			// faces.push_back(faceNumStart + i);
			faces[i] = faceNumStart + i;
		}

		break;
	case 2:

		// 2 for xy
		numOfFacesAlongX = IndicesEnd[0] - IndicesStart[0];
		numOfFacesAlongY = IndicesEnd[1] - IndicesStart[1];
		numOfFaces = numOfFacesAlongX * numOfFacesAlongY;

		faceNumStart = IndicesStart[0] + numCells[0] * IndicesStart[1] + numCells[0] * numCells[1] * IndicesStart[2];

		faces.resize(numOfFaces);

		for (i = 0; i < numOfFaces; i++)
		{
			// faces.push_back(faceNumStart + i);
			faces[i] = faceNumStart + i;
		}
		break;
	}
	delete[]IndicesStart;	IndicesStart = nullptr;
	delete[]IndicesEnd;		IndicesEnd = nullptr;
}



/*---------------------------------------------------------
	�������һ�����Ƿ��ǽ��
---------------------------------------------------------*/

bool Mesh3D::checkPointIsNode(double* Coordinate)
{
	//const int Dim = 3;
	bool IsNode = true;
	double Coords[Dimention];

	// �����������������Ĳ�ֵ
	for (int i = 0; i < Dimention; i++)  
		Coords[i] = Coordinate[i] - meshOrigin[i];

	// �жϸõ��Ƿ� ������
	if (Coords[0] > lateralSize[0] || Coords[1] > lateralSize[1] || Coords[2] > lateralSize[2])
	{
		IsNode = false;
	}		

	// �ж��Ƿ��� ��ֵ
	else if(Coords[0] < 0 || Coords[1] < 0 || Coords[2] < 0)
	{
		IsNode = false;
	}
	
	// �ж��Ƿ��� ���
	else
	{		
		double ii = fmod(Coords[0], lateralSize[0] / numCells[0]);
		double jj = fmod(Coords[1], lateralSize[1] / numCells[1]);
		double kk = fmod(Coords[2], lateralSize[2] / numCells[2]);

		double tolerance = 1e-8;
		if ((ii > tolerance && (lateralSize[0] / numCells[0]) - ii > tolerance)
			|| (jj > tolerance && (lateralSize[1] / numCells[1]) - jj > tolerance)
			|| (kk > tolerance && (lateralSize[2] / numCells[2]) - kk > tolerance))
			IsNode = false;
	}
	return IsNode;
}



int* Mesh3D::convertCoordinatesIntoIndices(double *Coordinate)
{
	double Coords[Dimention];
	// �����������������Ĳ�ֵ
	for (int i = 0; i < Dimention; i++)  Coords[i] = Coordinate[i] - meshOrigin[i];

	// ��ȫ�������ת��Ϊ Indices
	// ceil ��Ϊ��determine element index
	// round ��Ϊ��determine closest node to a point
	int *Indices = new int[Dimention];

	Indices[0] = int(round(Coords[0] * numCells[0] / lateralSize[0]));
	if (Indices[0] > numCells[0]) 
		Indices[0] = numCells[0];

	Indices[1] = int(round(Coords[1] * numCells[1] / lateralSize[1]));
	if (Indices[1] > numCells[1]) 
		Indices[1] = numCells[1];

	Indices[2] = int(round(Coords[2] * numCells[2] / lateralSize[2]));
	if (Indices[2] > numCells[2]) 
		Indices[2] = numCells[2];

	return Indices;
}


/*---------------------------------------------------------
	�ҵ��õ����ڵĵ�Ԫ
---------------------------------------------------------*/
int Mesh3D::findElementByPoint(double *Coords)
{
	double x, y, z;
	int numx, numy, numz;
	int ElementId;
	x = Coords[0] - meshOrigin[0];
	y = Coords[1] - meshOrigin[1];
	z = Coords[2] - meshOrigin[2];

	numx = int(ceil(x * numCells[0] / lateralSize[0]));
	if (numx > numCells[0]) numx = numCells[0];

	numy = int(ceil(y * numCells[1] / lateralSize[1]));
	if (numy > numCells[1]) numy = numCells[1];

	numz = int(ceil(z * numCells[2] / lateralSize[2]));
	if (numz > numCells[2]) numz = numCells[2];

	ElementId = numx-1 + (numy-1) * numCells[0] + (numz-1) * numCells[0] * numCells[1];

	return ElementId;
}








double* Mesh3D::getFaceVertexCoords(int faceId)
{
	double* FVC = new double[4* Mesh3D::Dimention]();
	int nodeId;

	for (int i = 0; i < 4; i++)
	{
		for (int j = 0; j < Mesh3D::Dimention; j++)
		{
			nodeId = faceVector_[faceId].nodesId[i];
			FVC[i * 3 + j] = nodeVector_[nodeId].coords_[j];
			//printf("%6.4f ", FVC[i * 4 + j]);
		}
	}
	return FVC;
}


double* Mesh3D::getEdgeVertexCoords(int edgeId)
{
	double* EVC = new double[2 * Mesh3D::Dimention];
	int nodeId;

	for (int i = 0; i < 2; i++)
	{
		for (int j = 0; j < Mesh3D::Dimention; j++)
		{
			nodeId = edgeVector_[edgeId].nodesId[i];
			EVC[i * 3 + j] = nodeVector_[nodeId].coords_[j];
		}
	}
	return EVC;
}




/*---------------------------------------------------------
	�� ���� �ľֲ�����ӳ�䵽ȫ������
---------------------------------------------------------*/
double* Mesh3D::mapLocalToGlobalForFace(double *Local, double *FVC)
{	
/* ----------���� ������----VC ������ŷ�λ��--------------
						x   y   z
		Vertices 0		0   1   2
		Vertices 1		3   4   5
		Vertices 2		6   7   8
		Vertices 3		9   10  11
-------------------------------------------------*/

	// ����� Local ָ���ǻ��ֵ�����꣬Ҳ���Ǿֲ�����
	// FVC means Face Vertices Coordinates
	double* global = new double[3]();

	const int &Dim = Mesh3D::Dimention;
	for (int i = 0; i < Dim; i++)
	{
		global[i] = 0.25*((1 - Local[0])*(1 - Local[1])*FVC[i]
			+ (1 + Local[0])*(1 - Local[1])*FVC[Dim + i]
			+ (1 + Local[0])*(1 + Local[1])*FVC[2 * Dim + i]
			+ (1 - Local[0])*(1 + Local[1])*FVC[3 * Dim + i]);
	}	

	return global;
}


double Mesh3D::calcFaceDetJ(int faceId)
{
	int edge1 = faceVector_[faceId].edgesId[0];
	int edge2 = faceVector_[faceId].edgesId[1];

	double detJ1 = calcEdgeDetJ(edge1);
	double detJ2 = calcEdgeDetJ(edge2);

	double faceFetJ = detJ1 * detJ2;
	return faceFetJ;
}




double Mesh3D::calcEdgeDetJ(int edgeId)
{

/* ----------���� ��----EVC ������ŷ�λ��--------------
						x   y   z
		Vertices 0		0   1   2
		Vertices 1		3   4   5
-------------------------------------------------*/

	double* EVC = getEdgeVertexCoords(edgeId);

	double x = EVC[3] - EVC[0];
	double y = EVC[4] - EVC[1];
	double z = EVC[5] - EVC[2];

	double detJ = 0.5*sqrt(x*x + y * y + z * z);
	return detJ;

}


