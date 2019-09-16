
#include <iostream>
#include "Mesh3D.h"
#include "Dof.h"
#include "Node.h"
#include "Edge.h"
#include "Face.h"
#include "Bulk.h"
#include "Element.h"



Mesh3D::Mesh3D()
{

}

Mesh3D::~Mesh3D()
{

}

/*---------------------------------------
			生成结点
  ---------------------------------------*/
void Mesh3D::creatNode()
{
	//int dCell[3];
	double xIncrement = (lateralSize[0]) / numCells[0];		// x 方向增量
	double yIncrement = (lateralSize[1]) / numCells[1];		// y 方向增量
	double zIncrement = (lateralSize[2]) / numCells[2];		// y 方向增量

	int numNodes[3];											// x,y,z 三个方向上的结点数
	for (int i = 0; i < 3; i++)
		numNodes[i] = numCells[i] + 1;

	totalNodes_ = numNodes[0] * numNodes[1] * numNodes[2];		//总结点数

	// Initialize nodeVector_
	for (int k = 0; k < numNodes[2]; k++)
	{
		for (int j = 0; j < numNodes[1]; j++)
		{
			for (int i = 0; i < numNodes[0]; i++)
			{
				Node node(dofDimension_);
				node.coords_[0] = MeshOrigin[0] + i * xIncrement;
				node.coords_[1] = MeshOrigin[1] + j * yIncrement;
				node.coords_[2] = MeshOrigin[2] + k * zIncrement;
				node.id = i + j * numNodes[0] + k * numNodes[0]* numNodes[1];
				nodeVector_.push_back(node);
			}
		}
	}

}
//		输出结点信息
void Mesh3D::printNode() const
{
	std::cout << "----------Print Nodes----------\n";
	int numNodes[3];											// x,y,z 三个方向上的结点数
	for (int i = 0; i < 3; i++)		numNodes[i] = numCells[i] + 1;

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
		printf("Node %4d : ( %5.1f, %5.1f, %5.1f )\n", 
			node.id, node.coords_[0], node.coords_[1], node.coords_[2]);
	}
}



/*---------------------------------------
			生成 边
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
	int numNodesInXY = (numCells[0] + 1)*(numCells[1] + 1);		//xy 平面内的结点数
	int count = 0;

	int leftNode, rightNode;
	for (int k = 0; k < numCells[2] + 1; k++)
	{
		for (int j = 0; j < numCells[1] + 1; j++)
		{
			for (int i = 0; i < numCells[0]; i++)
			{
				Edge edg(polynomialDegree_, dofDimension_);
				// 边的左右两个结点编号
				edg.id = count; 
				count++;
				leftNode = k * (numNodesInXY)+j * (numCells[0] + 1) + i;
				rightNode = leftNode + 1;
				edg.nodes[0] = &nodeVector_[leftNode];
				edg.nodes[1] = &nodeVector_[rightNode];
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
				// 边的左右两个结点编号
				edg.id = count;
				count++;
				leftNode = k * (numNodesInXY)+j * (numCells[0] + 1) + i;
				rightNode = leftNode + numCells[0] + 1;
				edg.nodes[0] = &nodeVector_[leftNode];
				edg.nodes[1] = &nodeVector_[rightNode];
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
				// 边的左右两个结点编号
				edg.id = count;
				count++;
				leftNode = k * (numNodesInXY)+j * (numCells[0] + 1) + i;
				rightNode = leftNode + numNodesInXY;
				edg.nodes[0] = &nodeVector_[leftNode];
				edg.nodes[1] = &nodeVector_[rightNode];
				edgeVector_.push_back(edg);
			}
		}
	}

}
//        输出 边 信息
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
		//printf("%4d\t%4d\n",edg->nodeId[0], edg->nodeId[1]);
		printf("Edge %d : %4.1f, %4.1f, %4.1f\n", edg->id, edg->nodes[0]->coords_[0],
			edg->nodes[0]->coords_[1], edg->nodes[0]->coords_[2]);
	}
	printf("\n");

}

/*---------------------------------------
			生成面
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
				fa.edges[0] = &edgeVector_[edge1];
				fa.edges[1] = &edgeVector_[edge2];
				fa.edges[2] = &edgeVector_[edge3];
				fa.edges[3] = &edgeVector_[edge4];

				fa.nodes[0] = fa.edges[0]->nodes[0];
				fa.nodes[1] = fa.edges[1]->nodes[0];
				fa.nodes[2] = fa.edges[2]->nodes[1];
				fa.nodes[3] = fa.edges[3]->nodes[1];

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
				fa.edges[0] = &edgeVector_[edge1];
				fa.edges[1] = &edgeVector_[edge2];
				fa.edges[2] = &edgeVector_[edge3];
				fa.edges[3] = &edgeVector_[edge4];

				fa.nodes[0] = fa.edges[0]->nodes[0];
				fa.nodes[1] = fa.edges[1]->nodes[0];
				fa.nodes[2] = fa.edges[2]->nodes[1];
				fa.nodes[3] = fa.edges[3]->nodes[1];

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
				fa.edges[0] = &edgeVector_[edge1];
				fa.edges[1] = &edgeVector_[edge2];
				fa.edges[2] = &edgeVector_[edge3];
				fa.edges[3] = &edgeVector_[edge4];

				fa.nodes[0] = fa.edges[0]->nodes[0];
				fa.nodes[1] = fa.edges[1]->nodes[0];
				fa.nodes[2] = fa.edges[2]->nodes[1];
				fa.nodes[3] = fa.edges[3]->nodes[1];

				//for (int n = 0; n < 4; n++)
				//	printf("-----node id %d: %4.1f, %4.1f, %4.1f\n", fa.nodes[n]->id, fa.nodes[n]->coords_[0], fa.nodes[n]->coords_[1],
				//		fa.nodes[n]->coords_[2]);;
				faceVector_.push_back(fa);

			}
		}
	}
}

/*---------------------------------------
			生成体
  ---------------------------------------*/
void Mesh3D::creatBulk()
{
	totalBulks_ = numCells[0] * numCells[1] * numCells[2];
	int NodeNum[8];
	int EdgeNum[12];
	int FaceNum[6];

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

				// 组成体的 结点 编号
				NodeNum[0] = i + j * (numCells[0] + 1) + k * numNodesInXY;
				NodeNum[1] = NodeNum[0] + 1;
				NodeNum[3] = NodeNum[0] + numCells[0] + 1;
				NodeNum[2] = NodeNum[3] + 1;
				NodeNum[4] = NodeNum[0] + numNodesInXY;
				NodeNum[5] = NodeNum[4] + 1;
				NodeNum[7] = NodeNum[4] + numCells[0] + 1;
				NodeNum[6] = NodeNum[7] + 1;
				for (int n = 0; n < 8; n++)
				{
					cube.nodes[n] = &nodeVector_[NodeNum[n]];
				}
				// 组成体的 边 编号
				EdgeNum[0] = i + j * numCells[0] + k * numCells[0] * (numCells[1] + 1);
				EdgeNum[2] = EdgeNum[0] + numCells[0];
				EdgeNum[3] = NumOfXEdges + j + i * numCells[1] + k * numCells[1] * (numCells[0] + 1);
				EdgeNum[1] = EdgeNum[3] + numCells[1];

				EdgeNum[8] = EdgeNum[0] + numCells[0] * (numCells[1] + 1);
				EdgeNum[10] = EdgeNum[8] + numCells[0];
				EdgeNum[11] = EdgeNum[3] + numCells[1] * (numCells[0] + 1);
				EdgeNum[9] = EdgeNum[11] + numCells[1];

				EdgeNum[4] = NumOfXEdges + NumOfYEdges + k + i * numCells[2] + j * numCells[2] * (numCells[0] + 1);
				EdgeNum[5] = EdgeNum[4] + numCells[2];
				EdgeNum[7] = EdgeNum[4] + numCells[2] * (numCells[0] + 1);
				EdgeNum[6] = EdgeNum[7] + numCells[2];
				for (int n = 0; n < 12; n++)
				{
					cube.edges[n] = &edgeVector_[EdgeNum[n]];
				}

				// 组成体的 面 编号
				FaceNum[0] = i + j * numCells[0] + k * numCells[0] * numCells[1];
				FaceNum[5] = FaceNum[0] + numCells[0] * numCells[1];
				FaceNum[4] = NumOfXYFaces + j + k * numCells[1] + i * numCells[1] * numCells[2];
				FaceNum[2] = FaceNum[4] + numCells[1] * numCells[2];
				FaceNum[1] = NumOfXYFaces + NumOfYZFaces + i + k * numCells[0] + j*numCells[0] * numCells[2];
				FaceNum[3] = FaceNum[1] + numCells[0] * numCells[2];
				
				//printf(" Face Num: ");
				for (int n = 0; n < 6; n++)
				{
					//printf(" %d ", FaceNum[n]);
					cube.faces[n] = &faceVector_[FaceNum[n]];
				}
				//printf("\n");
				bulkVector_.push_back(cube);

			}
		}
	}

}
void Mesh3D::printBulk()
{
	std::cout << "----------Print Bulk----------\n";

	for (int i = 0; i < bulkVector_.size(); i++) 
	{
		const Bulk* cube = &bulkVector_[i];
		printf("\nCube %d: \n", cube->id);
		for (int n = 0; n < 8; n++)
		{
			printf("   node %d: %4.1f, %4.1f, %4.1f\n", cube->nodes[n]->id,
				cube->nodes[n]->coords_[0], cube->nodes[n]->coords_[1],
				cube->nodes[n]->coords_[2]);
		}

		printf("   \n");
		for (int e = 0; e < 12; e++)
		{
			printf("   middle point of edge %d: %4.1f, %4.1f, %4.1f\n", cube->edges[e]->id,
				(cube->edges[e]->nodes[0]->coords_[0]+ cube->edges[e]->nodes[1]->coords_[0])/2,
				(cube->edges[e]->nodes[0]->coords_[1] + cube->edges[e]->nodes[1]->coords_[1]) / 2, 
				(cube->edges[e]->nodes[0]->coords_[2] + cube->edges[e]->nodes[1]->coords_[2]) / 2 );
		}

		printf("   \n");
		for (int f = 0; f < 6; f++)
		{
			printf("   center point of face %d: %4.1f, %4.1f, %4.1f\n", cube->faces[f]->id,
				(cube->faces[f]->nodes[0]->coords_[0] + cube->faces[f]->nodes[2]->coords_[0]) / 2,
				(cube->faces[f]->nodes[0]->coords_[1] + cube->faces[f]->nodes[2]->coords_[1]) / 2,
				(cube->faces[f]->nodes[0]->coords_[2] + cube->faces[f]->nodes[2]->coords_[2]) / 2);
		}

	}


}

/*---------------------------------------
			生成单元
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
		// 单元的体
		hex.bulks[0] = &bulkVector_[i];


		// 单元的结点
		//printf("\n    Element Node: ");
		for (int iNode = 0; iNode < 8; iNode++)
		{
			hex.nodes[iNode] = bulkVector_[i].nodes[iNode];
			//printf("%d  ",hex.nodes[iNode]->id );
		}
			
		
		// 单元的边
		//printf("\n    Element Edge: ");
		for (int iEdge = 0; iEdge < 12; iEdge++)
		{
			hex.edges[iEdge] = bulkVector_[i].edges[iEdge];
			//printf("%d  ", hex.edges[iEdge]->id);
		}
			

		// 单元的面
		//printf("\n    Element face: ");
		for (int iFace = 0; iFace < 6; iFace++)
		{
			hex.faces[iFace] = bulkVector_[i].faces[iFace];
			//printf("%d  ", hex.faces[iFace]->id);
		}
		//printf("\n");
		
		elementVector_.push_back(hex);
	}


}


/*---------------------------------------
			自由度编号
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
			totalDofCounter = assignNodesToEdges(p, totalDofCounter);	// 给边上结点编号
			//printf("total Dof Counter %d\n", totalDofCounter);
			totalDofCounter = assignNodesToFaces(p, totalDofCounter);	// 给面内结点编号
			//printf("total Dof Counter %d\n", totalDofCounter);
			totalDofCounter = assignNodesToBulks(p, totalDofCounter);
			//printf("total Dof Counter %d\n", totalDofCounter);
		}
	}


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





