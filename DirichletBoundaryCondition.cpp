
#include <iostream>
#include <vector>
#include <algorithm> 

#include "DirichletBoundaryCondition.h"



void FaceDirichletBoundaryCondition::modifyLinearSystem(Mesh3D* mesh, double* K, int* Kdiag, double* F)
{
	// find faces and constrain the face dofs
	std::vector<int> Faces;
	mesh->findFaces(Faces, faceStart, faceEnd);
	constrainFaces(mesh, Faces, K, Kdiag, F, Mesh3D::Dimention);


	// find edges and constrain the edge dofs
	std::vector<int> Edges(4* Faces.size());
	for (size_t i = 0; i < Faces.size(); i++)
	{
		for (size_t j = 0; j < 4; j++)
		{
			Edges[i * 4 + j] = mesh->faceVector_[Faces[i]].edgesId[j];
		}
	}
	std::sort(Edges.begin(), Edges.end());
	std::vector<int>::iterator iterEdge = std::unique(Edges.begin(), Edges.end());
	Edges.erase(iterEdge, Edges.end());
	
	constrainEdges(mesh, Edges, K, Kdiag, F, Mesh3D::Dimention);


	// find nodes and constrain the node dofs
	std::vector<int> Nodes(4 * Faces.size());
	for (size_t i = 0; i < Faces.size(); i++)
	{
		for (size_t j = 0; j < 4; j++)
		{
			Nodes[i * 4 + j] = mesh->faceVector_[Faces[i]].nodesId[j];
		}
	}
	std::sort(Nodes.begin(), Nodes.end());
	std::vector<int>::iterator iterNodes = std::unique(Nodes.begin(), Nodes.end());
	Nodes.erase(iterNodes, Nodes.end());

	constrainEdges(mesh, Nodes, K, Kdiag, F, Mesh3D::Dimention);

}


void DirichletBoundaryCondition::constrainFaces(Mesh3D* mesh, std::vector<int>& faces, double* K, int* Kdiag, double* F, const int& Dim)
{
	int perField = mesh->faceVector_[faces[0]].numberOfDofsPerFieldComponent;
	int ID = 0;

	for (size_t i = 0; i < faces.size(); i++)
	{		
		for (int p = 0; p < perField; p++)
		{
			for (int j = 0; j < Dim; j++)
			{
				if (direction[j] != 0)
				{
					ID = mesh->faceVector_[faces[i]].dofs[p * Dim + j].getId();
					PenaltyAlgorithm->modifyLinearSystem(ID, 0.0, K, Kdiag, F);
				}
			}			
		}
	}

}

void DirichletBoundaryCondition::constrainEdges(Mesh3D* mesh, std::vector<int>& edges, double* K, int* Kdiag, double* F, const int& Dim)
{
	int perField = mesh->edgeVector_[edges[0]].numberOfDofsPerFieldComponent;
	int ID = 0;

	for (size_t i = 0; i < edges.size(); i++)
	{
		for (int p = 0; p < perField; p++)
		{
			for (int j = 0; j < Dim; j++)
			{
				if (direction[j] != 0)
				{
					ID = mesh->edgeVector_[edges[i]].dofs[p * Dim + j].getId();
					PenaltyAlgorithm->modifyLinearSystem(ID, 0.0, K, Kdiag, F);
				}
			}
		}
	}
}

void DirichletBoundaryCondition::constrainNodes(Mesh3D* mesh, std::vector<int>& nodes, double* K, int* Kdiag, double* F, const int& Dim)
{
	int ID = 0;

	for (size_t i = 0; i < nodes.size(); i++)
	{
		for (int j = 0; j < Dim; j++)
		{
			if (direction[j] != 0)
			{
				ID = mesh->nodeVector_[nodes[i]].dofs[j].getId();
				PenaltyAlgorithm->modifyLinearSystem(ID, prescribedValue[j], K, Kdiag, F);
			}
		}
	}
}













