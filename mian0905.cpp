#include <iostream>
#include <vector>

#include "Mesh3D.h"
#include "Edge.h"
#include "Node.h"

int main()
{
	Mesh3D mesh;
	mesh.MeshOrigin[0] = 0.0;
	mesh.MeshOrigin[1] = 0.0;
	mesh.MeshOrigin[2] = 0.0;
	mesh.lateralSize[0] = 2.0;
	mesh.lateralSize[1] = 2.0;
	mesh.lateralSize[2] = 1.0;
	mesh.numCells[0] = 2;
	mesh.numCells[1] = 2;
	mesh.numCells[2] = 2;
	mesh.polynomialDegree_ = 3;
	mesh.dofDimension_ = 3;

	mesh.creatNode();
	//mesh.printNode();
	mesh.creatEdge();
	//mesh.printEdge();
	mesh.creatFace();

	mesh.creatBulk();
	//mesh.printBulk();
	mesh.creatElement();

	mesh.assignDofs();

	mesh.elementVector_[0].getLocationMatrix();


	getchar();
	return 0;
}


