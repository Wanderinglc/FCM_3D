#ifndef ELEMENT_H
#define ELEMENT_H

/*=======================================================================
					  Element
	 Define the element class base for FEM
  =======================================================================*/
#include <vector>

#include "Node.h"
#include "Edge.h"
#include "Face.h"
#include "Bulk.h"
#include "Material.h"
//#include "LegendrePolynomial.h"
class EmbeddedDomain;
class GaussIntegrator;

class Element
{
public:
	// Constructor
	Element() { }
	Element(const int polynomialDegree):polynomialDegree(polynomialDegree){}

	// Destructor
	virtual ~Element(){}



	/*-------------------------
		Basic Attributes
	  -------------------------*/
	int id;
	int polynomialDegree;

	AbsMaterial* material;

	std::vector<Node*> nodes;
	std::vector<Edge*> edges;
	std::vector<Face*> faces;
	std::vector<Bulk*> bulks;

	
};


class Hexahedral :public Element
{
	//  Vertices
	//
	//       8 +-------------------+ 7
	//        /|                  /|
	//       / |                 / |
	//      /  |                /  |
	//     /   |               /   |
	//    /    |            6 /    |
	// 5 +-----|-------------+     |          
	//   |     |             |     |
	//   |     +-------------|-----+ 3          
	//   |    / 4            |    /
	//   |   /               |   /
	//   |  /                |  / 
	//   | /                 | /
	//   |/			         |/
	//   +-------------------+
	//   1                   2
	//
	//
	//  Edges
	//         +---------11--------+  
	//        /|                  /|
	//       / |                 / |
	//     12  |               10  |
	//     /   8               /   7
	//    /    |              /    |
	//   +-----|---9---------+     |          
	//   |     |             |     |
	//   |     +---------3---|-----+           
	//   |    /              |    /
	//   5   /               6   /
	//   |  4                |  2 
	//   | /                 | /
	//   |/			         |/
	//   +---------1---------+


public:
	Hexahedral() { }
	Hexahedral(const int p):Element(p),
		locationArray(nullptr)
	{
		nodes.resize(NumNodes);
		edges.resize(NumEdges);
		faces.resize(NumFaces);
		bulks.resize(NumBulks);
	}

	~Hexahedral()
	{
		nodes.clear();
		edges.clear();
		faces.clear();
		bulks.clear();
	}

	// 联系数组
	int* getLocationMatrix();
	double calcPhi(int p, double Coord);				//高阶多项式 Phi
	double calcDerivOfPhi(int p, double Coord);		//高阶多项式 Phi的导数

	double* calcVertex_N(double *Coords, double &OneMinusR, double &OneMinusS, double &OneMinusT, double &OnePlusR, double &OnePlusS, double &OnePlusT);
	double* calcVertex_dN(double *Coords, double &OneMinusR, double &OneMinusS, double &OneMinusT, double &OnePlusR, double &OnePlusS, double &OnePlusT);

	double* calcEdge_N(double *Coords, double *PhiR, double *PhiS, double *PhiT, double &OneMinusR, double &OneMinusS, double &OneMinusT, double &OnePlusR, double &OnePlusS, double &OnePlusT);
	double* calcEdge_dN(double *Coords, double *PhiR, double *PhiS, double *PhiT, double *dPhiR, double *dPhiS, double *dPhiT, double &OneMinusR, double &OneMinusS, double &OneMinusT, double &OnePlusR, double &OnePlusS, double &OnePlusT);
	
	double* calcFace_N(double *Coords, double *PhiR, double *PhiS, double *PhiT, double &OneMinusR, double &OneMinusS, double &OneMinusT, double &OnePlusR, double &OnePlusS, double &OnePlusT);
	double* calcFace_dN(double *Coords, double *PhiR, double *PhiS, double *PhiT, double *dPhiR, double *dPhiS, double *dPhiT, double &OneMinusR, double &OneMinusS, double &OneMinusT, double &OnePlusR, double &OnePlusS, double &OnePlusT);
	
	double* calcBulk_N(double *Coords, double *PhiR, double *PhiS, double *PhiT, double &OneMinusR, double &OneMinusS, double &OneMinusT, double &OnePlusR, double &OnePlusS, double &OnePlusT);
	double* calcBulk_dN(double *Coords, double *PhiR, double *PhiS, double *PhiT, double *dPhiR, double *dPhiS, double *dPhiT, double &OneMinusR, double &OneMinusS, double &OneMinusT, double &OnePlusR, double &OnePlusS, double &OnePlusT);

	double* calcShapeN(double *Coords);
	void printShapeN(double *Coords);
	double* calcShape_dN(double *Coords);
	void printShape_dN(double *Coords);

	double* calcMatrixB(double *Coords);
	void printMatrixB(double *Coords);

	double* calcJacobiAtCenter();
	double calcDetJac(double *VC);
	double calcSubCellDetJac(const int &SubDepth);

	double* calcStiffnessMatrix(EmbeddedDomain *pEmDo, AbsMaterial *pMaterial, GaussIntegrator *pIntegrator, const int Depth);
	void calcUndividedStiffness(AbsMaterial *pMaterial, GaussIntegrator *pIntegrator, double *Kc, double CellDetJ);
	void CalcSubCellStifiness(EmbeddedDomain *pEmDo, AbsMaterial *pMaterial, GaussIntegrator *pIntegrator, double *GlobalVertexCoords, double *LocalVertexCoords,
		const int depthCounter, double *Kc, const double CellDetJ);

	void printStiffnessMatrix(double* K);

	double* getCellVertexCoords() const;
	void printCellVertexCoords();
	double* mapLocalToGlobal(double *Local, double *VerCoo);
	double* mapGlobalToLocal(double *Global);
	bool isOutside(EmbeddedDomain *pEmDo, double* VerCoo);
	bool isInside(EmbeddedDomain *pEmDo, double* VerCoo);
	bool isIntersected(EmbeddedDomain *pEmDo, double* VerCoo);
	void partitionRecursively(EmbeddedDomain *pEmDo, AbsMaterial *pMaterial, GaussIntegrator *pIntegrator, double *GlobalVerCoo, double *LocalVerCoo,
		int depthCounter, const int Depth, double *Kc, double CellDetJ);




public:
	static const int NumNodes = 8;
	static const int NumEdges = 12;
	static const int NumFaces = 6;
	static const int NumBulks = 1;
	//int nodesId[8];
	//int edgesId[12];
	//int facesId[6];
	//int bulksId[1];


	int* locationArray;			//联系数组
	//double* gaussPoints;		//高斯积分点坐标

};


/*--------------------------------------------------------------------------------------
		|  dN1_dx  dN2_x  dN3_x  ...  edge dN_dx  ...  face dN_dx  ...  bulk dN_dx  |
  dN =	|  dN1_dy  dN2_y  dN3_y  ...  edge dN_dy  ...  face dN_dy  ...  bulk dN_dy  |
		|  dN1_dz  dN2_z  dN3_z  ...  edge dN_dz  ...  face dN_dz  ...  bulk dN_dz  |
--------------------------------------------------------------------------------------*/

// 计算形函数
/*-----------------------------------------------------------------------------------------------------------------
		|  dN1_dx     0       0	    dN2_dx     0       0	...  edge dN  ...  face dN  ...  bulk dN  |
		|     0    dN1_dy     0	       0    dN2_dy     0	...  edge dN  ...  face dN  ...  bulk dN  |
	B =	|     0       0    dN1_dz	   0       0    dN2_dz  ...  edge dN  ...  face dN  ...  bulk dN  |
		|  dN1_dy  dN1_dx     0     dN2_dy  dN2_dx     0    ...  edge dN  ...  face dN  ...  bulk dN  |
		|     0    dN1_dz  dN1_dy      0    dN2_dz  dN2_dy  ...  edge dN  ...  face dN  ...  bulk dN  |
		|  dN1_dz    0     dN1_dx   dN2_dz    0     dN1_dx  ...  edge dN  ...  face dN  ...  bulk dN  |
-----------------------------------------------------------------------------------------------------------------*/

/*===================================================================
% Vertices
				+-------------------+-------------------+
			   /|24                /|25                /|26
			  / |                 / |                 / |
			 /  |                /  |                /  |
			/   |               /   |               /   |
		   /    |              /    |              /    |
		  +-------------------+-------------------+     |
		 /|21   |            /|22   |            /|23   |
		/ |     +-----------/-------+-----------/-|-----+
	   /  |    /|15        /  |    /|16        /  |    /|17
	  /   |   / |         /   |   / |         /   |   / |
	 /    |  /  |        /    |  /  |        /    |  /  |
	+-------------------+-------------------+     | /   |
	|18   |/    |       |19   |/    |       |20   |/    |
	|     +-------------|-----+-------------|-----+     |
	|	 /|12   |       |    /|13   |       |    /|14   |
	|   / |     +-------|---/-|-----+-------|---/-|-----+
	|  /  |    / 6      |  /  |    / 7      |  /  |    / 8
	| /   |   /         | /   |   /         | /   |   /
	|/    |  /          |/    |  /          |/    |  /
	+-------------------+-------------------+     | /
	|9    |/            |10   |/            |11   |/
	|     +-------------|-----+-------------|-----+
	|    / 3            |    / 4            |    / 5
	|   /               |   /               |   /
	|  /                |  /                |  /
	| /                 | /                 | /
	|/			        |/                  |/
	+-------------------+-------------------+
	0                   1                   2
===================================================================*/





#endif
