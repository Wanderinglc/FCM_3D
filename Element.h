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
//#include "LegendrePolynomial.h"

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
	Hexahedral(const int polynomialDegree):Element(polynomialDegree),
		locationArray(nullptr)
	{
		nodes.resize(8);
		edges.resize(12);
		faces.resize(6);
		bulks.resize(1);
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
	double* clacFace_dN(double *Coords, double *PhiR, double *PhiS, double *PhiT, double *dPhiR, double *dPhiS, double *dPhiT, double &OneMinusR, double &OneMinusS, double &OneMinusT, double &OnePlusR, double &OnePlusS, double &OnePlusT);
	
	double* calcBulk_N(double *Coords, double *PhiR, double *PhiS, double *PhiT, double &OneMinusR, double &OneMinusS, double &OneMinusT, double &OnePlusR, double &OnePlusS, double &OnePlusT);
	double* clacBulk_dN(double *Coords, double *PhiR, double *PhiS, double *PhiT, double *dPhiR, double *dPhiS, double *dPhiT, double &OneMinusR, double &OneMinusS, double &OneMinusT, double &OnePlusR, double &OnePlusS, double &OnePlusT);

	double* calcShapeN(double *Coords);
	double* calcMatrixB(double *Coords);

public:
	int* locationArray;			//联系数组
	//double* gaussPoints;		//高斯积分点坐标

};



#endif
