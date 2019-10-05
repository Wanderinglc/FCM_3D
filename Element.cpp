
#include <iostream>
#include <fstream>
#include <mkl.h>

#include "Element.h"
#include "LegendrePolynomial.h"
#include "Mesh3D.h"
#include "EmbeddedDomain.h"
#include "GaussIntegrator.h"


int* Hexahedral::getLocationMatrix()
{
	int numOfElementDofs = 3 * (polynomialDegree + 1)*(polynomialDegree + 1)*(polynomialDegree + 1);
	locationArray = new int[numOfElementDofs];

	int countDof = 0;
	for (int i = 0; i < 8; i++)
	{
		locationArray[3 * i] = nodes[i]->dofs[0].id;
		locationArray[3 * i + 1] = nodes[i]->dofs[1].id;
		locationArray[3 * i + 2] = nodes[i]->dofs[2].id;
	}
	countDof = 3 * 8;

	// 边上自由度编号
	int index = 0;
	for (int i = 0; i < 12; i++)
	{
		for (int p = 0; p < polynomialDegree - 1; p++)
		{
			index = countDof + 3 * i*(polynomialDegree - 1) + 3 * p;
			locationArray[index] = edges[i]->dofs[3 * p].id;
			locationArray[index + 1] = edges[i]->dofs[3 * p + 1].id;
			locationArray[index + 2] = edges[i]->dofs[3 * p + 2].id;
		}

	}	
	
	// 面上自由度编号
	countDof = countDof + 3 * 12 * (polynomialDegree - 1);
	int FaceNode = (polynomialDegree - 1)*(polynomialDegree - 1);
	for (int i = 0; i < 6; i++)
	{
		for (int p = 0; p < FaceNode; p++)
		{
			index = countDof + 3 * i*FaceNode + 3 * p;
			locationArray[index] = faces[i]->dofs[3 * p].id;
			locationArray[index + 1] = faces[i]->dofs[3 * p + 1].id;
			locationArray[index + 2] = faces[i]->dofs[3 * p + 2].id;
		}
	}	
	
	// 体内自由度编号
	countDof = countDof + 3 * 6 * (polynomialDegree - 1)*(polynomialDegree - 1);
	int BulkNode = (polynomialDegree - 1)*(polynomialDegree - 1)*(polynomialDegree - 1);

	for (int p = 0; p < BulkNode; p++)
	{
		index = countDof + 3 * p;
		locationArray[index] = bulks[0]->dofs[3 * p].id;
		locationArray[index + 1] = bulks[0]->dofs[3 * p + 1].id;
		locationArray[index + 2] = bulks[0]->dofs[3 * p + 2].id;
	}	

	//printf("\nLocation Matrix: \n");
	//for (int i = 0; i < numOfElementDofs; i++)
	//{
	//	if (i == 0)
	//		printf("\nNode Dof Num:\n");
	//	printf("%d ", locationArray[i]);
	//	if (i==23) 
	//		printf("\n\nEdge Dof Num:\n");
	//	if(i==23+ 3 * 12 * (polynomialDegree - 1))
	//		printf("\n\nFace Dof Num:\n");
	//	if (i == 23 + 3 * 12 * (polynomialDegree - 1)+ 3 * 6 * (polynomialDegree - 1)*(polynomialDegree - 1))
	//		printf("\n\nBulk Dof Num:\n");
	//}

	return locationArray;
}





//----------------------------高阶多项式 Phi--------------------------------
double Hexahedral::calcPhi(int p, double Coord)
{
	double Phi = 1 / sqrt(4 * p - 2);
	Phi = Phi * (getLegendrePolynomial(Coord, p) - getLegendrePolynomial(Coord, p - 2));
	return Phi;
}

//----------------------------高阶多项式 Phi的导数--------------------------------
double Hexahedral::calcDerivOfPhi(int p, double Coord)
{
	double dPhi = 1 / sqrt(4 * p - 2);
	dPhi = dPhi * (getLegendrePolynomialDerivative(Coord, p) - getLegendrePolynomialDerivative(Coord, p - 2));
	return dPhi;
}

//----------------------------角结点的形函数--------------------------------
double* Hexahedral::calcVertex_N(double *Coords, double &OneMinusR, double &OneMinusS, double &OneMinusT, double &OnePlusR, double &OnePlusS, double &OnePlusT)
{
	double *V_N = new double[8];

	V_N[0] = 0.125 * OneMinusR* OneMinusS* OneMinusT;		// N1
	V_N[1] = 0.125 * OnePlusR* OneMinusS* OneMinusT;		// N2
	V_N[2] = 0.125 * OnePlusR* OnePlusS* OneMinusT;			// N3
	V_N[3] = 0.125 * OneMinusR* OnePlusS* OneMinusT;		// N4
	V_N[4] = 0.125 * OneMinusR* OneMinusS* OnePlusT;		// N5
	V_N[5] = 0.125 * OnePlusR* OneMinusS* OnePlusT;			// N6
	V_N[6] = 0.125 * OnePlusR* OnePlusS* OnePlusT;			// N7
	V_N[7] = 0.125 * OneMinusR* OnePlusS* OnePlusT;			// N8
	return V_N;
}

//----------------------------角结点的形函数偏导数--------------------------------
// 先是对方向 1 的偏导数，然后是对方向 2 的偏导数，最后是对方向 3 的偏导数
double* Hexahedral::calcVertex_dN(double *Coords, double &OneMinusR, double &OneMinusS, double &OneMinusT, double &OnePlusR, double &OnePlusS, double &OnePlusT)
{
	double *V_dN = new double[8 * 3];

	V_dN[0] = -0.125 * OneMinusS* OneMinusT;		// dN1dr
	V_dN[1] = -V_dN[0];								// dN2dr
	V_dN[2] =  0.125 * OnePlusS* OneMinusT;			// dN3dr
	V_dN[3] = -V_dN[2];								// dN4dr
	V_dN[4] = -0.125 * OneMinusS* OnePlusT;			// dN5dr
	V_dN[5] = -V_dN[4];								// dN6dr
	V_dN[6] =  0.125 * OnePlusS*OnePlusT;			// dN7dr
	V_dN[7] = -V_dN[6];								// dN8dr


	V_dN[8]  = -0.125 * OneMinusR* OneMinusT;		// dN1ds
	V_dN[9]  = -0.125 * OnePlusR* OneMinusT;		// dN2ds
	V_dN[10] = -V_dN[9];							// dN3ds
	V_dN[11] = -V_dN[8];							// dN4ds
	V_dN[12] = -0.125 * OneMinusR* OnePlusT;		// dN5ds
	V_dN[13] = -0.125 * OnePlusR* OnePlusT;			// dN6ds
	V_dN[14] = -V_dN[13];							// dN7ds
	V_dN[15] = -V_dN[12];							// dN8ds

	V_dN[16] = -0.125 * OneMinusR* OneMinusS;		// dN1dt
	V_dN[17] = -0.125 * OnePlusR* OneMinusS;		// dN2dt
	V_dN[18] = -0.125 * OnePlusR* OnePlusS;			// dN3dt
	V_dN[19] = -0.125 * OneMinusR* OnePlusS;		// dN4dt
	V_dN[20] = -V_dN[16];							// dN5dt
	V_dN[21] = -V_dN[17];							// dN6dt
	V_dN[22] = -V_dN[18];							// dN7dt
	V_dN[23] = -V_dN[19];							// dN8dt

	return V_dN;
}

//----------------------------边结点的形函数--------------------------------
double* Hexahedral::calcEdge_N(double *Coords, double *PhiR, double *PhiS, double *PhiT, double &OneMinusR, double &OneMinusS, double &OneMinusT, double &OnePlusR, double &OnePlusS, double &OnePlusT)
{
	int PerEdge = polynomialDegree - 1;
	double *E_N = new double[12 * PerEdge];

	for (int p = 0; p < PerEdge; p++)
	{
		E_N[p] = 0.25*OneMinusS*OneMinusT*PhiR[p];						//Edge 1
		E_N[PerEdge + p] = 0.25*OnePlusR*OneMinusT*PhiS[p];				//Edge 2
		E_N[2 * PerEdge + p]  = 0.25*OnePlusS*OneMinusT*PhiR[p];		//Edge 3
		E_N[3 * PerEdge + p]  = 0.25*OneMinusR*OneMinusT*PhiS[p];		//Edge 4

		E_N[4 * PerEdge + p]  = 0.25*OneMinusR*OneMinusS*PhiT[p];		//Edge 5
		E_N[5 * PerEdge + p]  = 0.25*OnePlusR*OneMinusS*PhiT[p];		//Edge 6
		E_N[6 * PerEdge + p]  = 0.25*OnePlusR*OnePlusS*PhiT[p];			//Edge 7
		E_N[7 * PerEdge + p]  = 0.25*OneMinusR*OnePlusS*PhiT[p];		//Edge 8

		E_N[8 * PerEdge + p]  = 0.25*OneMinusS*OnePlusT*PhiR[p];		//Edge 9
		E_N[9 * PerEdge + p]  = 0.25*OnePlusR*OnePlusT*PhiS[p];			//Edge 10
		E_N[10 * PerEdge + p] = 0.25*OnePlusS*OnePlusT*PhiR[p];			//Edge 11
		E_N[11 * PerEdge + p] = 0.25*OneMinusR*OnePlusT*PhiS[p];		//Edge 12
	}
	return E_N;

}

//----------------------------边结点的形函数偏导数--------------------------------
// 先是对方向 1 的偏导数，然后是对方向 2 的偏导数，最后是对方向 3 的偏导数
double* Hexahedral::calcEdge_dN(double *Coords, double *PhiR, double *PhiS, double *PhiT, double *dPhiR, double *dPhiS, double *dPhiT, double &OneMinusR, double &OneMinusS, double &OneMinusT, double &OnePlusR, double &OnePlusS, double &OnePlusT)
{
	int PerEdge = polynomialDegree - 1;
	double *E_dN = new double[12 * PerEdge * 3];

	for (int p = 0; p < PerEdge; p++)
	{
		// edge modes derivatives with respect to r
		E_dN[p] = 0.25*OneMinusS*OneMinusT*dPhiR[p];						//Edge 1 wrt r
		E_dN[PerEdge + p] = 0.25*OneMinusT*PhiS[p];							//Edge 2 wrt r
		E_dN[2 * PerEdge + p] =  0.25*OnePlusS*OneMinusT*dPhiR[p];			//Edge 3 wrt r
		E_dN[3 * PerEdge + p] = -0.25*OneMinusT*PhiS[p];					//Edge 4 wrt r

		E_dN[4 * PerEdge + p] = -0.25*OneMinusS*PhiT[p];					//Edge 5 wrt r
		E_dN[5 * PerEdge + p] =  0.25*OneMinusS*PhiT[p];					//Edge 6 wrt r
		E_dN[6 * PerEdge + p] =  0.25*OnePlusS*PhiT[p];						//Edge 7 wrt r
		E_dN[7 * PerEdge + p] = -0.25*OnePlusS*PhiT[p];						//Edge 8 wrt r

		E_dN[8 * PerEdge + p] = 0.25*OneMinusS*OnePlusT*dPhiR[p];			//Edge 9 wrt r
		E_dN[9 * PerEdge + p] = 0.25*OnePlusT*PhiS[p];						//Edge 10 wrt r
		E_dN[10 * PerEdge + p] = 0.25*OnePlusS*OnePlusT*dPhiR[p];			//Edge 11 wrt r
		E_dN[11 * PerEdge + p] = -0.25*OnePlusT*PhiS[p];					//Edge 12 wrt r

		// edge modes derivatives with respect to s
		E_dN[12 * PerEdge + p] = -0.25*OneMinusT*PhiR[p];					//Edge 1 wrt s
		E_dN[13 * PerEdge + p] = 0.25*OnePlusR*OneMinusT*dPhiS[p];			//Edge 2 wrt s
		E_dN[14 * PerEdge + p] = 0.25*OneMinusT*PhiR[p];					//Edge 3 wrt s
		E_dN[15 * PerEdge + p] = 0.25*OneMinusR*OneMinusT*dPhiS[p];			//Edge 4 wrt s

		E_dN[16 * PerEdge + p] = -0.25*OneMinusR*PhiT[p];					//Edge 5 wrt s
		E_dN[17 * PerEdge + p] = -0.25*OnePlusR*PhiT[p];					//Edge 6 wrt s
		E_dN[18 * PerEdge + p] = 0.25*OnePlusR*PhiT[p];						//Edge 7 wrt s
		E_dN[19 * PerEdge + p] = 0.25*OneMinusR*PhiT[p];					//Edge 8 wrt s

		E_dN[20 * PerEdge + p] = -0.25*OnePlusT*PhiR[p];					//Edge 9 wrt s
		E_dN[21 * PerEdge + p] = 0.25*OnePlusR*OnePlusT*dPhiS[p];			//Edge 10 wrt s
		E_dN[22 * PerEdge + p] = 0.25*OnePlusT*PhiR[p];						//Edge 11 wrt s
		E_dN[23 * PerEdge + p] = 0.25*OneMinusR*OnePlusT*dPhiS[p];			//Edge 12 wrt s

		// edge modes derivatives with respect to t
		E_dN[24 * PerEdge + p] = - 0.25*OneMinusS*PhiR[p];					//Edge 1 wrt t
		E_dN[25 * PerEdge + p] = -0.25*OnePlusR*PhiS[p];					//Edge 2 wrt t
		E_dN[26 * PerEdge + p] = -0.25*OnePlusS*PhiR[p];					//Edge 3 wrt t
		E_dN[27 * PerEdge + p] = -0.25*OneMinusR*PhiS[p];					//Edge 4 wrt t

		E_dN[28 * PerEdge + p] = 0.25*OneMinusR*OneMinusS*dPhiT[p];			//Edge 5 wrt t
		E_dN[29 * PerEdge + p] = 0.25*OnePlusR*OneMinusS*dPhiT[p];			//Edge 6 wrt t
		E_dN[30 * PerEdge + p] = 0.25*OnePlusR*OnePlusS*dPhiT[p];			//Edge 7 wrt t
		E_dN[31 * PerEdge + p] = 0.25*OneMinusR*OnePlusS*dPhiT[p];			//Edge 8 wrt t

		E_dN[32 * PerEdge + p] = 0.25*OneMinusS*PhiR[p];					//Edge 9 wrt t
		E_dN[33 * PerEdge + p] = 0.25*OnePlusR*PhiS[p];						//Edge 10 wrt t
		E_dN[34 * PerEdge + p] = 0.25*OnePlusS*PhiR[p];						//Edge 11 wrt t
		E_dN[35 * PerEdge + p] = 0.25*OneMinusR*PhiS[p];					//Edge 12 wrt t

	}
	return E_dN;
}

//----------------------------面结点的形函数--------------------------------
double* Hexahedral::calcFace_N(double *Coords, double *PhiR, double *PhiS, double *PhiT, double &OneMinusR, double &OneMinusS, double &OneMinusT, double &OnePlusR, double &OnePlusS, double &OnePlusT)
{
	int PerFace = polynomialDegree - 1;
	int PerFaceSquare = PerFace * PerFace;
	
	double *F_N = new double[6 * PerFaceSquare];
	int k = 0;
	for (int p1 = 0; p1 < PerFace; p1++)
	{
		for (int p2 = 0; p2 < PerFace; p2++)
		{
			F_N[k] = 0.5 * OneMinusT*PhiR[p1] * PhiS[p2];							// Face 1
			F_N[PerFaceSquare + k] = 0.5 * OneMinusS*PhiR[p1] * PhiT[p2];			// Face 2
			F_N[2 * PerFaceSquare + k] = 0.5 *OnePlusR*PhiS[p1] * PhiT[p2];			// Face 3
			F_N[3 * PerFaceSquare + k] = 0.5 *OnePlusS*PhiR[p1] * PhiT[p2];			// Face 4
			F_N[4 * PerFaceSquare + k] = 0.5 *OneMinusR*PhiS[p1] * PhiT[p2];		// Face 5
			F_N[5 * PerFaceSquare + k] = 0.5 *OnePlusT*PhiR[p1] * PhiS[p2];			// Face 6
			k++;
		}
	}
	return F_N;
}

//----------------------------面结点的形函数偏导数--------------------------------
// 先是对方向 1 的偏导数，然后是对方向 2 的偏导数，最后是对方向 3 的偏导数
double* Hexahedral::calcFace_dN(double *Coords, double *PhiR, double *PhiS, double *PhiT, double *dPhiR, double *dPhiS, double *dPhiT, double &OneMinusR, double &OneMinusS, double &OneMinusT, double &OnePlusR, double &OnePlusS, double &OnePlusT)
{
	int PerFace = polynomialDegree - 1;
	int PerFaceSquare = PerFace * PerFace;

	double *F_dN = new double[6 * PerFaceSquare * 3];

	int k = 0;
	for (int p1 = 0; p1 < PerFace; p1++)
	{
		for (int p2 = 0; p2 < PerFace; p2++)
		{
			F_dN[k] = 0.5 * OneMinusT*dPhiR[p1] * PhiS[p2];							//Face 1 wrt r
			F_dN[PerFaceSquare + k] = 0.5 * OneMinusS*dPhiR[p1] * PhiT[p2];			//Face 2 wrt r
			F_dN[2 * PerFaceSquare + k] = 0.5 *PhiS[p1] * PhiT[p2];					//Face 3 wrt r
			F_dN[3 * PerFaceSquare + k] = 0.5 *OnePlusS*dPhiR[p1] * PhiT[p2];		//Face 4 wrt r
			F_dN[4 * PerFaceSquare + k] = -0.5 *PhiS[p1] * PhiT[p2];				//Face 5 wrt r
			F_dN[5 * PerFaceSquare + k] = 0.5 *OnePlusT*dPhiR[p1] * PhiS[p2];		//Face 6 wrt r

			F_dN[6 * PerFaceSquare + k] = 0.5 * OneMinusT*PhiR[p1] * dPhiS[p2];		//Face 1 wrt s
			F_dN[7 * PerFaceSquare + k] = -0.5 *PhiR[p1] * PhiT[p2];				//Face 2 wrt s
			F_dN[8 * PerFaceSquare + k] = 0.5 *OnePlusR*dPhiS[p1] * PhiT[p2];		//Face 3 wrt s
			F_dN[9 * PerFaceSquare + k] = 0.5 *PhiR[p1] * PhiT[p2];					//Face 4 wrt s
			F_dN[10 * PerFaceSquare + k] = 0.5 *OneMinusR*dPhiS[p1] * PhiT[p2];		//Face 5 wrt s
			F_dN[11 * PerFaceSquare + k] = 0.5 *OnePlusT*PhiR[p1] * dPhiS[p2];		//Face 6 wrt s

			F_dN[12 * PerFaceSquare + k] = -0.5 *PhiR[p1] * PhiS[p2];				//Face 1 wrt t
			F_dN[13 * PerFaceSquare + k] = 0.5 * OneMinusS*PhiR[p1] * dPhiT[p2];	//Face 2 wrt t
			F_dN[14 * PerFaceSquare + k] = 0.5 *OnePlusR*PhiS[p1] * dPhiT[p2];		//Face 3 wrt t
			F_dN[15 * PerFaceSquare + k] = 0.5 *OnePlusS*PhiR[p1] * dPhiT[p2];		//Face 4 wrt t
			F_dN[16 * PerFaceSquare + k] = 0.5 *OneMinusR*PhiS[p1] * dPhiT[p2];		//Face 5 wrt t
			F_dN[17 * PerFaceSquare + k] = 0.5 *PhiR[p1] * PhiS[p2];				//Face 6 wrt t

			k++;
		}
	}
	return F_dN;
}

//----------------------------体内结点的形函数--------------------------------
double* Hexahedral::calcBulk_N(double *Coords, double *PhiR, double *PhiS, double *PhiT, double &OneMinusR, double &OneMinusS, double &OneMinusT, double &OnePlusR, double &OnePlusS, double &OnePlusT)
{
	int PerBulk = polynomialDegree - 1;
	int PerBulkCube = PerBulk * PerBulk * PerBulk;

	double *B_N = new double[PerBulkCube];

	int k = 0;	// k just helps to have the value in the correct place
	for (int p1 = 0; p1 < PerBulk; p1++)
	{
		for (int p2 = 0; p2 < PerBulk; p2++)
		{
			for (int p3 = 0; p3 < PerBulk; p3++)
			{
				B_N[k] = PhiR[p1] * PhiS[p2] * PhiT[p3];
				k++;
			}
		}
	}
	return B_N;
}

//----------------------------体内结点的形函数偏导数--------------------------------
// 先是对方向 1 的偏导数，然后是对方向 2 的偏导数，最后是对方向 3 的偏导数
double* Hexahedral::calcBulk_dN(double *Coords, double *PhiR, double *PhiS, double *PhiT, double *dPhiR, double *dPhiS, double *dPhiT, double &OneMinusR, double &OneMinusS, double &OneMinusT, double &OnePlusR, double &OnePlusS, double &OnePlusT)
{
	int PerBulk = polynomialDegree - 1;
	int PerBulkCube = PerBulk * PerBulk * PerBulk;

	double *B_dN = new double[PerBulkCube * 3];

	int k = 0;	// k just helps to have the value in the correct place
	for (int p1 = 0; p1 < PerBulk; p1++)
	{
		for (int p2 = 0; p2 < PerBulk; p2++)
		{
			for (int p3 = 0; p3 < PerBulk; p3++)
			{
				B_dN[k] = dPhiR[p1] * PhiS[p2] * PhiT[p3];						// Bulk wrt r
				B_dN[PerBulkCube + k] = PhiR[p1] * dPhiS[p2] * PhiT[p3];		// Bulk wrt s
				B_dN[2 * PerBulkCube + k] = PhiR[p1] * PhiS[p2] * dPhiT[p3];	// Bulk wrt t
				k++;
			}
		}
	}
	return B_dN;
}

//----------------------------单元形函数--------------------------------
double* Hexahedral::calcShapeN(double *Coords)
{	
	int PerField = polynomialDegree - 1;
	double* phiR = new double[PerField]();
	double* phiS = new double[PerField]();
	double* phiT = new double[PerField]();

	// 计算高阶多项式
	for (int p = 0; p < PerField; p++)
	{
		phiR[p] = calcPhi(p + 2, Coords[0]);
		phiS[p] = calcPhi(p + 2, Coords[1]);
		phiT[p] = calcPhi(p + 2, Coords[2]);
	}

	double oneMinusR = (1 - Coords[0]);
	double onePlusR = (1 + Coords[0]);

	double oneMinusS = (1 - Coords[1]);
	double onePlusS = (1 + Coords[1]);

	double oneMinusT = (1 - Coords[2]);
	double onePlusT = (1 + Coords[2]);

	// 计算每种模式下的形函数
	double* vertexN = calcVertex_N(Coords, oneMinusR, oneMinusS, oneMinusT, onePlusR, onePlusS, onePlusT);
	double* edgeN = calcEdge_N(Coords, phiR, phiS, phiT, oneMinusR, oneMinusS, oneMinusT, onePlusR, onePlusS, onePlusT);
	double* faceN = calcFace_N(Coords, phiR, phiS, phiT, oneMinusR, oneMinusS, oneMinusT, onePlusR, onePlusS, onePlusT);
	double* bulkN = calcBulk_N(Coords, phiR, phiS, phiT, oneMinusR, oneMinusS, oneMinusT, onePlusR, onePlusS, onePlusT);

	// 每个胞元节点数 NCNodes
	// 每个胞元自由度数 NCDof
	int NCNodes = (polynomialDegree + 1)*(polynomialDegree + 1)*(polynomialDegree + 1);
	int NCDof = NCNodes * 3;
	double* shapeN = new double[NCDof*3]();			//后面加括号也表示调用默认构造函数，初始化为0
	//memset(shapeN, 0, NCDof*3 * sizeof(double));	//赋初值0


	// 计算形函数
	/*--------------------------------------------------------------------------------------
			|  N1 0 0  N2 0 0  N3 0 0  N4 0 0  ...  edge N  ...  face N  ...  bulk N  |
		N =	|  0 N1 0  0 N2 0  0 N3 0  0 N4 0  ...  edge N  ...  face N  ...  bulk N  |
			|  0 0 N1  0 0 N2  0 0 N3  0 0 N4  ...  edge N  ...  face N  ...  bulk N  |
	--------------------------------------------------------------------------------------*/
	int n = 0;
	int m = 0;

	// 角结点形函数
	for (n = 0; n < 8; n++)
	{
		m = 3 * n;
		shapeN[m] = vertexN[n];						// 矩阵 N 的第 1 行
		shapeN[NCDof + m + 1] = vertexN[n];			// 矩阵 N 的第 2 行
		shapeN[2 * NCDof + m + 2] = vertexN[n];		// 矩阵 N 的第 3 行
	}

	// 边结点形函数
	for (n = 0; n < 12* PerField; n++)
	{
		m = 3 * (n + 8);							
		shapeN[m] = edgeN[n];					// 矩阵 N 的第 1 行
		shapeN[NCDof + m + 1] = edgeN[n];		// 矩阵 N 的第 2 行
		shapeN[2 *NCDof + m + 2] = edgeN[n];	// 矩阵 N 的第 3 行
	}

	// 面结点形函数
	for (n = 0; n < 6 * PerField * PerField; n++)
	{
		m =3 * (n + 8 + 12 * PerField);
		shapeN[m] = faceN[n];					// 矩阵 N 的第 1 行
		shapeN[NCDof + m + 1] = faceN[n];		// 矩阵 N 的第 2 行
		shapeN[2 * NCDof + m + 2] = faceN[n];	// 矩阵 N 的第 3 行
	}

	// 体内结点形函数
	for (n = 0; n < PerField * PerField* PerField; n++)
	{
		m = 3 * (n + 8 + 12 * PerField + 6 * PerField * PerField);
		shapeN[m] = bulkN[n];					// 矩阵 N 的第 1 行
		shapeN[NCDof + m + 1] = bulkN[n];		// 矩阵 N 的第 2 行
		shapeN[2 * NCDof + m + 2] = bulkN[n];	// 矩阵 N 的第 3 行
	}

	delete[]phiR; phiR = nullptr;
	delete[]phiS; phiS = nullptr;
	delete[]phiS; phiT = nullptr;
	delete[]vertexN;	vertexN = nullptr;
	delete[]edgeN;	edgeN = nullptr;
	delete[]faceN;	faceN = nullptr;
	delete[]bulkN;  bulkN = nullptr;

	return shapeN;
}

//----------------------------输出单元形函数--------------------------------
void Hexahedral::printShapeN(double *Coords)
{
	double* N = calcShapeN(Coords);
	int perField = polynomialDegree - 1;
	int NCNodes = (polynomialDegree + 1)*(polynomialDegree + 1)*(polynomialDegree + 1);
	int NCDof = NCNodes * 3;
	printf("\nVertex Shape N:\n");
	for (int i = 0; i < NCNodes; i++)
	{
		printf("%-6.4f, ", N[3 * i]);
		if (i == 7)		printf("\n");
		if (i == (7 + 12 * perField))	printf("\n");
		if (i == (7 + 12 * perField + 6 * perField * perField))		printf("\n");
	}
	printf("\n");
}




//----------------------------单元形函数偏导数--------------------------------
double* Hexahedral::calcShape_dN(double *Coords)
{
	int PerField = polynomialDegree - 1;
	double* phiR = new double[PerField]();
	double* phiS = new double[PerField]();
	double* phiT = new double[PerField]();

	double* dphiR = new double[PerField]();
	double* dphiS = new double[PerField]();
	double* dphiT = new double[PerField]();

	// 计算高阶多项式
	for (int p = 0; p < PerField; p++)
	{
		phiR[p] = calcPhi(p + 2, Coords[0]);
		phiS[p] = calcPhi(p + 2, Coords[1]);
		phiT[p] = calcPhi(p + 2, Coords[2]);
		dphiR[p] = calcDerivOfPhi(p + 2, Coords[0]);
		dphiS[p] = calcDerivOfPhi(p + 2, Coords[1]);
		dphiT[p] = calcDerivOfPhi(p + 2, Coords[2]);
	}

	double oneMinusR = (1 - Coords[0]);
	double onePlusR = (1 + Coords[0]);

	double oneMinusS = (1 - Coords[1]);
	double onePlusS = (1 + Coords[1]);

	double oneMinusT = (1 - Coords[2]);
	double onePlusT = (1 + Coords[2]);

	// 计算每种模式下的形函数导数
	double* vertexdN = calcVertex_dN(Coords, oneMinusR, oneMinusS, oneMinusT, onePlusR, onePlusS, onePlusT);
	double* edgedN = calcEdge_dN(Coords, phiR, phiS, phiT, dphiR, dphiS, dphiT, oneMinusR, oneMinusS, oneMinusT, onePlusR, onePlusS, onePlusT);
	double* facedN = calcFace_dN(Coords, phiR, phiS, phiT, dphiR, dphiS, dphiT, oneMinusR, oneMinusS, oneMinusT, onePlusR, onePlusS, onePlusT);
	double* bulkdN = calcBulk_dN(Coords, phiR, phiS, phiT, dphiR, dphiS, dphiT, oneMinusR, oneMinusS, oneMinusT, onePlusR, onePlusS, onePlusT);

	// 每个胞元节点数 NCNodes
	// 每个胞元自由度数 NCDof
	int NCNodes = (polynomialDegree + 1)*(polynomialDegree + 1)*(polynomialDegree + 1);
	int NCDof = NCNodes * 3;
	double* dN = new double[NCNodes * 3]();
	// 计算形函数
/*--------------------------------------------------------------------------------------
		|  dN1_dx  dN2_x  dN3_x  ...  edge dN_dx  ...  face dN_dx  ...  bulk dN_dx  |
  dN =	|  dN1_dy  dN2_y  dN3_y  ...  edge dN_dy  ...  face dN_dy  ...  bulk dN_dy  |
		|  dN1_dz  dN2_z  dN3_z  ...  edge dN_dz  ...  face dN_dz  ...  bulk dN_dz  |
--------------------------------------------------------------------------------------*/
	int n = 0;
	int m = 0;

	// 角结点
	for (n = 0; n < 8; n++)
	{
		dN[n] = vertexdN[n];
		dN[NCNodes + n] = vertexdN[8 + n];
		dN[2 * NCNodes + n] = vertexdN[16 + n];
	}

	// 边
	for (n = 0; n < 12* PerField; n++)
	{
		dN[n + 8] = edgedN[n];
		dN[NCNodes + n + 8] = edgedN[12 * PerField + n];
		dN[2 * NCNodes + n + 8] = edgedN[24 * PerField + n];
	}

	// 面
	int FaceSquare = PerField * PerField;
	for (n = 0; n < 6 * FaceSquare; n++)
	{
		m = 8 + 12 * PerField;
		dN[n + m] = facedN[n];
		dN[NCNodes + n + m] = facedN[6 * FaceSquare + n];
		dN[2 * NCNodes + n + m] = facedN[12 * FaceSquare + n];
	}

	// 体
	int BulkCube = PerField * PerField * PerField;
	for (n = 0; n < BulkCube; n++)
	{
		m = 8 + 12 * PerField + 6 * FaceSquare;
		dN[n + m] = bulkdN[n];
		dN[NCNodes + n + m] = bulkdN[1 * BulkCube + n];
		dN[2 * NCNodes + n + m] = bulkdN[2 * BulkCube + n];
	}

	delete[]phiR; phiR = nullptr;
	delete[]phiS; phiS = nullptr;
	delete[]phiS; phiT = nullptr;
	delete[]dphiR; dphiR = nullptr;
	delete[]dphiS; dphiS = nullptr;
	delete[]dphiS; dphiT = nullptr;
	delete[]vertexdN;	vertexdN = nullptr;
	delete[]edgedN;	 edgedN = nullptr;
	delete[]facedN;	 facedN = nullptr;
	delete[]bulkdN;  bulkdN = nullptr;

	return dN;
}

void Hexahedral::printShape_dN(double *Coords)
{
	double* dN = calcShape_dN(Coords);
	int perField = polynomialDegree - 1;
	int NCNodes = (polynomialDegree + 1)*(polynomialDegree + 1)*(polynomialDegree + 1);
	int NCDof = NCNodes * 3;
	printf("\nVertex Shape dN:\n");
	for (int i = 0; i < NCNodes; i++)
	{
		for (int j = 0; j < 3; j++)
		{
			printf("%-6.4f, ", dN[j*NCNodes + i]);
		}
		printf("\n");
	}
	printf("\n");
}

double* Hexahedral::calcMatrixB(double *Coords)
{
	// 每个胞元节点数 NCNodes
	// 每个胞元自由度数 NCDof
	int NCNodes = (polynomialDegree + 1)*(polynomialDegree + 1)*(polynomialDegree + 1);
	int NCDof = NCNodes * 3;
	double* Bc = new double[6 * NCDof]();			//后面加括号也表示调用默认构造函数，初始化为0

	double* dN = calcShape_dN(Coords);
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

	double *jac = calcJacobiAtCenter();		//Jacobi矩阵

	int n = 0;
	int m = 0;

	for (int n = 0; n < NCNodes; n++)
	{
		m = 3 * n;
		Bc[m] = dN[n] / jac[0];
		Bc[NCDof + m + 1] = dN[NCNodes + n] / jac[4];
		Bc[2 * NCDof + m + 2] = dN[ 2 * NCNodes + n] / jac[8];

		Bc[3 * NCDof + m] = dN[NCNodes + n] / jac[4];
		Bc[3 * NCDof + m + 1] = dN[n] / jac[0];

		Bc[4 * NCDof + m + 1] = dN[2 * NCNodes + n] / jac[8];
		Bc[4 * NCDof + m + 2] = dN[NCNodes + n] / jac[4];

		Bc[5 * NCDof + m] = dN[2 * NCNodes + n] / jac[8];
		Bc[5 * NCDof + m + 2] = dN[n] / jac[0];
	}
	delete[]jac;  jac = nullptr;
	delete[]dN;  dN = nullptr;

	return Bc;
}


//----------------------------单元 应变-位移矩阵 B --------------------------------
//double* Hexahedral::calcMatrixB(double *Coords)
//{
//	int PerField = polynomialDegree - 1;
//	double* phiR = new double[PerField]();
//	double* phiS = new double[PerField]();
//	double* phiT = new double[PerField]();
//
//	double* dphiR = new double[PerField]();
//	double* dphiS = new double[PerField]();
//	double* dphiT = new double[PerField]();
//
//	// 计算高阶多项式
//	for (int p = 0; p < PerField; p++)
//	{
//		phiR[p] = calcPhi(p + 2, Coords[0]);
//		phiS[p] = calcPhi(p + 2, Coords[1]);
//		phiT[p] = calcPhi(p + 2, Coords[2]);
//		dphiR[p] = calcDerivOfPhi(p + 2, Coords[0]);
//		dphiS[p] = calcDerivOfPhi(p + 2, Coords[1]);
//		dphiT[p] = calcDerivOfPhi(p + 2, Coords[2]);
//	}
//
//	double oneMinusR = (1 - Coords[0]);
//	double onePlusR = (1 + Coords[0]);
//
//	double oneMinusS = (1 - Coords[1]);
//	double onePlusS = (1 + Coords[1]);
//
//	double oneMinusT = (1 - Coords[2]);
//	double onePlusT = (1 + Coords[2]);
//
//	// 计算每种模式下的形函数
//	double* vertexdN = calcVertex_dN(Coords, oneMinusR, oneMinusS, oneMinusT, onePlusR, onePlusS, onePlusT);
//	double* edgedN = calcEdge_dN(Coords, phiR, phiS, phiT, dphiR, dphiS, dphiT, oneMinusR, oneMinusS, oneMinusT, onePlusR, onePlusS, onePlusT);
//	double* facedN = calcFace_dN(Coords, phiR, phiS, phiT, dphiR, dphiS, dphiT, oneMinusR, oneMinusS, oneMinusT, onePlusR, onePlusS, onePlusT);
//	double* bulkdN = calcBulk_dN(Coords, phiR, phiS, phiT, dphiR, dphiS, dphiT, oneMinusR, oneMinusS, oneMinusT, onePlusR, onePlusS, onePlusT);
//
//	// 每个胞元节点数 NCNodes
//	// 每个胞元自由度数 NCDof
//	int NCNodes = (polynomialDegree + 1)*(polynomialDegree + 1)*(polynomialDegree + 1);
//	int NCDof = NCNodes * 3;
//	double* Bc = new double[6 * NCDof]();			//后面加括号也表示调用默认构造函数，初始化为0
//
//	// 计算形函数
//	/*-----------------------------------------------------------------------------------------------------------------
//			|  dN1_dx     0       0	    dN2_dx     0       0	...  edge dN  ...  face dN  ...  bulk dN  |
//			|     0    dN1_dy     0	       0    dN2_dy     0	...  edge dN  ...  face dN  ...  bulk dN  |
//		B =	|     0       0    dN1_dz	   0       0    dN2_dz  ...  edge dN  ...  face dN  ...  bulk dN  |
//			|  dN1_dy  dN1_dx     0     dN2_dy  dN2_dx     0    ...  edge dN  ...  face dN  ...  bulk dN  |
//			|     0    dN1_dz  dN1_dy      0    dN2_dz  dN2_dy  ...  edge dN  ...  face dN  ...  bulk dN  |
//			|  dN1_dz    0     dN1_dx   dN2_dz    0     dN1_dx  ...  edge dN  ...  face dN  ...  bulk dN  |
//	-----------------------------------------------------------------------------------------------------------------*/
//
//	double *jac = calcJacobiAtCenter();		//Jacobi矩阵
//
//	int n = 0;
//	int m = 0;
//	// 角结点形函数
//	for (n = 0; n < 8; n++)
//	{
//		m = 3 * n;
//		Bc[m] = vertexdN[n] / jac[0];
//		Bc[NCDof + m + 1] = vertexdN[8 + n] / jac[4];
//		Bc[2 * NCDof + m + 2] = vertexdN[16 + n] / jac[8];
//
//		Bc[3 * NCDof + m] = vertexdN[8 + n] / jac[4];
//		Bc[3 * NCDof + m + 1] = vertexdN[n] / jac[0];
//
//		Bc[4 * NCDof + m + 1] = vertexdN[16 + n] / jac[8];
//		Bc[4 * NCDof + m + 2] = vertexdN[8 + n] / jac[4];
//
//		Bc[5 * NCDof + m] = vertexdN[16 + n] / jac[8];
//		Bc[5 * NCDof + m + 2] = vertexdN[n] / jac[0];		
//	}
//
//	// 边结点形函数
//	for (n = 0; n < 12 * PerField; n++)
//	{
//		m = 3 * (n + 8);
//		Bc[m] = edgedN[n] / jac[0];
//		Bc[NCDof + m + 1] = edgedN[12 + n] / jac[4];
//		Bc[2 * NCDof + m + 2] = edgedN[24 + n] / jac[8];
//
//		Bc[3 * NCDof + m] = edgedN[12 + n] / jac[4];
//		Bc[3 * NCDof + m + 1] = edgedN[n] / jac[0];
//
//		Bc[4 * NCDof + m + 1] = edgedN[24 + n] / jac[8];
//		Bc[4 * NCDof + m + 2] = edgedN[12 + n] / jac[4];
//
//		Bc[5 * NCDof + m] = edgedN[24 + n] / jac[8];
//		Bc[5 * NCDof + m + 2] = edgedN[n] / jac[0];
//	}
//
//	// 面结点形函数
//	for (n = 0; n < 6 * PerField * PerField; n++)
//	{
//		m = 3 * (n + 8 + 12 * PerField);
//		Bc[m] = facedN[n] / jac[0];
//		Bc[NCDof + m + 1] = facedN[6 + n] / jac[4];
//		Bc[2 * NCDof + m + 2] = facedN[12 + n] / jac[8];
//
//		Bc[3 * NCDof + m] = facedN[6 + n] / jac[4];
//		Bc[3 * NCDof + m + 1] = facedN[n] / jac[0];
//
//		Bc[4 * NCDof + m + 1] = facedN[12 + n] / jac[8];
//		Bc[4 * NCDof + m + 2] = facedN[6 + n] / jac[4];
//
//		Bc[5 * NCDof + m] = facedN[12 + n] / jac[8];
//		Bc[5 * NCDof + m + 2] = facedN[n] / jac[0];
//	}
//
//	// 体内结点形函数
//	for (n = 0; n < PerField * PerField* PerField; n++)
//	{
//		m = 3 * (n + 8 + 12 * PerField + 6 * PerField * PerField);
//		Bc[m] = bulkdN[n] / jac[0];
//		Bc[NCDof + m + 1] = bulkdN[1 + n] / jac[4];
//		Bc[2 * NCDof + m + 2] = bulkdN[2 + n] / jac[8];
//
//		Bc[3 * NCDof + m] = bulkdN[1 + n] / jac[4];
//		Bc[3 * NCDof + m + 1] = bulkdN[n] / jac[0];
//
//		Bc[4 * NCDof + m + 1] = bulkdN[2 + n] / jac[8];
//		Bc[4 * NCDof + m + 2] = bulkdN[1 + n] / jac[4];
//
//		Bc[5 * NCDof + m] = bulkdN[2 + n] / jac[8];
//		Bc[5 * NCDof + m + 2] = bulkdN[n] / jac[0];
//	}
//
//	delete[]phiR; phiR = nullptr;
//	delete[]phiS; phiS = nullptr;
//	delete[]phiS; phiT = nullptr;
//	delete[]dphiR; dphiR = nullptr;
//	delete[]dphiS; dphiS = nullptr;
//	delete[]dphiS; dphiT = nullptr;
//	delete[]vertexdN;	vertexdN = nullptr;
//	delete[]edgedN;	 edgedN = nullptr;
//	delete[]facedN;	 facedN = nullptr;
//	delete[]bulkdN;  bulkdN = nullptr;
//	delete[]jac;  jac = nullptr;
//
//	return Bc;
//}

void Hexahedral::printMatrixB(double *Coords)
{
	double* B = calcMatrixB(Coords);
	int perField = polynomialDegree - 1;
	int NCNodes = (polynomialDegree + 1)*(polynomialDegree + 1)*(polynomialDegree + 1);
	int NCDof = NCNodes * 3;
	printf("\nElement B:-----------------------\n");
	//for (int i = 0; i < NCDof *6; i++)
	//{
	//	if (i % NCDof == 0)
	//	{
	//		printf("\n");
	//	}
	//	printf("%-6.4f, ", B[i]);

	//}

	int index = 0;
	for (int j = 0; j < NCDof; j++)
	{
		for (int i = 0; i < 6; i++)
		{
			index = i * NCDof + j;
			printf("%-6.4f, ", B[index]);
		}
		printf("\n");
	}

	printf("\n");
}




//----------------------------原始父胞元胞元的 Jacobian 矩阵--------------------------------
double* Hexahedral::calcJacobiAtCenter()
{
	// 简易版， 正规操作应该是把每个高斯点代入形函数矩阵再和结点坐标矩阵相乘
	double* Jac = new double[9]();			//Jacobian矩阵	

	double xc = nodes[1]->coords_[0] - nodes[0]->coords_[0];
	double yc = nodes[2]->coords_[1] - nodes[1]->coords_[1];
	double zc = nodes[4]->coords_[2] - nodes[0]->coords_[2];
	Jac[0] = 0.5 * xc;
	Jac[4] = 0.5 * yc;
	Jac[8] = 0.5 * zc;
	return Jac;
}

//----------------------------Jacobian 矩阵行列式--------------------------------
double Hexahedral::calcDetJac(double *VC)
{
	double xc = VC[3] - VC[0];
	double yc = VC[7] - VC[4];
	double zc = VC[14] - VC[2];
	double detJ = 0.5 * xc * 0.5 * yc *  0.5 * zc;
	return detJ;
}

//----------------------------子胞元的 Jacobian 矩阵行列式--------------------------------
double Hexahedral::calcSubCellDetJac(const int &SubDepth)
{
	double Jac11 = 1.0 / pow(2.0, SubDepth);
	double subDetJ = Jac11 * Jac11 * Jac11;
	return subDetJ;
}




//----------------------------计算胞元的 刚度矩阵(包括子胞元)  对每个单元遍历循环的总函数--------------------------------
double* Hexahedral::calcStiffnessMatrix(EmbeddedDomain *pEmDo, AbsMaterial *pMaterial, GaussIntegrator *pIntegrator, const int Depth)
{

	double* vertexCoords = getCellVertexCoords();		// 胞元顶点的坐标

	int NCDof = (polynomialDegree + 1)*(polynomialDegree + 1)*(polynomialDegree + 1) * 3;
	double* Kc = new double[NCDof * NCDof]();			// 胞元单刚

	double cellDetJ = calcDetJac(vertexCoords);			// 胞元Jacobi矩阵行列式

	if (isInside(pEmDo, vertexCoords))
	{
		printf("\nElement %d is in physic domain\n", this->id);
		calcUndividedStiffness(pMaterial, pIntegrator, Kc, cellDetJ);
	}
	else if (isIntersected(pEmDo, vertexCoords))
	{
		printf("\nElement %d is intersected\n", this->id);
		int depthCounter = 0;
		double LocalVerCoo[24] = { -1.0,-1.0,-1.0, 1.0,-1.0,-1.0, 1.0, 1.0,-1.0, -1.0, 1.0,-1.0, -1.0,-1.0,1.0, 1.0,-1.0,1.0, 1.0, 1.0,1.0, -1.0, 1.0,1.0, };
		partitionRecursively(pEmDo, pMaterial, pIntegrator, vertexCoords, LocalVerCoo, depthCounter, Depth, Kc, cellDetJ);
	}
	else
	{
		printf("\nElement %d is outside the physic domain\n", this->id);
	}

	delete[]vertexCoords; vertexCoords = nullptr;
	return Kc;
}



//----------------------------八叉树细分计算 子胞元的 刚度矩阵--------------------------------
void Hexahedral::partitionRecursively(EmbeddedDomain *pEmDo, AbsMaterial *pMaterial, GaussIntegrator *pIntegrator, double *GlobalVerCoo, double *LocalVerCoo,
	int depthCounter, const int Depth, double *Kc, double CellDetJ)
{
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

	if (isIntersected(pEmDo, GlobalVerCoo) && depthCounter < Depth)
	{
		depthCounter += 1;		// 深度 +1

		double** GlobalCoords = new double*[27];		// 子胞元的 全局坐标
		double** LocalCoords = new double*[27];			// 子胞元相对于父胞元的 局部坐标

		double SplitLocalCoords[27][3] = {
			{-1.0,-1.0,-1.0}, {0.0,-1.0,-1.0}, {1.0,-1.0,-1.0}, {-1.0, 0.0,-1.0}, {0.0, 0.0,-1.0}, {1.0, 0.0,-1.0}, {-1.0, 1.0,-1.0}, {0.0, 1.0,-1.0}, {1.0, 1.0,-1.0},
			{-1.0,-1.0, 0.0}, {0.0,-1.0, 0.0}, {1.0,-1.0, 0.0}, {-1.0, 0.0, 0.0}, {0.0, 0.0, 0.0}, {1.0, 0.0, 0.0}, {-1.0, 1.0, 0.0}, {0.0, 1.0, 0.0}, {1.0, 1.0, 0.0},
			{-1.0,-1.0, 1.0}, {0.0,-1.0, 1.0}, {1.0,-1.0, 1.0}, {-1.0, 0.0, 1.0}, {0.0, 0.0, 1.0}, {1.0, 0.0, 1.0}, {-1.0, 1.0, 1.0}, {0.0, 1.0, 1.0}, {1.0, 1.0, 1.0} };
		
		for (int i = 0; i < 27; i++)
		{
			GlobalCoords[i] = mapLocalToGlobal(SplitLocalCoords[i], GlobalVerCoo);		// 将 子胞元 顶点坐标映射到 父胞元 中的 全局坐标
			LocalCoords[i] = mapLocalToGlobal(SplitLocalCoords[i], LocalVerCoo);		// 将 子胞元 顶点坐标映射到 父胞元 中的 局部坐标
		}

		int num1[8] = { 0, 1, 4, 3, 9, 10, 13, 12 };
		int num2[8] = { 1, 2, 5, 4, 10, 11, 14, 13 };
		int num3[8] = { 4, 5, 8, 7, 13, 14, 17, 16 };
		int num4[8] = { 3, 4, 7, 6, 12, 13, 16, 15 };
		int num5[8] = { 9, 10, 13, 12, 18, 19, 22, 21 };
		int num6[8] = { 10, 11, 14, 13, 19, 20, 23, 22 };
		int num7[8] = { 13, 14, 17, 16, 22, 23, 26, 25 };
		int num8[8] = { 12, 13, 16, 15, 21, 22, 25, 24 };

		double globalSubCell1[24], globalSubCell2[24], globalSubCell3[24], globalSubCell4[24],
			   globalSubCell5[24], globalSubCell6[24], globalSubCell7[24], globalSubCell8[24];

		double localSubCell1[24], localSubCell2[24], localSubCell3[24], localSubCell4[24],
			localSubCell5[24], localSubCell6[24], localSubCell7[24], localSubCell8[24];

		int m = 0;
		for (int i = 0; i < 8; i++)
		{
			for (int j = 0; j < 3; j++)
			{
				m = 3 * i + j;
				globalSubCell1[m] = GlobalCoords[num1[i]][j];
				globalSubCell2[m] = GlobalCoords[num2[i]][j];
				globalSubCell3[m] = GlobalCoords[num3[i]][j];
				globalSubCell4[m] = GlobalCoords[num4[i]][j];
				globalSubCell5[m] = GlobalCoords[num5[i]][j];
				globalSubCell6[m] = GlobalCoords[num6[i]][j];
				globalSubCell7[m] = GlobalCoords[num7[i]][j];
				globalSubCell8[m] = GlobalCoords[num8[i]][j];

				localSubCell1[m] = LocalCoords[num1[i]][j];
				localSubCell2[m] = LocalCoords[num2[i]][j];
				localSubCell3[m] = LocalCoords[num3[i]][j];
				localSubCell4[m] = LocalCoords[num4[i]][j];
				localSubCell5[m] = LocalCoords[num5[i]][j];
				localSubCell6[m] = LocalCoords[num6[i]][j];
				localSubCell7[m] = LocalCoords[num7[i]][j];
				localSubCell8[m] = LocalCoords[num8[i]][j];
			}
		}

		for (int i = 0; i < 27; i++)	delete[]GlobalCoords[i];
		delete[]GlobalCoords; 		GlobalCoords = nullptr;

		for (int i = 0; i < 27; i++)	delete[]LocalCoords[i];
		delete[]LocalCoords; 		LocalCoords = nullptr;

		partitionRecursively(pEmDo, pMaterial, pIntegrator, globalSubCell1, localSubCell1, depthCounter, Depth, Kc, CellDetJ);
		partitionRecursively(pEmDo, pMaterial, pIntegrator, globalSubCell2, localSubCell2, depthCounter, Depth, Kc, CellDetJ);
		partitionRecursively(pEmDo, pMaterial, pIntegrator, globalSubCell3, localSubCell3, depthCounter, Depth, Kc, CellDetJ);
		partitionRecursively(pEmDo, pMaterial, pIntegrator, globalSubCell4, localSubCell4, depthCounter, Depth, Kc, CellDetJ);

		partitionRecursively(pEmDo, pMaterial, pIntegrator, globalSubCell5, localSubCell5, depthCounter, Depth, Kc, CellDetJ);
		partitionRecursively(pEmDo, pMaterial, pIntegrator, globalSubCell6, localSubCell6, depthCounter, Depth, Kc, CellDetJ);
		partitionRecursively(pEmDo, pMaterial, pIntegrator, globalSubCell7, localSubCell7, depthCounter, Depth, Kc, CellDetJ);
		partitionRecursively(pEmDo, pMaterial, pIntegrator, globalSubCell8, localSubCell8, depthCounter, Depth, Kc, CellDetJ);
	}

	else
	{
		//printf("\n-----After ifififififif depth Counter = %d-----\n", depthCounter);
		//printf("\nSub Cell LocalCoords:\n");
		//for (int ii = 0; ii < 8; ii++)
		//{
		//	for (int jj = 0; jj < 3; jj++)
		//	{
		//		printf("%-5.3f, ", LocalVerCoo[3 * ii + jj]);
		//	}
		//	printf("\n");			
		//}

		//printf("\nSub Cell GlobalCoords:\n");
		//for (int ii = 0; ii < 8; ii++)
		//{
		//	for (int jj = 0; jj < 3; jj++)
		//	{
		//		printf("%-5.3f, ", GlobalVerCoo[3 * ii + jj]);
		//	}
		//	printf("\n");
		//}

		CalcSubCellStifiness(pEmDo, pMaterial, pIntegrator, GlobalVerCoo, LocalVerCoo, depthCounter, Kc, CellDetJ);
	}

}


//----------------------------计算 子胞元的 刚度矩阵--------------------------------
void Hexahedral::CalcSubCellStifiness(EmbeddedDomain *pEmDo, AbsMaterial *pMaterial, GaussIntegrator *pIntegrator, double *GlobalVerCoo, double *LocalVerCoo,
	const int depthCounter, double *Kc, const double CellDetJ)
{

	double* GaussPoints = pIntegrator->getGaussCoordinates();	// 这个指针不用删,留给类的析构函数释放
	double* GaussWeights = pIntegrator->getGaussWeights();		// 这个指针不用删,留给类的析构函数释放
	double* MatD = pMaterial->getMaterialMatrix();				// 这个指针不用删,留给类的析构函数释放
	

	double subDetJ = calcSubCellDetJac(depthCounter);			// 子胞元的的Jacobi矩阵行列式
	double scalingFactor;										// 定义缩放因子，在物理域内就 =1，在物理域外就 =alpha

	int NCDof = (polynomialDegree + 1)*(polynomialDegree + 1)*(polynomialDegree + 1) * 3;
	int m = NCDof;
	int n =6, k = 6;		//应变位移矩阵的大小
	double alpha1 = 1.0, beta1 = 0.0, alpha2, beta2 = 1.0;
	int lda = NCDof, ldb = 6, ldc = 6;


	int nip = pIntegrator->NIP;					// 单方向积分点个数
	double* Bc;									// 应变位移矩阵，不用分配大小，函数 “calcMatrixB()” 会分配大小，所以用完后要 delete[]Bc
	double* BTD = new double[NCDof * 6]();		// B的转置乘以D = BTD, 在所有积分点遍历后释放

	double subCellIntegratePoint[3];			// 积分点在 子胞元 中的坐标
	double* cellIntegratePoint;					// 积分点在 大胞元 中的坐标
	double* globalPoint;						// 积分点的 全局坐标，用来判断是否在物理域内

	for (int kk = 0; kk < nip; kk++)
	{
		for (int jj = 0; jj < nip; jj++)
		{
			for (int ii = 0; ii < nip; ii++)
			{
				// 积分点三个方向的坐标
				subCellIntegratePoint[0] = GaussPoints[ii];
				subCellIntegratePoint[1] = GaussPoints[jj];
				subCellIntegratePoint[2] = GaussPoints[kk];

				globalPoint = mapLocalToGlobal(subCellIntegratePoint, GlobalVerCoo);			// 将高斯点映射到全局坐标
				scalingFactor = pEmDo->getMaterialAtPoint(globalPoint);							// 判断该高斯积分点位置，并返回材料惩罚参数 1 或 Alpha
				cellIntegratePoint = mapLocalToGlobal(subCellIntegratePoint, LocalVerCoo);		// 将SubCell的高斯点映射到Cell中的局部坐标	

				Bc = calcMatrixB(cellIntegratePoint);

				// BTD = BT*D;
				cblas_dgemm(CblasRowMajor, CblasTrans, CblasNoTrans, m, n, k, alpha1, Bc, lda, MatD, ldb, beta1, BTD, ldc);

				// Kc = BT*D*B + Kc;
				alpha2 = CellDetJ * subDetJ * GaussWeights[ii] * GaussWeights[jj] * GaussWeights[kk] * scalingFactor;
				cblas_dgemm(CblasRowMajor, CblasNoTrans, CblasNoTrans, NCDof, NCDof, 6, alpha2, BTD, 6, Bc, NCDof, beta2, Kc, NCDof);

				delete[]Bc;
				delete[]cellIntegratePoint;
				delete[]globalPoint;
			}
		}
	}
	delete[]BTD;
}


//----------------------------计算完全在物理域内的胞元的 刚度矩阵--------------------------------
void Hexahedral::calcUndividedStiffness(AbsMaterial *pMaterial, GaussIntegrator *pIntegrator, double *Kc, double CellDetJ)
{
	
	double* GaussPoints = pIntegrator->getGaussCoordinates();	// 这个指针不用删,留给类的析构函数释放
	double* GaussWeights = pIntegrator->getGaussWeights();		// 这个指针不用删,留给类的析构函数释放
	double* MatD = pMaterial->getMaterialMatrix();				// 这个指针不用删,留给类的析构函数释放

	double integratePoint[3];

	int m, n, k;
	double alpha1 = 1.0, beta1 = 0.0, alpha2, beta2 = 1.0;
	int lda, ldb, ldc;

	//C = alpha * A*B + beta * C;
	// cblas_dgemm()--------参数说明
	// CblasRowMajor--------按行存储
	// CblasTrans-----------表示A的转置,后面用op表示,op = CblasTrans
	// CblasNoTrans---------表示B不转置,后面用op表示,op = blasNoTrans
	// m, n, k -------------op(A):[m*k]   op(B):[k*n]   C:[m*n]
	// alpha ---------------实数， alpha*A*B 
	// A -------------------矩阵A
	// lda -----------------矩阵A的列数
	// B -------------------矩阵B
	// ldb -----------------矩阵B的列数
	// beta ----------------实数， beta*C
	// C -------------------矩阵C
	// ldc -----------------矩阵C的列数

	int NCDof = (polynomialDegree + 1)*(polynomialDegree + 1)*(polynomialDegree + 1) * 3;
	m = NCDof;
	k = 6;	n = 6;		//应变位移矩阵的大小
	lda = NCDof;	ldb = 6;  ldc = 6;
	
	int nip = pIntegrator->NIP;				// 单方向积分点个数
	double* Bc;								// 应变位移矩阵，不用分配大小，函数 “calcMatrixB()” 会分配大小，所以用完后要 delete[]Bc
	double* BTD = new double[NCDof*6]();			// B的转置乘以D = BTD, 在所有积分点遍历后释放

	for (int kk = 0; kk < nip; kk++)
	{
		for (int jj = 0; jj < nip; jj++)
		{
			for (int ii = 0; ii < nip; ii++)
			{
				// 积分点三个方向的坐标
				integratePoint[0] = GaussPoints[ii];
				integratePoint[1] = GaussPoints[jj];
				integratePoint[2] = GaussPoints[kk];
				//printf("[ %-6.4f, %-6.4f, %-6.4f ]\n", integratePoint[0], integratePoint[1], integratePoint[2]);
				Bc = calcMatrixB(integratePoint);

				// BTD = BT*D;
				cblas_dgemm(CblasRowMajor, CblasTrans, CblasNoTrans, m, n, k, alpha1, Bc, lda, MatD, ldb, beta1, BTD, ldc);

				// Kc = BT*D*B + Kc;
				alpha2 = CellDetJ * GaussWeights[ii] * GaussWeights[jj] * GaussWeights[kk];
				cblas_dgemm(CblasRowMajor, CblasNoTrans, CblasNoTrans, NCDof, NCDof, 6, alpha2, BTD, 6, Bc, NCDof, beta2, Kc, NCDof);

				delete[]Bc;
			}
		}
	}
	delete[]BTD;
}


//----------------------------判断胞元是否在完全在 物理域外--------------------------------
bool Hexahedral::isOutside(EmbeddedDomain *pEmDo, double* VerCoo)
{
	const int numOfSeedPoints = 4;
	double StepSizeX = (VerCoo[3] - VerCoo[0]) / (numOfSeedPoints - 1);
	double StepSizeY = (VerCoo[7] - VerCoo[4]) / (numOfSeedPoints - 1);
	double StepSizeZ = (VerCoo[12] - VerCoo[0]) / (numOfSeedPoints - 1);

	double point[3];
	int currentDomainIndex = 1;
	bool check = true;
	int i, j, k;

	for (k = 0; k < numOfSeedPoints; k++)
	{
		for (j = 0; j < numOfSeedPoints; j++)
		{
			for (i = 0; i < numOfSeedPoints; i++)
			{
				point[0] = VerCoo[0] + i * StepSizeX;
				point[1] = VerCoo[1] + j * StepSizeY;
				point[2] = VerCoo[2] + k * StepSizeZ;

				currentDomainIndex = pEmDo->getDomainIndex(point);
				if (currentDomainIndex == 2)
				{
					check = false;		//一旦判断顶点的所在位置为 域内，就返回false
					goto BREAKLOOP1;
				}
			}
		}
	}
BREAKLOOP1:;
	return check;
}

//----------------------------判断胞元是否在完全在 物理域内--------------------------------
bool Hexahedral::isInside(EmbeddedDomain *pEmDo, double* VerCoo)
{
	const int numOfSeedPoints = 4;
	double StepSizeX = (VerCoo[3] - VerCoo[0]) / (numOfSeedPoints - 1);
	double StepSizeY = (VerCoo[7] - VerCoo[4]) / (numOfSeedPoints - 1);
	double StepSizeZ = (VerCoo[12] - VerCoo[0]) / (numOfSeedPoints - 1);

	double point[3];
	int currentDomainIndex = 1;
	bool check = true;
	int i, j, k;

	for (k = 0; k < numOfSeedPoints; k++)
	{
		for (j = 0; j < numOfSeedPoints; j++)
		{
			for (i = 0; i < numOfSeedPoints; i++)
			{
				point[0] = VerCoo[0] + i * StepSizeX;
				point[1] = VerCoo[1] + j * StepSizeY;
				point[2] = VerCoo[2] + k * StepSizeZ;

				currentDomainIndex = pEmDo->getDomainIndex(point);
				if (currentDomainIndex != 2)
				{
					check = false;		//一旦判断顶点的所在位置为 域外，就返回false
					goto BREAKLOOP2;
				}
			}
		}
	}
BREAKLOOP2:;
	return check;
}

//----------------------------判断胞元是否 被边界分割--------------------------------
bool Hexahedral::isIntersected(EmbeddedDomain *pEmDo, double* VerCoo)
{
	const int numOfSeedPoints = 4;
	double StepSizeX = (VerCoo[3] - VerCoo[0]) / (numOfSeedPoints - 1);
	double StepSizeY = (VerCoo[7] - VerCoo[4]) / (numOfSeedPoints - 1);
	double StepSizeZ = (VerCoo[12] - VerCoo[0]) / (numOfSeedPoints - 1);

	double point[3];
	int currentDomainIndex = 1;
	bool check = false;
	int i, j, k;

	point[0] = VerCoo[0];
	point[1] = VerCoo[1];
	point[2] = VerCoo[2];
	int initialDomainIndex = pEmDo->getDomainIndex(point);

	for (k = 0; k < numOfSeedPoints; k++)
	{
		for (j = 0; j < numOfSeedPoints; j++)
		{
			for (i = 0; i < numOfSeedPoints; i++)
			{
				point[0] = VerCoo[0] + i * StepSizeX;
				point[1] = VerCoo[1] + j * StepSizeY;
				point[2] = VerCoo[2] + k * StepSizeZ;

				currentDomainIndex = pEmDo->getDomainIndex(point);
				if (currentDomainIndex != initialDomainIndex)
				{
					check = true;		//一旦顶点的所在域的标号不同，就返回true
					goto BREAKLOOP3;
				}
			}
		}
	}
BREAKLOOP3:;
	return check;


}



//----------------------------返回胞元顶点坐标--------------------------------
double* Hexahedral::getCellVertexCoords() const
{
	double* VertexCoords = new double[NumNodes * 3];
	for (int i = 0; i < NumNodes; i++)
	{
		VertexCoords[3 * i] = nodes[i]->coords_[0];
		VertexCoords[3 * i + 1] = nodes[i]->coords_[1];
		VertexCoords[3 * i + 2] = nodes[i]->coords_[2];
	}
	return VertexCoords;
}
//----------------------------输出胞元顶点坐标--------------------------------
void Hexahedral::printCellVertexCoords()
{
	double* VC = getCellVertexCoords();
	printf("\nCell Vertex Coords: \n");
	for (int i = 0; i < 8; i++)
	{
		printf("Veterx %d : %-5.2f, %-5.2f, %-5.2f \n", i, VC[3 * i], VC[3 * i + 1], VC[3 * i + 2]);
	}

	delete[] VC; VC = nullptr;
}


//----------------------------将胞元内的局部坐标映射到全局坐标--------------------------------
double* Hexahedral::mapLocalToGlobal(double *Local, double *VC)
{
	// 这里的 Local 指的是积分点的坐标，也就是局部坐标
	// VC menas Vertices Coordinates
	double* global = new double[3];
	/*--------------------对于 普通六面体--------------------
	for (int i = 0; i < 3; i++)
	{
		global[i] = 0.125*((1 - Local[0])*(1 - Local[1])*(1 - Local[2])*VC[i]
			+ (1 + Local[0])*(1 - Local[1])*(1 - Local[2])*VC[3 + i]
			+ (1 + Local[0])*(1 + Local[1])*(1 - Local[2])*VC[6 + i]
			+ (1 - Local[0])*(1 + Local[1])*(1 - Local[2])*VC[9 + i]
			+ (1 - Local[0])*(1 - Local[1])*(1 + Local[2])*VC[12 + i]
			+ (1 + Local[0])*(1 - Local[1])*(1 + Local[2])*VC[15 + i]
			+ (1 + Local[0])*(1 + Local[1])*(1 + Local[2])*VC[18 + i]
			+ (1 - Local[0])*(1 + Local[1])*(1 + Local[2])*VC[21 + i]);
	}
	*/

	/* --------------------对于 长方体--------------------
	------------VC 数组的排放位置
				    x   y   z
	Vertices 1		0   1   2
	Vertices 2		3   4   5
	Vertices 3		6   7   8
	Vertices 4		9   10  11
	Vertices 5		12  13  14
	Vertices 6		15  16  17
	Vertices 7		18  19  20
	Vertices 8		21  22  23
	-------------------------------*/
	double lengthX = VC[3] - VC[0];
	double lengthY = VC[7] - VC[4];
	double lengthZ = VC[14] - VC[2];

	global[0] = 0.5*lengthX*(Local[0] + 1) + VC[0];
	global[1] = 0.5*lengthY*(Local[1] + 1) + VC[1];
	global[2] = 0.5*lengthZ*(Local[2] + 1) + VC[2];
	return global;
}


//----------------------------将全局坐标映射到胞元内的局部坐标--------------------------------
double* Hexahedral::mapGlobalToLocal(double *Global)
{
	double *VC = getCellVertexCoords();
	double* local = new double[3];

	double lengthX = VC[3] - VC[0];
	double lengthY = VC[7] - VC[4];
	double lengthZ = VC[14] - VC[2];

	local[0] = (Global[0] - VC[0]) * 2 / lengthX - 1;
	local[1] = (Global[1] - VC[1]) * 2 / lengthY - 1;
	local[2] = (Global[2] - VC[2]) * 2 / lengthZ - 1;
	
	delete[]VC; VC = nullptr;
	return local;
}



void  Hexahedral::printStiffnessMatrix(double* Kc)
{
	std::ofstream outfile("Stiffness.dat");
	if (!outfile)
	{
		std::cerr << "open Stiffness.dat error!" << std::endl;
		exit(1);
	}
	

	int NCNodes = (polynomialDegree + 1)*(polynomialDegree + 1)*(polynomialDegree + 1);
	int NCDof = NCNodes * 3;

	//----------输出单刚 Kc -----------
	printf("\n--------------------------Cell Stiffness--------------------------------\n");
	for (int i = 0; i < NCDof*NCDof; i++)
	{
		if ( (i) % NCDof == 0)
		{
			outfile << "\n";
			//printf("\n**%d**  ",i);
		}
		outfile << Kc[i] << " ";
		//printf("%8.3f  ", Kc[i]);
	}
	
	outfile.close();
	printf("\n");

}








