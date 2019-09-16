
#include "Element.h"
#include "LegendrePolynomial.h"

int* Hexahedral::getLocationMatrix()
{
	int numOfElementDofs = 3 * (polynomialDegree + 1)*(polynomialDegree + 1)
		*(polynomialDegree + 1);
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
double* Hexahedral::clacFace_dN(double *Coords, double *PhiR, double *PhiS, double *PhiT, double *dPhiR, double *dPhiS, double *dPhiT, double &OneMinusR, double &OneMinusS, double &OneMinusT, double &OnePlusR, double &OnePlusS, double &OnePlusT)
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
double* Hexahedral::clacBulk_dN(double *Coords, double *PhiR, double *PhiS, double *PhiT, double *dPhiR, double *dPhiS, double *dPhiT, double &OneMinusR, double &OneMinusS, double &OneMinusT, double &OnePlusR, double &OnePlusS, double &OnePlusT)
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
	double* faceN = calcEdge_N(Coords, phiR, phiS, phiT, oneMinusR, oneMinusS, oneMinusT, onePlusR, onePlusS, onePlusT);
	double* bulkN = calcEdge_N(Coords, phiR, phiS, phiT, oneMinusR, oneMinusS, oneMinusT, onePlusR, onePlusS, onePlusT);

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
		m = 3 * n + 8;							
		shapeN[m] = edgeN[n];					// 矩阵 N 的第 1 行
		shapeN[NCDof + m + 1] = edgeN[n];		// 矩阵 N 的第 2 行
		shapeN[2 *NCDof + m + 2] = edgeN[n];	// 矩阵 N 的第 3 行
	}

	// 面结点形函数
	for (n = 0; n < 6 * PerField * PerField; n++)
	{
		m =3 * n + 8 + 12 * PerField;
		shapeN[m] = faceN[n];					// 矩阵 N 的第 1 行
		shapeN[NCDof + m + 1] = faceN[n];		// 矩阵 N 的第 2 行
		shapeN[2 * NCDof + m + 2] = faceN[n];	// 矩阵 N 的第 3 行
	}

	// 体内结点形函数
	for (n = 0; n < PerField * PerField* PerField; n++)
	{
		m = 3 * n + 8 + 12 * PerField + 6 * PerField * PerField;
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
}



double* Hexahedral::calcMatrixB(double *Coords)
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

	// 计算每种模式下的形函数
	double* vertexdN = calcVertex_dN(Coords, oneMinusR, oneMinusS, oneMinusT, onePlusR, onePlusS, onePlusT);
	double* edgedN = calcEdge_dN(Coords, phiR, phiS, phiT, dphiR, dphiS, dphiT, oneMinusR, oneMinusS, oneMinusT, onePlusR, onePlusS, onePlusT);
	double* facedN = calcEdge_dN(Coords, phiR, phiS, phiT, dphiR, dphiS, dphiT, oneMinusR, oneMinusS, oneMinusT, onePlusR, onePlusS, onePlusT);
	double* bulkdN = calcEdge_dN(Coords, phiR, phiS, phiT, dphiR, dphiS, dphiT, oneMinusR, oneMinusS, oneMinusT, onePlusR, onePlusS, onePlusT);


	// 每个胞元节点数 NCNodes
	// 每个胞元自由度数 NCDof
	int NCNodes = (polynomialDegree + 1)*(polynomialDegree + 1)*(polynomialDegree + 1);
	int NCDof = NCNodes * 3;
	double* shapeN = new double[NCDof * 6]();			//后面加括号也表示调用默认构造函数，初始化为0

	// 计算形函数
	/*-----------------------------------------------------------------------------------------------------------------
			|  dN1_dx     0       0	    dN2_dx     0       0	...  edge dN  ...  face dN  ...  bulk dN  |
			|     0    dN1_dy     0	       0    dN2_dy     0	...  edge dN  ...  face dN  ...  bulk dN  |
		B =	|     0       0    dN1_dz	   0       0    dN2_dz  ...  edge dN  ...  face dN  ...  bulk dN  |
			|  dN1_dy  dN1_dx     0     dN2_dy  dN2_dx     0    ...  edge dN  ...  face dN  ...  bulk dN  |
			|     0    dN1_dz  dN1_dy      0    dN2_dz  dN2_dy  ...  edge dN  ...  face dN  ...  bulk dN  |
			|  dN1_dz    0     dN1_dx   dN2_dz    0     dN1_dx  ...  edge dN  ...  face dN  ...  bulk dN  |
	-----------------------------------------------------------------------------------------------------------------*/

	double *Jacobi;		//Jacobi矩阵









	delete[]phiR; phiR = nullptr;
	delete[]phiS; phiS = nullptr;
	delete[]phiS; phiT = nullptr;
	delete[]dphiR; dphiR = nullptr;
	delete[]dphiS; dphiS = nullptr;
	delete[]dphiS; dphiT = nullptr;
	delete[]vertexdN;	vertexdN = nullptr;
	delete[]edgedN;	edgedN = nullptr;
	delete[]facedN;	facedN = nullptr;
	delete[]bulkdN;  bulkdN = nullptr;



}













