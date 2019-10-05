#include <vector>
#include <algorithm> 
#include <iostream>
#include "FiniteCellManager.h"


void FiniteCellManager::initialize()
{
	/*������YAML Node ���г�ʼ��*/



}



void FiniteCellManager::generateMesh()
{
	mesh->generateStructuredMesh();
}


double* FiniteCellManager::calcKsparse()
{
	const int &NEQ = mesh->totalDofs;								// �ܷ�����
	printf("NEQ =  %d\n", NEQ);
	int *MHT = new int[NEQ];										// �п����飬�ܸ���ÿһ�еĿ��
	int* Kdiag = new int[NEQ + 1];									// �Խ�Ԫ��һά�ܸ��е�λ��


	std::vector<int> columns;										// ����Ԫ�����ܸ��е����±�

	calcMHTandColumns(MHT, columns, mesh->polynomialDegree_);		// �γ� �п����� �����±�����

	calcKdiag(Kdiag, MHT, NEQ);										// ����Խ�Ԫλ������

	double* Kv = new double[NWK]();			// �ܸ�����
	double* Kc = nullptr;					// ��������
	int* LM = nullptr;
	// ����Ԫѭ�������� Kc, Ȼ����� assembleKc_To_Kv() �ѵ����鼯���ܸյĶ�Ӧλ����

	for (int iele = 0; iele < mesh->totalCells_; iele++)
	{		
		LM = mesh->elementVector_[iele].getLocationMatrix();
		Kc = mesh->elementVector_[iele].calcStiffnessMatrix(domain, material, integrator, spaceTreeDepth);
		assemble_Kc_To_Kv(Kc, Kv, Kdiag, columns, LM, mesh->NCDof);
		delete[]LM;
		delete[]Kc;
	}

	printStiffness(Kv, Kdiag, NEQ);

	return Kv;

}


void FiniteCellManager::printStiffness(const double* Kv, const int* Kdiag, const int NEQ) const
{
	printf("\n");
	int j;

	//----------����ܸ� Kv -----------
	for (int i = 0; i < NEQ; i++)
	{
		j = Kdiag[i];
		printf("K[%d,%d] = %-10.3f  ", i+1, i+1, Kv[j]);
		std::cout << std::endl;
	}

}



void FiniteCellManager::assemble_Kc_To_Kv(double *Kc, double *Kv, int *Kdiag, std::vector<int> &Columns, int *LM, const int NCDof)
{

	// Kv ---Global Stiffness (1D Compressed Sparse Row(CSR) storage)
	// kc ---Element Stiffness
	int i, ii, mi, ni, j, jj, kk;

	kk = 0;
	int* p = nullptr;
	int dist;
	//--- ����Ԫ�նȾ���� ��ѭ��
	for (i = 0; i < NCDof; i++)
	{
		// i----�����кţ�  ii----�����к�����Ӧ�����ܸ��е��к�
		ii = LM[i];		
		// mi----�Խ�Ԫ�����ڵ�λ��
		mi = Kdiag[ii];	
		for (j = i; j < NCDof; j++)
		{
			//  j----�����кţ�  	jj----�����к�����Ӧ�����ܸ��е��к�
			jj = LM[j];		
			//  �ܸ���һά�洢����ֻ�������ǣ����д洢
			//  ��� jj >= ii��˵���Ƕ�Ӧ���ܸյ�������

			if (jj >= ii)
			{
				p = std::find(&Columns[mi], &Columns[Kdiag[ii + 1] - 1], jj);
				dist = int (std::distance(&Columns[mi], p));
				
				kk = mi + dist;					// �ҵ��Խ�Ԫ����λ�ú�kk---������Ԫ�ؾ����ܸ��жԽ�Ԫ�ľ��룬		
				//printf("kk = %d\n",kk);
				Kv[kk] += Kc[NCDof*i + j];		// ע�ⵥ��Ҳ��һά�洢
			}

			//  ��� ii >= jj��˵���Ƕ�Ӧ���ܸյ������ǣ���ʱ��Ҫ�ҵ���Ӧ�������ǵ�λ�ã��Խ�Ԫ����λ�þ�Ӧ��ȡ�кź��к��нϴ����һ����
			else
			{
				// ni----�Խ�Ԫ�����ڵ�λ��
				ni = Kdiag[jj];		
				p = std::find(&Columns[ni], &Columns[Kdiag[jj + 1] - 1], ii);
				dist = int (std::distance(&Columns[ni], p));

				kk = ni + dist;					// �ҵ��Խ�Ԫ����λ�ú�kk---������Ԫ�ؾ����ܸ��жԽ�Ԫ�ľ��룬
				//printf("kk = %d\n", kk);
				Kv[kk] += Kc[NCDof*i + j];		// ע�ⵥ��Ҳ��һά�洢
			}
		}
	}
}




void FiniteCellManager::calcKdiag(int *Kdiag, int *MHT, const int neq)
{
	Kdiag[0] = 0;
	//printf("Kdiag[%d] = %d\n", 0, Kdiag[0]);
	for (int i = 0; i < neq; i++)
	{
		Kdiag[i + 1] = Kdiag[i] + MHT[i] + 1;
		//printf("Kdiag[%d] = %d\n", i+1, Kdiag[i+1]);
	}
	NWK = Kdiag[neq] - Kdiag[0];		// �ܸվ����е�Ԫ��������Ҳ����һά������ܸյĿռ��С
	printf("NWK = %d\n", NWK);
}



void FiniteCellManager::calcMHTandColumns(int *MHT, std::vector<int>&columnVec, const int& PD)
{
	//const int& PD = mesh->polynomialDegree_;
	int cellDofs = mesh->NCDof;

	//std::vector<int> columnVec;
	std::vector<int> columnx;
	size_t nc = 0;
	
	int x, count;			// x��ʾx�����ϵ����ɶȱ��

	for (int iNode = 0; iNode < mesh->totalNodes_; iNode++)
	{
		nc = mesh->nodeVector_[iNode].shareCellsId.size();								// �õ��ý��������cell�ĸ���
		const std::vector<int> &cellId = mesh->nodeVector_[iNode].shareCellsId;			// �õ��ý��������cell�ı��
		x = mesh->nodeVector_[iNode].dofs[0].id;

		int* shareLM = new int[nc * cellDofs];											//
		getShareCellLM(shareLM, cellId, cellDofs);										// �õ��ý��������cell�� ��ϵ����

		int m = 0, n = 0;
		std::sort(shareLM, shareLM + nc * cellDofs);									// ���ý�����ڹ���Ԫ�ĵ���ϵ�����������

		count = 0;
		for (m = 0; m < nc * cellDofs; m++)												// ������֮���ҵ��ȵ�ǰ���ɶȱ�� �� �����ɶȱ�ŷ���columnx��
		{
			if (shareLM[m] >= x)
			{
				if (shareLM[m] != shareLM[m - 1])
				{
					count++;
					columnx.push_back(shareLM[m]);
				}
			}
		}
		columnVec.insert(columnVec.end(), columnx.begin(), columnx.end());				// x�����ϵ��и�����		
		columnVec.insert(columnVec.end(), columnx.begin() + 1, columnx.end());			// y�����ϵ��и������x������1		
		columnVec.insert(columnVec.end(), columnx.begin() + 2, columnx.end());			// z�����ϵ��и������x������2
		columnx.clear();																// ��columnx���

		MHT[x] = count - 1;			
		MHT[x + 1] = count - 2;		
		MHT[x + 2] = count - 3;
		delete[]shareLM;

	}
	
	int perField = PD - 1;
	for (int p = 0; p < perField; p++)
	{
		// ���ߵĵ�һ��ѭ��
		for (int iEdge = 0; iEdge < mesh->totalEdges_; iEdge++)
		{
			nc = mesh->edgeVector_[iEdge].shareCellsId.size();								// �õ��ý��������cell�ĸ���
			const std::vector<int> &cellId = mesh->edgeVector_[iEdge].shareCellsId;			// �õ��ý��������cell�ı��
			x = mesh->edgeVector_[iEdge].dofs[3 * p].id;

			int* shareLM = new int[nc * cellDofs];											//
			getShareCellLM(shareLM, cellId, cellDofs);										// �õ��ý��������cell�� ��ϵ����

			size_t m = 0, n = 0;
			std::sort(shareLM, shareLM + nc * cellDofs);									// ���ý�����ڹ���Ԫ�ĵ���ϵ�����������

			count = 0;
			// ��Ϊǰ8������Ĺ�24�����ɶȣ�����һ��С�ڱ��ϵ����ɶȣ����� m = nc * 24; 
			for (m = nc * 24; m < nc * cellDofs; m++)										// ������֮���ҵ��ȵ�ǰ���ɶȱ�� �� �����ɶȱ�ŷ���columnx��
			{
				if (shareLM[m] >= x)
				{
					if (shareLM[m] != shareLM[m - 1])
					{
						count++;
						columnx.push_back(shareLM[m]);
					}
				}
			}
			columnVec.insert(columnVec.end(), columnx.begin(), columnx.end());				// x�����ϵ��и�����		
			columnVec.insert(columnVec.end(), columnx.begin() + 1, columnx.end());			// y�����ϵ��и������x������1		
			columnVec.insert(columnVec.end(), columnx.begin() + 2, columnx.end());			// z�����ϵ��и������x������2
			columnx.clear();																// ��columnx���

			MHT[x] = count - 1;			
			MHT[x + 1] = count - 2;		
			MHT[x + 2] = count - 3;
			delete[]shareLM;
		}

		// �� �� �ĵ�һ��ѭ��
		for (int iFace = 0; iFace < mesh->totalFaces_; iFace++)
		{
			nc = mesh->faceVector_[iFace].shareCellsId.size();								// �õ��ý��������cell�ĸ���
			const std::vector<int> &cellId = mesh->faceVector_[iFace].shareCellsId;			// �õ��ý��������cell�ı��

			int* shareLM = new int[nc * cellDofs];											//
			getShareCellLM(shareLM, cellId, cellDofs);										// �õ��ý��������cell�� ��ϵ����

			size_t m = 0, n = 0;
			std::sort(shareLM, shareLM + nc * cellDofs);									// ���ý�����ڹ���Ԫ�ĵ���ϵ�����������

			int start = p * p;
			int end = (p + 1)*(p + 1);

			for (int j = start; j < end; j++)
			{
				count = 0;
				x = mesh->faceVector_[iFace].dofs[3 * j].id;
				// ��Ϊǰ8������Ĺ�24�����ɶȣ�12����36�����ɶȣ�����һ��С�����ϵ����ɶȣ����� m = nc*(24+36); 
				for (m = nc * 60; m < nc * cellDofs; m++)												// ������֮���ҵ��ȵ�ǰ���ɶȱ�� �� �����ɶȱ�ŷ���columnx��
				{
					if (shareLM[m] >= x)
					{
						if (shareLM[m] != shareLM[m - 1])
						{
							count++;
							columnx.push_back(shareLM[m]);
						}
					}
				}
				columnVec.insert(columnVec.end(), columnx.begin(), columnx.end());				// x�����ϵ��и�����		
				columnVec.insert(columnVec.end(), columnx.begin() + 1, columnx.end());			// y�����ϵ��и������x������1		
				columnVec.insert(columnVec.end(), columnx.begin() + 2, columnx.end());			// z�����ϵ��и������x������2
				columnx.clear();																// ��columnx���

				MHT[x] = count - 1;			
				MHT[x + 1] = count - 2;		
				MHT[x + 2] = count - 3;				
			}
			delete[]shareLM;
		}

		// �� �� �ĵ�һ��ѭ��
		for (int iBulk = 0; iBulk < mesh->totalBulks_; iBulk++)
		{
			nc = 1;																				// �ý��������cell�ĸ���
			const int &cellId = mesh->bulkVector_[iBulk].shareCellsId;							// �õ��ý��������cell�ı��
								
			int* shareLM = mesh->elementVector_[cellId].getLocationMatrix();					// �õ��ý��������cell�� ��ϵ����

			size_t m = 0, n = 0;
			std::sort(shareLM, shareLM + nc * cellDofs);										// ���ý�����ڹ���Ԫ�ĵ���ϵ�����������

			int start = p * p * p;
			int end = (p + 1)*(p + 1)*(p + 1);

			for (int j = start; j < end; j++)
			{
				count = 0;
				x = mesh->bulkVector_[iBulk].dofs[3 * j].id;
				// ��Ϊǰ8������Ĺ�24�����ɶȣ�12��������36�����ɶȣ�6����������18�����ɶȣ�����һ��С�����ڵ����ɶȣ����� m = nc*(24+36+18); 
				for (m = nc * 78; m < nc * cellDofs; m++)										// ������֮���ҵ��ȵ�ǰ���ɶȱ�� �� �����ɶȱ�ŷ���columnx��
				{
					if (shareLM[m] >= x)
					{
						if (shareLM[m] != shareLM[m - 1])
						{
							count++;
							columnx.push_back(shareLM[m]);
						}
					}
				}
				columnVec.insert(columnVec.end(), columnx.begin(), columnx.end());				// x�����ϵ��и�����		
				columnVec.insert(columnVec.end(), columnx.begin() + 1, columnx.end());			// y�����ϵ��и������x������1		
				columnVec.insert(columnVec.end(), columnx.begin() + 2, columnx.end());			// z�����ϵ��и������x������2
				columnx.clear();																// ��columnx���

				MHT[x] = count - 1;			//printf("MHT[%d] =  %d \n", x, MHT[x - 1]);
				MHT[x + 1] = count - 2;		//printf("MHT[%d] =  %d \n", y, MHT[y - 1]);
				MHT[x + 2] = count - 3;
			}
			delete[]shareLM;
		}

	}

}



void FiniteCellManager::getShareCellLM(int* shareLM, const std::vector<int> &cellId, const int cellDofs)
{	
	size_t nc = cellId.size();

	int *LM = nullptr;
	for (int n = 0; n < nc; n++)
	{
		LM = mesh->elementVector_[cellId[n]].getLocationMatrix();
		for (int i = 0; i < cellDofs; i++)
		{
			shareLM[n*cellDofs + i] = LM[i];
		}
		delete[]LM;
	}

}



double* FiniteCellManager::calcLoadVector()
{
	const int &NEQ = mesh->totalDofs;								// �ܷ�����
	double* Fv = new double[NEQ]();									// ���غ�����

	NeumannBoundaryCondition* NeuBC;								// �غɱ߽�����ָ��
	for (size_t nBC = 0; nBC < NeumannBC.size(); nBC++)
	{
		NeuBC = NeumannBC[nBC];
		
		NeuBC->calcLoadVector(mesh, Fv);
	}



}

