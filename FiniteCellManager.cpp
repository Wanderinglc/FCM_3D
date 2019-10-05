#include <vector>
#include <algorithm> 
#include <iostream>
#include "FiniteCellManager.h"


void FiniteCellManager::initialize()
{
	/*后续用YAML Node 进行初始化*/



}



void FiniteCellManager::generateMesh()
{
	mesh->generateStructuredMesh();
}


double* FiniteCellManager::calcKsparse()
{
	const int &NEQ = mesh->totalDofs;								// 总方程数
	printf("NEQ =  %d\n", NEQ);
	int *MHT = new int[NEQ];										// 行宽数组，总刚中每一行的宽度
	int* Kdiag = new int[NEQ + 1];									// 对角元在一维总刚中的位置


	std::vector<int> columns;										// 非零元素在总刚中的列下标

	calcMHTandColumns(MHT, columns, mesh->polynomialDegree_);		// 形成 行宽数组 和列下标数组

	calcKdiag(Kdiag, MHT, NEQ);										// 计算对角元位置数组

	double* Kv = new double[NWK]();			// 总刚向量
	double* Kc = nullptr;					// 单刚向量
	int* LM = nullptr;
	// 按单元循环，计算 Kc, 然后调用 assembleKc_To_Kv() 把单刚组集到总刚的对应位置上

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

	//----------输出总刚 Kv -----------
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
	//--- 按胞元刚度矩阵的 行循环
	for (i = 0; i < NCDof; i++)
	{
		// i----单刚行号，  ii----单刚行号所对应的在总刚中的行号
		ii = LM[i];		
		// mi----对角元所在在的位置
		mi = Kdiag[ii];	
		for (j = i; j < NCDof; j++)
		{
			//  j----单刚列号，  	jj----单刚列号所对应的在总刚中的列号
			jj = LM[j];		
			//  总刚是一维存储，且只存上三角，按行存储
			//  如果 jj >= ii，说明是对应到总刚的上三角

			if (jj >= ii)
			{
				p = std::find(&Columns[mi], &Columns[Kdiag[ii + 1] - 1], jj);
				dist = int (std::distance(&Columns[mi], p));
				
				kk = mi + dist;					// 找到对角元所在位置后，kk---单刚中元素距离总刚中对角元的距离，		
				//printf("kk = %d\n",kk);
				Kv[kk] += Kc[NCDof*i + j];		// 注意单刚也是一维存储
			}

			//  如果 ii >= jj，说明是对应到总刚的下三角，此时就要找到对应的上三角的位置，对角元所在位置就应该取行号和列号中较大的那一个。
			else
			{
				// ni----对角元所在在的位置
				ni = Kdiag[jj];		
				p = std::find(&Columns[ni], &Columns[Kdiag[jj + 1] - 1], ii);
				dist = int (std::distance(&Columns[ni], p));

				kk = ni + dist;					// 找到对角元所在位置后，kk---单刚中元素距离总刚中对角元的距离，
				//printf("kk = %d\n", kk);
				Kv[kk] += Kc[NCDof*i + j];		// 注意单刚也是一维存储
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
	NWK = Kdiag[neq] - Kdiag[0];		// 总刚矩阵中的元素总数，也就是一维变带宽总刚的空间大小
	printf("NWK = %d\n", NWK);
}



void FiniteCellManager::calcMHTandColumns(int *MHT, std::vector<int>&columnVec, const int& PD)
{
	//const int& PD = mesh->polynomialDegree_;
	int cellDofs = mesh->NCDof;

	//std::vector<int> columnVec;
	std::vector<int> columnx;
	size_t nc = 0;
	
	int x, count;			// x表示x方向上的自由度编号

	for (int iNode = 0; iNode < mesh->totalNodes_; iNode++)
	{
		nc = mesh->nodeVector_[iNode].shareCellsId.size();								// 得到该结点所共享cell的个数
		const std::vector<int> &cellId = mesh->nodeVector_[iNode].shareCellsId;			// 得到该结点所共享cell的编号
		x = mesh->nodeVector_[iNode].dofs[0].id;

		int* shareLM = new int[nc * cellDofs];											//
		getShareCellLM(shareLM, cellId, cellDofs);										// 得到该结点所共享cell的 联系数组

		int m = 0, n = 0;
		std::sort(shareLM, shareLM + nc * cellDofs);									// 将该结点所在共享单元的的联系数组进行排序

		count = 0;
		for (m = 0; m < nc * cellDofs; m++)												// 排序完之后，找到比当前自由度编号 大 的自由度编号放入columnx中
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
		columnVec.insert(columnVec.end(), columnx.begin(), columnx.end());				// x方向上的列高数组		
		columnVec.insert(columnVec.end(), columnx.begin() + 1, columnx.end());			// y方向上的列高数组比x方向少1		
		columnVec.insert(columnVec.end(), columnx.begin() + 2, columnx.end());			// z方向上的列高数组比x方向少2
		columnx.clear();																// 把columnx清空

		MHT[x] = count - 1;			
		MHT[x + 1] = count - 2;		
		MHT[x + 2] = count - 3;
		delete[]shareLM;

	}
	
	int perField = PD - 1;
	for (int p = 0; p < perField; p++)
	{
		// 按边的第一阶循环
		for (int iEdge = 0; iEdge < mesh->totalEdges_; iEdge++)
		{
			nc = mesh->edgeVector_[iEdge].shareCellsId.size();								// 得到该结点所共享cell的个数
			const std::vector<int> &cellId = mesh->edgeVector_[iEdge].shareCellsId;			// 得到该结点所共享cell的编号
			x = mesh->edgeVector_[iEdge].dofs[3 * p].id;

			int* shareLM = new int[nc * cellDofs];											//
			getShareCellLM(shareLM, cellId, cellDofs);										// 得到该结点所共享cell的 联系数组

			size_t m = 0, n = 0;
			std::sort(shareLM, shareLM + nc * cellDofs);									// 将该结点所在共享单元的的联系数组进行排序

			count = 0;
			// 因为前8个顶点的共24个自由度，并且一定小于边上的自由度，所以 m = nc * 24; 
			for (m = nc * 24; m < nc * cellDofs; m++)										// 排序完之后，找到比当前自由度编号 大 的自由度编号放入columnx中
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
			columnVec.insert(columnVec.end(), columnx.begin(), columnx.end());				// x方向上的列高数组		
			columnVec.insert(columnVec.end(), columnx.begin() + 1, columnx.end());			// y方向上的列高数组比x方向少1		
			columnVec.insert(columnVec.end(), columnx.begin() + 2, columnx.end());			// z方向上的列高数组比x方向少2
			columnx.clear();																// 把columnx清空

			MHT[x] = count - 1;			
			MHT[x + 1] = count - 2;		
			MHT[x + 2] = count - 3;
			delete[]shareLM;
		}

		// 按 面 的第一阶循环
		for (int iFace = 0; iFace < mesh->totalFaces_; iFace++)
		{
			nc = mesh->faceVector_[iFace].shareCellsId.size();								// 得到该结点所共享cell的个数
			const std::vector<int> &cellId = mesh->faceVector_[iFace].shareCellsId;			// 得到该结点所共享cell的编号

			int* shareLM = new int[nc * cellDofs];											//
			getShareCellLM(shareLM, cellId, cellDofs);										// 得到该结点所共享cell的 联系数组

			size_t m = 0, n = 0;
			std::sort(shareLM, shareLM + nc * cellDofs);									// 将该结点所在共享单元的的联系数组进行排序

			int start = p * p;
			int end = (p + 1)*(p + 1);

			for (int j = start; j < end; j++)
			{
				count = 0;
				x = mesh->faceVector_[iFace].dofs[3 * j].id;
				// 因为前8个顶点的共24个自由度，12条边36个自由度，并且一定小于面上的自由度，所以 m = nc*(24+36); 
				for (m = nc * 60; m < nc * cellDofs; m++)												// 排序完之后，找到比当前自由度编号 大 的自由度编号放入columnx中
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
				columnVec.insert(columnVec.end(), columnx.begin(), columnx.end());				// x方向上的列高数组		
				columnVec.insert(columnVec.end(), columnx.begin() + 1, columnx.end());			// y方向上的列高数组比x方向少1		
				columnVec.insert(columnVec.end(), columnx.begin() + 2, columnx.end());			// z方向上的列高数组比x方向少2
				columnx.clear();																// 把columnx清空

				MHT[x] = count - 1;			
				MHT[x + 1] = count - 2;		
				MHT[x + 2] = count - 3;				
			}
			delete[]shareLM;
		}

		// 按 体 的第一阶循环
		for (int iBulk = 0; iBulk < mesh->totalBulks_; iBulk++)
		{
			nc = 1;																				// 该结点所共享cell的个数
			const int &cellId = mesh->bulkVector_[iBulk].shareCellsId;							// 得到该结点所共享cell的编号
								
			int* shareLM = mesh->elementVector_[cellId].getLocationMatrix();					// 得到该结点所共享cell的 联系数组

			size_t m = 0, n = 0;
			std::sort(shareLM, shareLM + nc * cellDofs);										// 将该结点所在共享单元的的联系数组进行排序

			int start = p * p * p;
			int end = (p + 1)*(p + 1)*(p + 1);

			for (int j = start; j < end; j++)
			{
				count = 0;
				x = mesh->bulkVector_[iBulk].dofs[3 * j].id;
				// 因为前8个顶点的共24个自由度，12条边至少36个自由度，6个面上至少18个自由度，并且一定小于体内的自由度，所以 m = nc*(24+36+18); 
				for (m = nc * 78; m < nc * cellDofs; m++)										// 排序完之后，找到比当前自由度编号 大 的自由度编号放入columnx中
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
				columnVec.insert(columnVec.end(), columnx.begin(), columnx.end());				// x方向上的列高数组		
				columnVec.insert(columnVec.end(), columnx.begin() + 1, columnx.end());			// y方向上的列高数组比x方向少1		
				columnVec.insert(columnVec.end(), columnx.begin() + 2, columnx.end());			// z方向上的列高数组比x方向少2
				columnx.clear();																// 把columnx清空

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
	const int &NEQ = mesh->totalDofs;								// 总方程数
	double* Fv = new double[NEQ]();									// 总载荷向量

	NeumannBoundaryCondition* NeuBC;								// 载荷边界条件指针
	for (size_t nBC = 0; nBC < NeumannBC.size(); nBC++)
	{
		NeuBC = NeumannBC[nBC];
		
		NeuBC->calcLoadVector(mesh, Fv);
	}



}

