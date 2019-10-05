#ifndef StrongPenaltyAlgorithm_H
#define StrongPenaltyAlgorithm_H


class StrongPenaltyAlgorithm
{
public:
	StrongPenaltyAlgorithm(double Penalty):penaltyValue(Penalty) {	}
	
	~StrongPenaltyAlgorithm()
	{

	}

	void modifyLinearSystem(int dofId, double prescribedValues, double* K, int* Kdiag, double* F)
	{
		K[Kdiag[dofId]] += penaltyValue;
		F[dofId] += penaltyValue * prescribedValues;
	}


public:
	double penaltyValue;

};






#endif