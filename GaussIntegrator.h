#ifndef GaussIntegrator_H
#define GaussIntegrator_H


class GaussIntegrator
{
public:
	GaussIntegrator(int PD):NIP(PD+1) 	
	{
		GaussPoints = new double[NIP];
		GaussWeights = new double[NIP];
		formGaussCoordinates(GaussPoints);
		formGaussWeights(GaussWeights);
	}
	~GaussIntegrator()
	{	
		delete[]GaussPoints; GaussPoints = nullptr;
		delete[]GaussWeights; GaussWeights = nullptr;
	}

	void formGaussCoordinates(double* GaussPoints);
	void formGaussWeights(double* GaussWeights);

	double* getGaussCoordinates() { return GaussPoints; }
	double* getGaussWeights() { return GaussWeights; }



public:
	int NIP;	// Num of Integration Points
	double* GaussPoints;
	double* GaussWeights;
};





#endif
