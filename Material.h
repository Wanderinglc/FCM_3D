#ifndef MATERIAL_H
#define MATERIAL_H


class AbsMaterial
{
public:
	virtual double *getMaterialMatrix() = 0;

};


class Hooke3D :public AbsMaterial
{
public:
	Hooke3D(double E, double Poisson)
	{

	}

public:
	double MaterialMatrix[36];
	double Density;
	double ScalingFactor;
};



#endif
