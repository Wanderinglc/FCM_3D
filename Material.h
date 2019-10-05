#ifndef MATERIAL_H
#define MATERIAL_H

class AbsMaterial
{
public:
	AbsMaterial() {}
	virtual ~AbsMaterial()
	{
		if (MaterialMatrix != nullptr)
		{
			delete MaterialMatrix;
			MaterialMatrix = nullptr;
		}
	}


	virtual double *getMaterialMatrix()
	{	return MaterialMatrix;	}

public:
	double* MaterialMatrix;
};


class Hooke3D :public AbsMaterial
{
public:

	Hooke3D() { }

	Hooke3D(double E, double v, double d)
		:Density(d)
	{
		MaterialMatrix = new double[36]();

		MaterialMatrix[0] = E / (1 + v) / (1 - 2 * v) * (1 - v);
		MaterialMatrix[1] = E / (1 + v) / (1 - 2 * v) * v;
		MaterialMatrix[2] = MaterialMatrix[1];
		MaterialMatrix[6] = MaterialMatrix[1];
		MaterialMatrix[7] = MaterialMatrix[0];
		MaterialMatrix[8] = MaterialMatrix[1];
		MaterialMatrix[12] = MaterialMatrix[1];
		MaterialMatrix[13] = MaterialMatrix[1];
		MaterialMatrix[14] = MaterialMatrix[0];
		MaterialMatrix[21] = E / (1 + v) / 2;
		MaterialMatrix[28] = MaterialMatrix[21];
		MaterialMatrix[35] = MaterialMatrix[21];
	}

	virtual ~Hooke3D()
	{

	}

	void printMaterialMatrix() const
	{
		printf("\nMaterial Matrix D: ");
		for (int i = 0; i < 36; i++)
		{
			if (i % 6 == 0)  printf("| \n |  ");
			printf("%-7.2f  ", MaterialMatrix[i]);			
		}
		printf("| \n");
	}

public:	
	double Density;
	//double ScalingFactor;
};



#endif
