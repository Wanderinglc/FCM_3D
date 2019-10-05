#ifndef EmbeddedDomain_H
#define EmbeddedDomain_H

#include <cmath>

class EmbeddedDomain
{
public:
	
	EmbeddedDomain(double a):alpha(a) { }		//  构造函数
	virtual ~EmbeddedDomain(){}				//  析构函数

	virtual int getDomainIndex(double* Coords) = 0;
	virtual double getMaterialAtPoint(double* Coords) = 0;

public:
	double alpha;		// ScalingFactor 
};


class ThickPlateWithAHole :public EmbeddedDomain
{
public:
	//  构造函数
	ThickPlateWithAHole(double a):EmbeddedDomain(a){}

	ThickPlateWithAHole(double* start, double* length, double* center, double radius, double a)
		: holeRadius(radius), EmbeddedDomain(a)
	{
		origin[0] = start[0];
		origin[1] = start[1];
		origin[2] = start[2];

		lengths[0] = length[0];
		lengths[1] = length[1];
		lengths[2] = length[2];
		
		holeCenter[0] = center[0];
		holeCenter[1] = center[1];
	}

	//  析构函数
	virtual ~ThickPlateWithAHole() {}
	
	
	//  指示函数，用来判断点的位置，是在物理域外还是物理域内
	virtual int getDomainIndex(double* Coords)
	{
		// domaiIndex = 1 : Outside
		// domaiIndex = 2 : InsidePlate
		double tolerance = 1E-15;
		int domainIndex = 0;
		
		double diffx = Coords[0] - origin[0] + tolerance;
		double diffy = Coords[1] - origin[1] + tolerance;
		double diffz = Coords[2] - origin[2] + tolerance;

		double end[3];
		end[0] = origin[0] + lengths[0] + tolerance;
		end[1] = origin[1] + lengths[1] + tolerance;
		end[2] = origin[2] + lengths[2] + tolerance;

		//  Check if inside box
		if (diffx > 0 && Coords[0] < end[0] && diffy>0 && Coords[1] < end[1]
			&& diffz>0 && Coords[2] < end[2])
		{
			//  Check if outside hole
			double distance = hypot(Coords[0] - holeCenter[0], Coords[1] - holeCenter[1]);
			if (distance >= holeRadius)
				domainIndex = 2;	
			else
				domainIndex = 1;
		}
		else
			domainIndex = 1;
		return domainIndex;
	}


	//  用来判断点的位置，并返回这一点的材料惩罚系数
	virtual double getMaterialAtPoint(double* Coords)
	{
		double scalingFactor;
		int domainIndex = getDomainIndex(Coords);
		if (domainIndex == 2)
			scalingFactor = 1.0;
		else
			scalingFactor = alpha;
		return scalingFactor;
	}


public:
	double origin[3];
	double lengths[3];
	double holeCenter[2];
	double holeRadius;
};



#endif
