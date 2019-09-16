#ifndef DOF_H
#define DOF_H

/*=======================================================================
						 Dof
  Define the attributes of degrees of freedom for finite element method.
 =======================================================================*/

class Dof
{
	// ��Ա����
public:
	// constructor
	Dof() :id(0), value(0) { }
	Dof(int Id) :id(Id), value(0) { }
	Dof(int Id, double Value) :id(Id), value(Value) { }
	// Destructor
	virtual ~Dof()
	{
		// nothing for now 
	}

	int getId() const { return id; }
	int setId(int Id) { id = Id; }

	double getValue() const { return value; }
	double setValue(double a) { value = a; }

	// ��Ա����
public:
	int id;
	double value;
};


//Dof::Dof() { id = 0; value = 0.0; }
//Dof::Dof(int Id) { id = Id; value = 0.0; }
//Dof::Dof(int Id, double Value) :id(Id), value(Value) {};


#endif