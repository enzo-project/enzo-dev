
#ifndef ____SuperNova__
#define ____SuperNova__

#include <stdio.h>


typedef struct {
	FLOAT dbx;		// x component of magnetic flux
	FLOAT dby;		// y component of magnetic flux
	FLOAT dbz;		// z component of magnetic flux
	FLOAT bx;
	FLOAT by;
	FLOAT bz;
	FLOAT dUb;
	FLOAT dUtot;
	
} snsf_source_terms;

class SuperNova
{
public:
	SuperNova();
	~SuperNova() {}
	void setValues(FLOAT phi_x, FLOAT phi_y, FLOAT phi_z,  FLOAT x, FLOAT y, FLOAT z, \
				   FLOAT radius, FLOAT time_started, FLOAT duration, float energy, float sigma_sn);
	snsf_source_terms getSourceTerms(double dx, double dy, double dz, double enzoTime);
	
        FLOAT* getPosition();
private:
	FLOAT location[3];
	FLOAT zHat[3];
	FLOAT characteristicLength;
	FLOAT timeStarted;
	FLOAT characteristicTime;
	
	float totalEnergy;
	float sigma;
};



#endif /* defined(____SuperNova__) */
