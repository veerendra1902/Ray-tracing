

#ifndef VECTOROPERATION_H_
#define VECTOROPERATION_H_

#include "vector3D.h"

	void printVector(vector3D v);
	vector3D add(vector3D v1,vector3D v2);
	vector3D subtract(vector3D v1,vector3D v2);
	float dotProduct(vector3D v1,vector3D v2);
	float angle(vector3D v1,vector3D v2);
	vector3D crossProduct(vector3D v1,vector3D v2);
	vector3D constProduct(vector3D v1,float t);
	float magnitudeSqr(vector3D v);
	float magnitude(vector3D v);
	vector3D normalise(vector3D v);
	float angx(vector3D v);//gives the cos(angle) with x axis
	float angy(vector3D v);
	float angz(vector3D v);
	vector3D rotatevectorLeft(vector3D v,vector3D n,float angle);
	vector3D rotatevectorRight(vector3D v,vector3D n,float angle);///  it will return a unit vector



#endif /* VECTOROPERATION_H_ */
