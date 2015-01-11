/*
 * vectorOperation.cpp
 *
 *  Created on: Aug 24, 2014
 *      Author: Veerendra
 */
#include<Math.h>
#include<iostream>
#include<stdio.h>
#include "vector3D.h"
#include "vectorOperation.h"

using namespace std;


const float pi=3.14;


float dotProduct(vector3D v1, vector3D v2){
		float dot=(v1.x*v2.x)+(v1.y*v2.y)+(v1.z*v2.z);
		return dot;
	}

	vector3D crossProduct(vector3D v1, vector3D v2){



		return vector3D(v1.y*v2.z-v2.y*v1.z,v2.x*v1.z-v1.x*v2.z,v1.x*v2.y-v2.x*v1.y);
	}

	float magnitudeSqr(vector3D v){

		float result=((v.x)*(v.x)+(v.y)*(v.y)+(v.z)*(v.z));
		return result;
	}

	float magnitude(vector3D v){

			float result=sqrt((v.x)*(v.x)+(v.y)*(v.y)+(v.z)*(v.z));
			return result;
		}

	vector3D add(vector3D v1,vector3D v2){



		return vector3D(v1.x+v2.x,v1.y+v2.y,v1.z+v2.z);
	}

	vector3D subtract(vector3D v1,vector3D v2){
		return vector3D(v1.x-v2.x,v1.y-v2.y,v1.z-v2.z);

		}

	vector3D normalise(vector3D v){

		float r=sqrt((v.x)*(v.x)+(v.y)*(v.y)+(v.z)*(v.z));

vector3D n(v.x/r,v.y/r,v.z/r);
		return n;


	}

	void printVector(vector3D v){

		cout<<"x coord:"<< v.x<<endl; ;
		cout<<"y coord:"<< v.y<<endl;;
		cout<< "z coord:"<<v.z<<endl;
		cout<<endl;
	}

	float angle(vector3D v1,vector3D v2){
		float a=magnitude(v1);
		float b=magnitude(v2);
		float d=dotProduct(v1,v2);
		return  d/(a*b);

	}


	float angx(vector3D v){
		float r=sqrt((v.x)*(v.x)+(v.y)*(v.y)+(v.z)*(v.z));
		float a=v.x/r;
		return a;
	}


	float angy(vector3D v){
			float r=sqrt((v.x)*(v.x)+(v.y)*(v.y)+(v.z)*(v.z));
			float a=v.y/r;
			return a;
		}

	float angz(vector3D v){
			float r=sqrt((v.x)*(v.x)+(v.y)*(v.y)+(v.z)*(v.z));
			float a=v.z/r;
			return a;
		}


	vector3D constProduct(vector3D v,float t){

	//	vector3D n(v.x,v.y,v.z);

		return vector3D(v.x*t,v.y*t,v.z*t);
		//return n;
	}

	vector3D rotateVectorLeft(vector3D v,vector3D n,float angle){
		vector3D i=normalise(v);
	//	v.Printvector3D(i);
		float a= (pi/180)*angle;
	//	cout<<cos(a)<<endl;
	//	v.Printvector3D(v.ConstProduct(i,cos(a)));
		return subtract(constProduct(i,cos(a)),crossProduct(i,constProduct(n,sin(a))));

	}

	vector3D rotateVectorRight(vector3D v,vector3D n,float angle){
			vector3D i=normalise(v);
		//	v.Printvector3D(i);
			float a= (pi/180)*angle;
		//	cout<<cos(a)<<endl;
//			v.Printvector3D(v.ConstProduct(i,cos(a)));
//			v.Printvector3D(i);
//			v.Printvector3D(v.CrossProduct(i,v.ConstProduct(n,-1)));
			return subtract(constProduct(i,cos(a)),crossProduct(i,constProduct(n,-sin(a))));


	}

