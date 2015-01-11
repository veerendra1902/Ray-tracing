#include <iostream>
#include <windows.h>
#include <GL/glut.h>
#include <math.h>
#include <vector>
#include <bits/stdc++.h>
#include "vector3D.h"
#include "vectorOperation.h"
#include "matrixFunctions.h"

using namespace std;

GLfloat a=0;
vector3D normal;//(-10,0,-10);
vector3D xyzUp;//(0.0,1.0,0.0);
vector3D n;
vector3D u;
vector3D v;
vector3D origin;//(2,0,2);
vector3D vcsEye;//(0,0,-3);
vector3D nearPlane;//(0,0,-2);
vector3D farPlane;//(0,0,15);
vector3D viewC1;//(1,1,0);
vector3D viewC2;//(-1,1,0);
vector3D viewC3;//(-1,-1,0);
vector3D viewC4;//(1,-1,0);
vector3D viewPlane[4];//={viewC1,viewC2,viewC3,viewC4};

GLfloat M[4][4];
GLfloat matrixOrigin[4][4];

GLfloat xyzTuvn[16];
GLfloat wcsTvcs[4][4];
//GLfloat uvnTxyz[16];
GLfloat inv[4][4];
GLfloat perspect[4][4];
GLfloat perspective[16];




struct light{
	vector3D position;
	vector3D intensity;

};



struct objects{
		vector3D center;
		float radius;
		string type;
		int id;

		vector3D normal;
		float D;

		float kaR;
		float kaG;
		float kaB;

		float kdR;
		float kdG;
		float kdB;

		float ksR;
		float ksG;
		float ksB;

		float reflectionCoff;
		float refractionCoff;

		float n;
};

struct point{
	objects obj;
	bool intersected;
	float t;
};

struct pixel{
	vector3D color;
	vector3D position;
};

vector <vector3D> variables ;
vector< objects > object;
vector<light > light_sources ;
vector3D ambient;

void init(void)
{
glClearColor(0,0,0,0);

}



void matrixMult(GLfloat matrixA[4][4],GLfloat matrixB[4][4],GLfloat matrixC[4][4],int mA,int nA,int nB){
	for(int i=0;i<mA;i++){
	for(int j=0;j<nB;j++){
	        matrixC[i][j]=0;
	        for(int k=0;k<nA;k++){
	        	matrixC[i][j]=matrixC[i][j]+(matrixA[i][k] * matrixB[k][j]);
	}

	}
	}

}

void matrixAdd(GLfloat matrixA[4][4],GLfloat matrixB[4][4],GLfloat matrixC[4][4],int n,int m){
	for(int i=0;i<n;i++){
		for(int j=0;j<m;j++){
			matrixC[i][j]=matrixA[i][j]+matrixB[i][j];
		}
	}

}



vector3D convert(vector3D v){

	GLfloat m1[4][4];
		m1[0][0]=v.x;
		m1[0][1]=v.y;
		m1[0][2]=v.z;
	GLfloat m3[4][4];
	GLfloat m4[4][4];
	vector3D result;
	matrixMult(m1,M,m3,1,3,3);

	matrixAdd(m3,matrixOrigin,m4,1,3);

		result.x=m4[0][0];
		result.y=m4[0][1];
		result.z=m4[0][2];

return result;

}

void drawPlane(vector3D* p){
	vector3D start=convert(p[0]);
		for(int i=1;i<4;i++){
			vector3D end=convert(p[i]);
			glBegin(GL_LINES);
			glVertex3f(start.x,start.y,start.z);
			glVertex3f(end.x,end.y,end.z);
			glEnd();

			start=end;

		}
		vector3D end=convert(p[0]);
		glBegin(GL_LINES);
			glVertex3f(start.x,start.y,start.z);
			glVertex3f(end.x,end.y,end.z);
		glEnd();

}

void drawFrustam(vector3D* view){
	vector3D e=convert(vcsEye);
	drawPlane(view);
	vector3D nearP[4];
	float nDist=(float)(abs(vcsEye.z)-abs(nearPlane.z))/(float)abs(vcsEye.z);
	float fDist=(float)(abs(vcsEye.z)+abs(farPlane.z))/(float)abs(vcsEye.z);
	for(int i=0;i<4;i++){
		vector3D nearX1(view[i].x*nDist,view[i].y*nDist,nearPlane.z);
		nearP[i]=nearX1;
	}
//	drawPlane(nearP);
	vector3D farP[4];
		for(int i=0;i<4;i++){
			vector3D farX1(view[i].x*fDist,view[i].y*fDist,farPlane.z);
			vector3D line=convert(farX1);
			   glBegin(GL_LINES);
			       glVertex3f(e.x,e.y,e.z);
			       glVertex3f(line.x,line.y,line.z);
			   glEnd();
			farP[i]=farX1;
		}
//	drawPlane(farP);
}


void read(){
	ifstream myfile("input.txt");
	string line_str ;
	if (myfile.is_open())
	{
		while(myfile.good()){
		getline(myfile,line_str);
		if(line_str.compare("object")==0){
			getline(myfile,line_str);
				int n = atoi(line_str.c_str());
				for(int i=0;i<n;i++){
		        	vector<vector3D> object_info ;
		        	objects object_single ;
					for(int ij=0;ij<18;ij++){
					if(ij<2){
					getline(myfile,line_str);
					string s = line_str;
		        	string delimiter = ",";
		        	size_t pos = 0;
		      	  	string token;
		        	vector<float> v;
		        	float co_ord  ;
		        while ((pos = s.find(delimiter)) != std::string::npos) {
		        	token = s.substr(0, pos);
					co_ord = atof(token.c_str()) ;
					v.push_back(co_ord);
		        	s.erase(0, pos + delimiter.length());
		    	}

		    	co_ord = atof(s.c_str()) ;
				v.push_back(co_ord);
				vector3D entry(v[0],v[1],v[2]);
				v.clear();
				if(ij==0){
					object_single.center = entry ;
				}
				else if(ij==1){
					object_single.normal = entry ;
				}
				}
				else if(ij==2){
					getline(myfile,line_str) ;
					object_single.radius = atof(line_str.c_str()) ;

					}
				else if(ij==3){
						getline(myfile,line_str) ;
						object_single.id = atoi(line_str.c_str()) ;

				}
				else if(ij==4){
										getline(myfile,line_str) ;
										object_single.D = atof(line_str.c_str()) ;

								}
				else if(ij==5){
										getline(myfile,line_str) ;
										object_single.kaR = atof(line_str.c_str());

								}
				else if(ij==6){
														getline(myfile,line_str) ;
														object_single.kaG = atof(line_str.c_str());

								}
				else if(ij==7){
														getline(myfile,line_str) ;
														object_single.kaB = atof(line_str.c_str());

												}
				else if(ij==8){
														getline(myfile,line_str) ;
														object_single.kdR = atof(line_str.c_str());

												}
				else if(ij==9){
														getline(myfile,line_str) ;
														object_single.kdG = atof(line_str.c_str());

												}
				else if(ij==10){
														getline(myfile,line_str) ;
														object_single.kdB = atof(line_str.c_str());

												}
				else if(ij==11){
														getline(myfile,line_str) ;
														object_single.ksR = atof(line_str.c_str());

												}
				else if(ij==12){
														getline(myfile,line_str) ;
														object_single.ksG = atof(line_str.c_str());

												}
				else if(ij==13){
														getline(myfile,line_str) ;
														object_single.ksB = atof(line_str.c_str());

												}
				else if(ij==14){
														getline(myfile,line_str) ;
														object_single.reflectionCoff = atof(line_str.c_str());

												}
				else if(ij==15){
														getline(myfile,line_str) ;
														object_single.refractionCoff = atof(line_str.c_str());

												}
				else if(ij==16){
														getline(myfile,line_str) ;
														object_single.type = line_str.c_str();

												}
				else if(ij==17){
														getline(myfile,line_str) ;
														object_single.n = atof(line_str.c_str());

												}

					}
				object.push_back(object_single) ;

				}
				}
		else if((line_str.compare("light")==0)){
			getline(myfile,line_str) ;
			int light_no = atoi(line_str.c_str());
			for(int li=0;li<light_no;li++){
				light light_sin ;
				for(int ly=0;ly<2;ly++){
				getline(myfile,line_str) ;
				string s = line_str;
				string delimiter = ",";
				size_t pos = 0;
				string token;
				vector<float> v;
				float co_ord  ;
				while ((pos = s.find(delimiter)) != std::string::npos) {
					token = s.substr(0, pos);

					co_ord = atof(token.c_str()) ;
					v.push_back(co_ord);
					s.erase(0, pos + delimiter.length());
				}

				co_ord = atof(s.c_str()) ;

				v.push_back(co_ord);
				v.clear();
				vector3D entry(v[0],v[1],v[2]);
				if(ly==0){

					light_sin.position = entry ;
				}
				else{

					light_sin.intensity = entry ;
				}
				}
				light_sources.push_back(light_sin) ;


			}

		}
		else{

				size_t found = line_str.find(":") ;
				string s = line_str.substr(found+1,line_str.length());
			    string delimiter = ",";
			    size_t pos = 0;
			    string token;
			    vector<float> v;
			    float co_ord  ;
			    while ((pos = s.find(delimiter)) != std::string::npos) {
			    	token = s.substr(0, pos);
					co_ord = atof(token.c_str()) ;
					v.push_back(co_ord);
			    	s.erase(0, pos + delimiter.length());
				}
			    co_ord = atof(s.c_str()) ;

				v.push_back(co_ord);
				v.clear();
				vector3D entry(v[0],v[1],v[2]);
				variables.push_back(entry);


	}
		}
		cout << variables.size() << endl ;
		cout << object[0].id << endl ;
		cout << object[0].type << endl ;
		myfile.close();
	}

	else{
	cout<< "not opende" << endl;
	}
}



point findIntersection(vector3D start,vector3D dir,vector<objects> items){
		float t=numeric_limits<float>::max();
		objects obj;
		bool interaction=false;
		bool plane;
		for(int i=0;i<items.size();i++){
			if(items[i].type=="sphere"){

				float A=1;
				float B=2*(dir.x*(start.x-items[i].center.x)+dir.y*(start.y-items[i].center.y)+dir.z*(start.z-items[i].center.z));
				float C=pow((start.x-items[i].center.x),2)+pow((start.y-items[i].center.y),2)+pow((start.z-items[i].center.z),2)-pow(items[i].radius,2);

				if((B*B-4*A*C)>=0){


					float t0= (-B-sqrt(B*B-4*A*C))/2.0*A;
					float t1= (-B+sqrt(B*B-4*A*C))/2.0*A;
					float tempT;
							if(t0>0.0 && t1>0.0){
								if(t0>0.0001 && t1>0.001){
									tempT=min(t0,t1);
								}
								else if(t0<0.0001 && t1>0.001){
									tempT=t1;
								}
								else if(t0>0.0001 && t1<0.001){
									tempT=t0;
								}
								else{
									tempT=numeric_limits<float>::max();
								}

							}
							else if(t0>0.0001 && t1<0){
								tempT=t0;
							}
							else if(t0<0 && t1>0.0001){
								tempT=t1;
							}
							else{
								tempT=numeric_limits<float>::max();
							}
							if(t>tempT){
								interaction=true;
								plane=false;
								t=tempT;
								obj=items[i];
							}
				}
			}
			else if(items[i].type=="plane"){
				float vd=dotProduct(dir,items[i].normal);
					if(vd<0){
						float tempT=-1.0*(dotProduct(items[i].normal,start)+items[i].D)/vd;
						if(t>tempT && tempT>0.001){
							interaction=true;
							plane=true;
							t=tempT;
							obj=items[i];
						}
					}
			}

		}
		point p;
		p.obj=obj;
		p.t=t;
		if(interaction){
			p.intersected=true;
		}
		else{
			p.intersected=false;
		}

		return p;

}




vector3D findColor(vector3D start,vector3D dir,float nMedium, vector<objects> items,vector<light> sun,vector3D amb, int level){
	if(level>4){
		vector3D newColor(0.529,0.80,0.921);
		return newColor;
	}
	else{
		point p=findIntersection(start,dir,items);
		if(!(p.intersected)){
			vector3D newColor(0.529,0.80,0.921);
			return newColor;
		}
		else{
		vector3D intersection(start.x+dir.x*p.t,start.y+dir.y*p.t,start.z+dir.z*p.t);
		vector3D intersectionNormal;
		if(p.obj.type=="plane"){
			intersectionNormal=p.obj.normal;
			double mult=1.2;
			if(((fmod(intersection.x*mult,2.0)<1.0 && fmod(intersection.x*mult,2.0)>=0.0) || (fmod(intersection.x*mult,2.0)>-1.0 && fmod(intersection.x*mult,2.0)<=0.0))&&((fmod(intersection.z*mult,2.0)<1.0 && fmod(intersection.z*mult,2.0)>=0.0) || (fmod(intersection.z*mult,2.0)>-1.0 && fmod(intersection.z*mult,2.0)<=0.0))){
				p.obj.kdR=0.00/255.0;
				p.obj.kdG=102.0/255.0;
				p.obj.kdB=102.0/255.0;

				p.obj.kaR=0.00/255.0;
				p.obj.kaG=102.0/255.0;
				p.obj.kaB=102.0/255.0;
			}
			else if(((fmod(intersection.x*mult,2.0)>1.0) || (fmod(intersection.x*mult,2.0)<-1.0))&&((fmod(intersection.z*mult,2.0)>1.0) || (fmod(intersection.z*mult,2.0)<-1.0))){
				p.obj.kdR=0.0;
				p.obj.kdG=102.0/255.0;
				p.obj.kdB=102.0/255.0;

				p.obj.kaR=0.0/255.0;
				p.obj.kaG=102.0/255.0;
				p.obj.kaB=102.0/255.0;
				}
			else{
				p.obj.kdR=1.0;
				p.obj.kdG=1.0;
				p.obj.kdB=1.0;

				p.obj.kaR=1.0;
				p.obj.kaG=1.0;
				p.obj.kaB=1.0;
			}

		}
		else{
			vector3D intersectNorm((intersection.x-p.obj.center.x)/p.obj.radius,(intersection.y-p.obj.center.y)/p.obj.radius,(intersection.z-p.obj.center.z)/p.obj.radius);
			intersectionNormal=normalise(intersectNorm);
		}
		float cosT = angle(constProduct(dir,-1),intersectionNormal);
		vector3D colorReflection(0.0,0.0,0.0);
		vector3D colorRefraction(0.0,0.0,0.0);
		float kReflection=0.25;
		bool reflection=false;
		bool refraction=false;
		if(p.obj.reflectionCoff>0.01){
			reflection=true;
		vector3D nextReflection=subtract(constProduct(intersectionNormal,2*cosT),normalise(constProduct(dir,-1)));

		colorReflection=constProduct(findColor(intersection,nextReflection,nMedium,items,sun,amb,level+1),2);
		kReflection=0.25;
		}
		else{
			vector3D color1(0.0,0.0,0.0);
			colorReflection=color1;
		}
		if(p.obj.refractionCoff>0.01){
				refraction=true;
				if(nMedium==1.0){
					float sinThetaI=sqrt(1.0-pow(cosT,2.0));
					float sinThetaT=sinThetaI*1.0/p.obj.refractionCoff;
					float cosThetaT=sqrt(1.0-pow(sinThetaT,2.0));

					vector3D nextRefraction= subtract(constProduct(add(dir,constProduct(intersectionNormal,cosT)),(1.0/p.obj.refractionCoff)),constProduct(intersectionNormal,cosThetaT));
					colorRefraction=constProduct(findColor(intersection,nextRefraction,p.obj.refractionCoff,items,sun,amb,level+1),2);

				}
				else{
					float sinThetaI=sqrt(1.0-pow(abs(cosT),2.0));
					float sinThetaT=sinThetaI*p.obj.refractionCoff;
					float cosThetaT=sqrt(1.0-pow(sinThetaT,2.0));
					vector3D nextRefraction= subtract(constProduct(add(dir,constProduct(intersectionNormal,-1.0*abs(cosT))),(p.obj.refractionCoff)),constProduct(intersectionNormal,-1.0*cosThetaT));
					colorRefraction=constProduct(findColor(intersection,nextRefraction,1.0,items,sun,amb,level+1),2);

				}

			}
			else{
					vector3D color1(0.0,0.0,0.0);
					colorRefraction=color1;
			}
		if(reflection && refraction){
			kReflection=0.25;
		}
		else if(reflection && !refraction){
			kReflection=1.0;
		}
		else if(refraction && !reflection){
			kReflection=0.0;
		}
		else{
			kReflection=0.5;
		}

		float intensityR=p.obj.kaR*amb.x;
		float intensityG=p.obj.kaG*amb.y;
		float intensityB=p.obj.kaB*amb.z;


		intensityR=intensityR+kReflection*p.obj.kdR*colorReflection.x*cosT+kReflection*p.obj.ksR*colorReflection.x+(1.0-kReflection)*p.obj.kdR*colorRefraction.x*cosT+(1.0-kReflection)*p.obj.ksR*colorRefraction.x;
		intensityG=intensityG+kReflection*p.obj.kdG*colorReflection.y*cosT+kReflection*p.obj.ksG*colorReflection.y+(1.0-kReflection)*p.obj.kdG*colorRefraction.y*cosT+(1.0-kReflection)*p.obj.ksG*colorRefraction.y;
		intensityB=intensityB+kReflection*p.obj.kdB*colorReflection.z*cosT+kReflection*p.obj.ksB*colorReflection.z+(1.0-kReflection)*p.obj.kdB*colorRefraction.z*cosT+(1.0-kReflection)*p.obj.ksB*colorRefraction.z;

		for(int i=0;i<sun.size();i++){
			vector3D lightDir=subtract(sun[i].position,intersection);
			vector3D lightDirection=normalise(lightDir);

			float cosTheta = angle(lightDirection,intersectionNormal);

			vector3D reflectionDir=subtract(constProduct(intersectionNormal,2*cosTheta),normalise(lightDir));
//			float cosAlpha=angle(intersectionNormal,normalise(add(lightDir,constProduct(dir,-1))));
			float cosAlpha = angle(constProduct(dir,-1),reflectionDir);
			if(cosAlpha<=0){
					cosAlpha=0;
			}
			if(cosTheta<=0){
				cosTheta=0;
				cosAlpha=0;
			}
			bool obstacle=false;
			for(int j=0;j<items.size();j++){
				if(p.obj.id!=items[j].id){
				if(items[j].type=="sphere"){
					float A=1;
					float B=2*(lightDirection.x*(intersection.x-items[j].center.x)+lightDirection.y*(intersection.y-items[j].center.y)+lightDirection.z*(intersection.z-items[j].center.z));
					float C=pow((intersection.x-items[j].center.x),2)+pow((intersection.y-items[j].center.y),2)+pow((intersection.z-items[j].center.z),2)-pow(items[j].radius,2);
					if((B*B-4*A*C)>=0){

						float t0= (-B-sqrt(B*B-4*A*C))/2.0*A;
						float t1= (-B+sqrt(B*B-4*A*C))/2.0*A;

						if(t0>=0 || t1>=0){
							obstacle=true;
							break;
						}

					}
				}
				}
			}

			if(obstacle){
				cosTheta=0;
				cosAlpha=0;
			}
			intensityR=intensityR+p.obj.kdR*sun[i].intensity.x*cosTheta+p.obj.ksR*sun[i].intensity.x*(pow(cosAlpha,p.obj.n));
			intensityG=intensityG+p.obj.kdG*sun[i].intensity.y*cosTheta+p.obj.ksG*sun[i].intensity.y*(pow(cosAlpha,p.obj.n));
			intensityB=intensityB+p.obj.kdB*sun[i].intensity.z*cosTheta+p.obj.ksB*sun[i].intensity.z*(pow(cosAlpha,p.obj.n));


		}

		vector3D newColor(intensityR,intensityG,intensityB);
		return newColor;
	  }
	}
}
void raytracer(vector3D dir,vector<objects> items,vector<light> sun,vector3D amb){
	vector3D wEye=convert(vcsEye);
	point p=findIntersection(wEye,dir,items);
	if((p.intersected)){

	vector3D intersection(wEye.x+dir.x*p.t,wEye.y+dir.y*p.t,wEye.z+dir.z*p.t);
	vector3D color=findColor(wEye,dir,1.0,items,sun,amb,0);
	glBegin(GL_POINTS );
					glColor3f(color.x,color.y,color.z);
					vector3D st=intersection;


					      			GLfloat start[4][4];

					      			   start[0][0]=st.x;
					      			   start[0][1]=st.y;
					      			   start[0][2]=st.z;
					      			   start[0][3]=1;


					      			GLfloat startR[4][4];

					      			matrixMult(start,wcsTvcs,startR,1,4,4);


					      			matrixMult(startR,perspect,start,1,4,4);

					      			GLfloat divide1=start[0][3];


					      			for(int j=0;j<4;j++){
					      				start[0][j]=start[0][j]/divide1;

					      			}

					      			matrixMult(start,inv,startR,1,4,4);

					      			 glVertex3f(startR[0][0],startR[0][1],startR[0][2]);

				glEnd();
	}

}



void generate(void)
{
     glMatrixMode(GL_MODELVIEW);
     glClearColor (0.529,0.80,0.921,0);
     glClear(GL_COLOR_BUFFER_BIT);
     glLoadIdentity();
     vector3D myeye=convert(vcsEye);
     printVector(myeye);
     float s=-3;
     gluLookAt(origin.x+n.x*s,origin.y,origin.z+n.z*s,origin.x,origin.y,origin.z,0.0,1.0,0.0);

	for(float x=-1.0;x<1.0;x=x+0.005){
		for(float y=1.0;y>-1.0;y=y-0.005){
			vector3D vcsD(x,y,3.0);
			raytracer(normalise(convert(vcsD)),object,light_sources,ambient);
		}
	}

	glFlush();
}






void reshape(int x, int y)
{
    if (y == 0 || x == 0) return;  //Nothing is visible then, so return
    glMatrixMode(GL_PROJECTION);
    glLoadIdentity();
   //    gluPerspective(90.0,(GLdouble)x/(GLdouble)y,0.5,200.0);
    glOrtho(-1.0,1.0,-1.0,1.0,-1000.0,1000.0);
    glMatrixMode(GL_MODELVIEW);
    glViewport(0,0,x,y);  //Use the whole window for rendering

}




int main(int argc, char** argv){


glutInit(&argc, argv);
read();
normal=variables[0];
printVector(normal);
xyzUp=variables[1];
printVector(xyzUp);
n=normalise(normal);
v= normalise(subtract(xyzUp,constProduct(n,dotProduct(xyzUp,n))));
u=crossProduct(v,n);
M[0][0]=u.x;
M[0][1]=u.y;
M[0][2]=u.z;
M[1][0]=v.x;
M[1][1]=v.y;
M[1][2]=v.z;
M[2][0]=n.x;
M[2][1]=n.y;
M[2][2]=n.z;

origin=variables[2];
vcsEye=variables[3];
nearPlane=variables[4];
farPlane=variables[5];

viewC1=variables[6];
viewC2=variables[7];
viewC3=variables[8];
viewC4=variables[9];
ambient=variables[10];

viewPlane[0]=viewC1;
viewPlane[1]=viewC2;
viewPlane[2]=viewC3;
viewPlane[3]=viewC4;

matrixOrigin[0][0]=origin.x;
matrixOrigin[0][1]=origin.y;
matrixOrigin[0][2]=origin.z;

xyzTuvn[0]=u.x;
xyzTuvn[1]=v.x;
xyzTuvn[2]=n.x;
xyzTuvn[3]=0;
xyzTuvn[4]=u.y;
xyzTuvn[5]=v.y;
xyzTuvn[6]=n.y;
xyzTuvn[7]=0;
xyzTuvn[8]=u.z;
xyzTuvn[9]=v.z;
xyzTuvn[10]=n.z;
xyzTuvn[11]=0;
xyzTuvn[12]= (-1)*dotProduct(origin,u);
xyzTuvn[13]=(-1)*dotProduct(origin,v);
xyzTuvn[14]=(-1)*dotProduct(origin,n);
xyzTuvn[15]=1;

int row=-1;
int col=0;
for(int i=0;i<16;i++){
	if(i%4==0){
		row++;
		col=0;
	}
	wcsTvcs[row][col]=xyzTuvn[i];
	col++;
}


for(int i=0;i<16;i++){
	if(i%5==0){
		perspective[i]=1;
	}
	else{
		perspective[i]=0;
	}

}
perspective[11]=-1.0/(float)(vcsEye.z);
perspective[10]=0;

 row=-1;
 col=0;
for(int i=0;i<16;i++){
	if(i%4==0){
		row++;
		col=0;
	}
	perspect[row][col]=perspective[i];
	col++;
}

	GLfloat out[16];
	gluInvertMatrix(xyzTuvn,out);

	 row=-1;
	 col=0;
	for(int i=0;i<16;i++){
		if(i%4==0){
			row++;
			col=0;
		}
		inv[row][col]=out[i];
		col++;
	}


glutInitDisplayMode(GLUT_SINGLE|GLUT_RGB);
glutInitWindowPosition(100, 100);
glutCreateWindow(argv[0]);
init();
glutDisplayFunc(generate);
glutReshapeFunc(reshape);
glutMainLoop();
return 0;
}




