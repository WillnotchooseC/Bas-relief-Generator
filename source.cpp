#define _CRT_SECURE_NO_DEPRECATE
#include <stdio.h>
#include <vector>
#include <cstdlib>
#include <algorithm>
#include <cmath> 
#include <string.h>
#include <stdarg.h>
#include <iostream>
#include <time.h>
#include <math.h>
#include <cstdio> 
#include <cmath> 
#include <fstream> 
#include <vector> 
#include <cassert>
#include <iostream>
#include <Eigen/Dense>
#include <Eigen/SparseCore>
#include <Eigen/Sparse>
#include <ctime> 
#include <stdlib.h> 
#include <cassert> 
#include <fstream>
#include <windows.h>
#include <gdiplus.h>
#include <Eigen/Sparse>
#define M_PI 3.1415
using namespace Eigen;
using namespace std;
using Eigen::MatrixXd;
typedef Eigen::SparseMatrix<double> SpMat;
typedef Eigen::Triplet<double> T;
using namespace std;
double static u = 0.5;//when u going up, height going down.
double static theta = 0.5;//threshold for all heights
double static canshu = 1.0;//coefficient of all LN
#define BUF_SIZE 256

#define INFINITE 10000
#define min(x,y) ((x) < (y) ? (x) : (y))
#define max(x,y) ((x) < (y) ? (y) : (x))
double static **red;
double static **blue;
double static **green;
double static **fred;
double static **fblue;
double static **fgreen;
double static **ln;
double static **alpha;
double static weight[5][5];
double static xinterval;
double static yinterval;
double static original;
double static minln;
static int height = 500;
static int width = 400;

//bool isPoint = false, isLine = false, isFill = false, isBoth = true, isOn = true, isFlat = true;
struct NORMAL {
	double x;
	double y;
	double z;
	void normalize();
};

void NORMAL::normalize()
{
	double factor = sqrt(pow((double)this->x, 2) + pow((double)this->y, 2) + pow((double)this->z, 2));
	this->x /= factor;
	this->y /= factor;
	this->z /= factor;
}



struct vector2{
	double x;
	double y;
	double z;

};

struct INT3VECT {
	int a;
	int b;
	int c;
	bool canuse;
};

struct FLTVECT
{
	double x;
	double y;
	double z;
};

struct COLORVECT {
	int red;
	int green;
	int blue;
	int alpha;

};

struct PIXEL {
	double u, v;
	bool checked;
	int facen;
	double deep;
	NORMAL *normalvalue;
};

struct SurFacemesh {
	int nv;
	int nf;
	double minx, miny, maxx, maxy, maxz;
	FLTVECT *vertex;
	INT3VECT *face;
	//COLORVECT *color;
	NORMAL *normal;
	NORMAL *vertexNormal;
};

SurFacemesh* readPolygon()
{
	int num, n, m;
	int a, b, c, d;
	double x, y, z;
	int red, green, blue, alpha;
	SurFacemesh *surfmesh;
	char line[256];
	FILE *fin;

	char* fileName = "elephant.off";

	if ((fin = fopen(fileName, "r")) == NULL) {
		printf("1 read error...\n");
		exit(0);
	};

	/* OFF format */
	while (fgets(line, 256, fin) != NULL) {
		if (line[0] == 'O' && line[1] == 'F' && line[2] == 'F')
			break;
	}
	fscanf(fin, "%d %d %d\n", &m, &n, &num);

	// allocate size for surface mesh
	surfmesh = (SurFacemesh*)malloc(sizeof(SurFacemesh));
	surfmesh->nv = m;
	surfmesh->nf = n;
	surfmesh->vertex = (FLTVECT *)malloc(sizeof(FLTVECT)*surfmesh->nv);
	surfmesh->vertexNormal = (NORMAL *)malloc(sizeof(NORMAL)*surfmesh->nv);
	surfmesh->face = (INT3VECT *)malloc(sizeof(INT3VECT)*surfmesh->nf);
	surfmesh->normal = (NORMAL *)malloc(sizeof(NORMAL)*surfmesh->nf);
	//surfmesh->color = (COLORVECT *)malloc(sizeof(COLORVECT)*surfmesh->nv);

	// store vertices
	for (n = 0; n < surfmesh->nv; n++) {
		fscanf(fin, "%lf %lf %lf \n", &x, &y, &z);
		surfmesh->vertex[n].x = (double)x*1000;
		surfmesh->vertex[n].y = (double)y*1000;
		surfmesh->vertex[n].z = (double)z*1000;


		// initialize vertex normal variables
		surfmesh->vertexNormal[n].x = surfmesh->vertexNormal[n].y = surfmesh->vertexNormal[n].z = 0;
	}
	//printf("%d\n", surfmesh->color[265].blue);
	// store faces
	for (n = 0; n < surfmesh->nf; n++) {
		fscanf(fin, "%d %d %d %d\n", &a, &b, &c, &d);
		surfmesh->face[n].a = b;
		surfmesh->face[n].b = c;
		surfmesh->face[n].c = d;
		//printf("%d %d %d \n", surfmesh->face[n].a, surfmesh->face[n].b, surfmesh->face[n].c);
		if (a != 3)
			printf("Errors: reading surfmesh .... \n");
	}



	for (n = 0; n < surfmesh->nf; n++)
	{
		vector2 *v1 = new vector2, *v2 = new vector2;

		v1->x = (double)surfmesh->vertex[surfmesh->face[n].b].x - surfmesh->vertex[surfmesh->face[n].a].x;
		v1->y = (double)surfmesh->vertex[surfmesh->face[n].b].y - surfmesh->vertex[surfmesh->face[n].a].y;
		v1->z = (double)surfmesh->vertex[surfmesh->face[n].b].z - surfmesh->vertex[surfmesh->face[n].a].z;

		v2->x = (double)surfmesh->vertex[surfmesh->face[n].c].x - surfmesh->vertex[surfmesh->face[n].a].x;
		v2->y = (double)surfmesh->vertex[surfmesh->face[n].c].y - surfmesh->vertex[surfmesh->face[n].a].y;
		v2->z = (double)surfmesh->vertex[surfmesh->face[n].c].z - surfmesh->vertex[surfmesh->face[n].a].z;

		// calculate face normals
		double normx = (v1->y * v2->z) - (v2->y * v1->z),
			normy = (v1->z * v2->x) - (v2->z * v1->x),
			normz = (v1->x * v2->y) - (v2->x * v1->y);

		surfmesh->normal[n].x = normx;
		surfmesh->normal[n].y = normy;
		surfmesh->normal[n].z = normz;
		surfmesh->normal[n].normalize();

		int verta = surfmesh->face[n].a,
			vertb = surfmesh->face[n].b,
			vertc = surfmesh->face[n].c;

		// for each face add its normal to the vertices that define it
		surfmesh->vertexNormal[verta].x += normx;
		surfmesh->vertexNormal[vertb].x += normx;
		surfmesh->vertexNormal[vertc].x += normx;
		surfmesh->vertexNormal[verta].y += normy;
		surfmesh->vertexNormal[vertb].y += normy;
		surfmesh->vertexNormal[vertc].y += normy;
		surfmesh->vertexNormal[verta].z += normz;
		surfmesh->vertexNormal[vertb].z += normz;
		surfmesh->vertexNormal[vertc].z += normz;

		delete[] v1, v2;
	}

	for (int i = 0; i < surfmesh->nv; i++)
	{
		surfmesh->vertexNormal[i].normalize();
		
	}
	surfmesh->minx = surfmesh->miny = INFINITE;
	surfmesh->maxx = surfmesh->maxy = surfmesh->maxz = -INFINITE;
	for (int i = 0; i < surfmesh->nv; i++) {
		if (surfmesh->vertex[i].x < surfmesh->minx) {
			surfmesh->minx = surfmesh->vertex[i].x;
		}
		if (surfmesh->vertex[i].y < surfmesh->miny) {
			surfmesh->miny = surfmesh->vertex[i].y;
		}
		if (surfmesh->vertex[i].x > surfmesh->maxx) {
			surfmesh->maxx = surfmesh->vertex[i].x;
		}
		if (surfmesh->vertex[i].y > surfmesh->maxy) {
			surfmesh->maxy = surfmesh->vertex[i].y;
		}
		if (surfmesh->vertex[i].z > surfmesh->maxz) {
			surfmesh->maxz = surfmesh->vertex[i].z;
		}
	}

	fclose(fin);

	return surfmesh;
}
//this function judge whether the point is in the face or not
bool currentfaceinclude(double curx, double cury, double ax, double ay, double bx, double by, double cx, double cy) {
	double cax = curx - ax;
	double cay = cury - ay;
	double cbx = curx - bx;
	double cby = cury - by;
	double ccx = curx - cx;
	double ccy = cury - cy;
	return (((cax*cby) - (cay*cbx))*((cbx*ccy) - (ccx*cby)) >= 0 &&
		((cbx*ccy) - (ccx*cby))*((ccx*cay) - (ccy*cax)) >= 0);
}
//Gaussian filter template construction
void gaussianIni(double weight[][5],int sigma) {
	double r, s = 2.0 * sigma * sigma;
	double sum = 0.0;

	for (int x = -2; x <= 2; x++)
	{
		for (int y = -2; y <= 2; y++)
		{
			r = sqrt(x*x + y*y);
			weight[x + 2][y + 2] =
				(exp(-(r*r) / s)) / (sqrt(2*M_PI) * sigma);
			sum += weight[x + 2][y + 2];
		}
	}

	// normalising the Kernel
	for (int i = 0; i < 5; ++i)
		for (int j = 0; j < 5; ++j)
			weight[i][j] /= sum;
}
//gaussain Filter the red color, namely the x value from normal
void gaussianFilterred() {
	fred = new double*[height];
	/*double weight2[5][5] = {
		{ 0.11, 0.11, 0.11 } ,
		{ 0.11, 0.11, 0.11 } ,
		{ 0.11, 0.11, 0.11 }
	};*/
	for (int i = 0; i < height; i++) {
		fred[i] = new double[width];
		for (int j = 0; j < width; j++) {
			if (alpha[i][j] > 0.0&&i>1&&j>1&&i<height-2&&j<width-2) {
				fred[i][j] = 0.0;
				for (int m = 0; m < 5; m++) {
					for (int n = 0; n < 5; n++) {
						fred[i][j] += red[i - 2 + m][j - 2 + n] * weight[m][n];
					}
				}
				//red[i][j] + 1.0) / 2.0 * 255
				//fred[i][j] = fred[i][j] / 255 * 2.0 - 1.0;
			}
			else {
				fred[i][j] = 0.0;
			}
		}
	}
}
//gaussian filter the green color
void gaussianFiltergreen() {
	fgreen = new double*[height];
	/*double weight2[5][5] = {
		{ 0.11, 0.11, 0.11 } ,
		{ 0.11, 0.11, 0.11 } ,
		{ 0.11, 0.11, 0.11 }
	};*/
	for (int i = 0; i < height; i++) {
		fgreen[i] = new double[width];
		for (int j = 0; j < width; j++) {
			if (alpha[i][j] > 0.0&&i>1 && j>1 && i<height - 2 && j<width - 2) {
				fgreen[i][j] = 0.0;
				for (int m = 0; m < 5; m++) {
					for (int n = 0; n < 5; n++) {
						fgreen[i][j] += green[i - 2 + m][j - 2 + n] * weight[m][n];
					}
				}
				//fgreen[i][j] = fgreen[i][j] / 255 * 2.0 - 1.0;
			}
			else {
				fgreen[i][j] = 0.0;
			}
		}
	}

	
}

void gaussianFilterblue() {
	fblue = new double*[height];
	/*double weight2[5][5] = {
	{ 0.11, 0.11, 0.11 } ,
	{ 0.11, 0.11, 0.11 } ,
	{ 0.11, 0.11, 0.11 }
	};*/
	for (int i = 0; i < height; i++) {
		fblue[i] = new double[width];
		for (int j = 0; j < width; j++) {
			if (alpha[i][j] > 0.0&&i>1 && j>1 && i<height - 2 && j<width - 2) {
				fblue[i][j] = 0.0;
				/*for (int m = 0; m < 5; m++) {
					for (int n = 0; n < 5; n++) {
						fblue[i][j] += blue[i - 2 + m][j - 2 + n] * weight[m][n];
					}
				}*/
				//fblue[i][j] = blue[i][j] / 255 * 2.0 - 1.0;
			}
			else {
				fblue[i][j] = 1.0;
			}
		}
	}


}
//calculate the distance from a point to a line
double distanceLine(double linesx, double linesy, double linedx, double linedy, double u, double v) {
	int linex = linedx - linesx;
	int liney = linedy - linesy;

	int psx = linesx - u;
	int psy = linesy - v;
	int dot = psx * linex + psy* liney;

	int distancex = psx - dot * linex;
	int distancey = psy - dot * liney;
	double output = sqrt((double)(distancex * distancex + distancey * distancey));
	return output;
}
//filter out the maximum value among three number
double maxoftriple(double x, double y, double z) {
	if (x > y) {
		if (x > z) {
			return x;
		}
		else {
			return z;
		}
	}
	else {
		if (y > z) {
			return y;
		}
		else {
			return z;
		}
	}
}

double minoftriple(double x, double y, double z) {
	if (x < y) {
		if (x < z) {
			return x;
		}
		else {
			return z;
		}
	}
	else {
		if (y < z) {
			return y;
		}
		else {
			return z;
		}
	}
}

PIXEL* makeNormalMap(SurFacemesh* meshoff) {
	/*double meshw = meshoff->maxx - meshoff->minx;
	double meshh = meshoff->maxy - meshoff->miny;
	double meshz = meshoff->maxz;
	xinterval = meshw / width;
	yinterval = meshh / height;*/
	double xminbound = INFINITE;
	double yminbound = INFINITE;
	double xmaxbound = -INFINITE;
	double ymaxbound = -INFINITE;
	double meshw = meshoff->maxx - meshoff->minx;
	double meshh = meshoff->maxy - meshoff->miny;
	double meshz = meshoff->maxz;
	meshw = meshw / 5.0*6.0;
	meshh = meshh / 5.0*7.0;
	meshoff->minx = meshoff->minx - meshw / 10.0;
	meshoff->maxy = meshoff->maxy + meshh / 5.0;
	xinterval = meshw / width;
	yinterval = meshh / height;
	//vector selectFaces = { 0,0,1 };
	cout << "xinterval is: " << xinterval << endl;
	cout << "yinterval is: " << yinterval << endl;
	PIXEL* pixels;
	pixels = (PIXEL *)malloc(sizeof(PIXEL)*height*width);
	pixels->normalvalue = (NORMAL *)malloc(sizeof(NORMAL)*height*width);
	for (int i = 0; i < height; i++) {
		for (int j = 0; j < width; j++) {
			pixels[i*width + j].u = meshoff->minx + xinterval*j;
			pixels[i*width + j].v = meshoff->maxy - yinterval*i;
		}
	}
	for (int i = 0; i < meshoff->nf; i++) {
		if (meshoff->normal[i].z < 0) {
			meshoff->face[i].canuse = false;
		}
		else {
			meshoff->face[i].canuse = true;
		}
	}
	for (int k = 0; k < height*width; k++) {
		pixels[k].checked = false;
		pixels[k].deep = -INFINITE;
	}
	double facez = 0.0;
	for (int i = 0; i < meshoff->nf; i++) {
		if (meshoff->face[i].canuse) {
			xmaxbound = maxoftriple(meshoff->vertex[meshoff->face[i].a].x, meshoff->vertex[meshoff->face[i].b].x, meshoff->vertex[meshoff->face[i].c].x);
			xminbound = minoftriple(meshoff->vertex[meshoff->face[i].a].x, meshoff->vertex[meshoff->face[i].b].x, meshoff->vertex[meshoff->face[i].c].x);
			ymaxbound = maxoftriple(meshoff->vertex[meshoff->face[i].a].y, meshoff->vertex[meshoff->face[i].b].y, meshoff->vertex[meshoff->face[i].c].y);
			yminbound = minoftriple(meshoff->vertex[meshoff->face[i].a].y, meshoff->vertex[meshoff->face[i].b].y, meshoff->vertex[meshoff->face[i].c].y);
			for (int j = 0; j < height*width; j++) {
				if (pixels[j].u > xminbound&&pixels[j].u < xmaxbound&&pixels[j].v<ymaxbound&&pixels[j].v>yminbound) {
					facez = (meshoff->vertex[meshoff->face[i].a].z + meshoff->vertex[meshoff->face[i].b].z + meshoff->vertex[meshoff->face[i].c].z) / 3.0;
					if (pixels[j].deep < facez) {
						//surfmesh->vertex[surfmesh->face[n].c].x
						if (currentfaceinclude(pixels[j].u, pixels[j].v, meshoff->vertex[meshoff->face[i].a].x,
							meshoff->vertex[meshoff->face[i].a].y, meshoff->vertex[meshoff->face[i].b].x,
							meshoff->vertex[meshoff->face[i].b].y, meshoff->vertex[meshoff->face[i].c].x,
							meshoff->vertex[meshoff->face[i].c].y)) {
							pixels[j].facen = i;
							pixels[j].checked = true;
							pixels[j].deep = facez;
						}
					}
				}
			}
		}
	}
	/*cout << "now let's output the normal of vertex." << endl;
	for (int i = 0; i < meshoff->nv; i++) {
		cout << meshoff->vertexNormal[i].x << endl;
		cout << meshoff->vertexNormal[i].y << endl;
	}
	*/
	double ax, ay, bx, by, cx, cy;
	double da, db, dc;
	double la, lb, lc;
	for (int i = 0; i < height; i++) {
		for (int j = 0; j < width; j++) {
			if (!pixels[i*width + j].checked) {
				pixels->normalvalue[i*width + j].x = 0;
				pixels->normalvalue[i*width + j].y = 0;
				pixels->normalvalue[i*width + j].z = 1;
			}
			//this is calcluate the point position in the triangle, more accurate, however, it will easily cause noise.
			/*else {
				ax = meshoff->vertex[meshoff->face[pixels[i*width + j].facen].a].x;
				ay = meshoff->vertex[meshoff->face[pixels[i*width + j].facen].a].y;
				bx = meshoff->vertex[meshoff->face[pixels[i*width + j].facen].b].x;
				by = meshoff->vertex[meshoff->face[pixels[i*width + j].facen].b].y;
				cx = meshoff->vertex[meshoff->face[pixels[i*width + j].facen].c].x;
				cy = meshoff->vertex[meshoff->face[pixels[i*width + j].facen].c].y;
				//float linesx, float linesy, float linedx, float linedy, float u, float v
				da = distanceLine(ax, ay, bx, by, pixels[i*width + j].u, pixels[i*width + j].v);
				db = distanceLine(bx, by, cx, cy, pixels[i*width + j].u, pixels[i*width + j].v);
				dc = distanceLine(cx, cy, ax, ay, pixels[i*width + j].u, pixels[i*width + j].v);
				la = sqrt((double)((ax - bx) * (ax - bx) + (ay - by) * (ay - by)));
				lb = sqrt((double)((bx - cx) * (bx - cx) + (by - cy) * (by - cy)));
				lc = sqrt((double)((cx - ax) * (cx - ax) + (cy - ay) * (cy - ay)));
				//surfmesh->vertexNormal[surfmesh->face[n].a].x
				pixels->normalvalue[i*width + j].x = (meshoff->vertexNormal[meshoff->face[pixels[i*width + j].facen].a].x *
					lb*db + meshoff->vertexNormal[meshoff->face[pixels[i*width + j].facen].b].x * lc * dc +
					meshoff->vertexNormal[meshoff->face[pixels[i*width + j].facen].c].x * la * da) / (lb*db + lc*dc + la*da);
				pixels->normalvalue[i*width + j].y = (meshoff->vertexNormal[meshoff->face[pixels[i*width + j].facen].a].y *
					lb*db + meshoff->vertexNormal[meshoff->face[pixels[i*width + j].facen].b].y * lc * dc +
					meshoff->vertexNormal[meshoff->face[pixels[i*width + j].facen].c].y * la * da) / (lb*db + lc*dc + la*da);
				pixels->normalvalue[i*width + j].z = (meshoff->vertexNormal[meshoff->face[pixels[i*width + j].facen].a].z *
					lb*db + meshoff->vertexNormal[meshoff->face[pixels[i*width + j].facen].b].z * lc * dc +
					meshoff->vertexNormal[meshoff->face[pixels[i*width + j].facen].c].z * la * da) / (lb*db + lc*dc + la*da);
				pixels->normalvalue[i*width + j].normalize();
			}*/
			//this is simple give triangle's normal to pixel, it is a little unaccurate, but it has less noise.
			else {

				pixels->normalvalue[i*width + j].x = meshoff->normal[pixels[i*width + j].facen].x;
				pixels->normalvalue[i*width + j].y = meshoff->normal[pixels[i*width + j].facen].y;
				pixels->normalvalue[i*width + j].z = meshoff->normal[pixels[i*width + j].facen].z;
				pixels->normalvalue[i*width + j].normalize();
			}
		}
		

	}


	return pixels;
}

void treatOriginalData(PIXEL* pixels) {
	red = new double*[height];
	green = new double*[height];
	blue = new double*[height];
	alpha = new double*[height];
	for (int i = 0; i < height; i++) {
		red[i] = new double[width];
		green[i] = new double[width];
		blue[i] = new double[width];
		alpha[i] = new double[width];
		for (int j = 0; j < width; j++) {
			if (pixels[i*width + j].checked) {
				//red[i][j]+1.0)/2.0 * 255
				red[i][j] = pixels->normalvalue[i*width + j].x ;
				green[i][j] = pixels->normalvalue[i*width + j].y;
				blue[i][j] = pixels->normalvalue[i*width + j].z ;
				alpha[i][j] = 1.0;
			}
			else {
				alpha[i][j] = 0.0;
				red[i][j] = 0.0;
				green[i][j] = 0.0;
				blue[i][j] = 1.0;
			}
			/*if (!(red[i][j] == 0 && green[i][j] == 0 && blue[i][j] == 1.0)) {
				alpha[i][j] = 1.0;
			}
			else {
				alpha[i][j] = 0.0;
			}*/
		}
	}
	free(pixels);
}
/*
double median(double *p) {
	int size = sizeof(p) / sizeof(double);
	sort(p, p + size);
	if (size % 2 == 0) {
		return (p[size / 2] + p[size / 2 + 1]) / 2.0;
	}
	else {
		return p[size / 2];
	}
}

void medianfilter(double **color,double **rcolor) {
	int radius = 1.0;
	rcolor = new double*[height];
	for (int i = 0; i < height; i++) {
		int top = max(i - radius, 0);
		int bottom = min(i + radius, height - 1);
		rcolor[i] = new double[width];
		for (int j = 0; j < width; j++) {
			int left = max(j - radius, 0);
			int right = min(j + radius, width - 1);
			int size = (bottom - top + 1)*(right-left+1);
			double* p = new double[size];
			int i = 0;
			for (int v = top; v < bottom + 1; v++) {
				for (int u = left; u < right + 1; u++) {
					p[i++] = color[v][u];
				}
			}
			rcolor[i][j] = median(p);

		}
	}
}*/

void initDataTreat(double **red, double **green, double **blue) {
	gaussianIni(weight, 1.5);
	gaussianFilterred();
	gaussianFiltergreen();
	//gaussianFilterblue();
	ln = new double*[height];
	double **xtemp = new double*[height];//right - left x 
	double **ytemp = new double*[height];//top - bottom y 
	for (int i = 0; i < height; i++) {
		xtemp[i] = new double[width];
		for (int j = 0; j < width; j++) {
			if (j == 0) {
				xtemp[i][j] = (red[i][j] / blue[i][j] - red[i][j+1] / blue[i][j+1]);
			}
			else if (j == width - 1) {
				xtemp[i][j] = (red[i][j-1] / blue[i][j-1] - red[i][j] / blue[i][j]);
			}
			else {
				xtemp[i][j] = (red[i][j-1] / blue[i][j-1] - red[i][j+1] / blue[i][j+1]) / (2.0);
			}
		}
	}
	for (int i = 0; i < height; i++) {
		ytemp[i] = new double[width];
		for (int j = 0; j < width; j++) {
			if (i == 0) {
				ytemp[i][j] = (green[i+1][j] / blue[i+1][j] - green[i][j] / blue[i][j]);
			}
			else if (i == height - 1) {
				ytemp[i][j] = (green[i][j] / blue[i][j] - green[i-1][j] / blue[i-1][j]);
			}
			else {
				ytemp[i][j] = (green[i + 1][j] / blue[i + 1][j] - green[i - 1][j] / blue[i - 1][j]) / (2.0);
			}
		}
	}
	/*for (int i = 0; i < height; i++) {
		xtemp[i] = new double[width];
		ytemp[i] = new double[width];
		for (int j = 0; j < width; j++) {
			
			xtemp[i][j] = red[i][j] / blue[i][j];
			ytemp[i][j] = green[i][j] / blue[i][j];
		}
	}*/
	/*gaussianFilter(red, 0.1);
	gaussianFilter(green, 0.1);
	gaussianFilter(blue, 0.1);*/
	/*for (int i = 0; i < height; i++) {
		ln[i] = new double[width];
		for (int j = 0; j < width; j++) {
			if (i == 0 && j == 0) {
				ln[i][j] = (-green[i + 1][j] / blue[i + 1][j]) / (2) + (red[i][j + 1] / blue[i][j + 1]) / (2);
				continue;
			}
			else if (i == height - 1 && j == 0) {
				ln[i][j] = (green[i - 1][j] / blue[i - 1][j]) / (2) + (red[i][j + 1] / blue[i][j + 1]) / (2);
				continue;
			}
			else if (i == 0 && j == width - 1) {
				ln[i][j] = (-green[i + 1][j] / blue[i + 1][j]) / (2) + (-red[i][j - 1] / blue[i][j - 1]) / (2);
				continue;
			}
			else if (i == height - 1 && j == width - 1) {
				ln[i][j] = (green[i - 1][j] / blue[i - 1][j]) / (2) + (-red[i][j - 1] / blue[i][j - 1]) / (2);
				continue;
			}
			else if (i == 0) {

				ln[i][j] = (-green[i + 1][j] / blue[i + 1][j]) / (2)
					+ (red[i][j + 1] / blue[i][j + 1] - red[i][j - 1] / blue[i][j - 1]) / (2);
				continue;
			}
			else if (j == 0) {
				ln[i][j] = (green[i - 1][j] / blue[i - 1][j] - green[i + 1][j] / blue[i + 1][j]) / (2)
					+ (red[i][j + 1] / blue[i][j + 1]) / (2);
				continue;
			}
			else if (i == height - 1) {
				ln[i][j] = (green[i - 1][j] / blue[i - 1][j]) / (2)
					+ (red[i][j + 1] / blue[i][j + 1] - red[i][j - 1] / blue[i][j - 1]) / (2);
				continue;
			}
			else if (j == width - 1) {
				ln[i][j] = (green[i - 1][j] / blue[i - 1][j] - green[i + 1][j] / blue[i + 1][j]) / (2)
					+ (-red[i][j - 1] / blue[i][j - 1]) / (2);
				continue;
			}
			else {
				//green[i + 1][j]/blue[i+1][j] second
				ln[i][j] = (green[i - 1][j] / blue[i - 1][j] - green[i + 1][j] / blue[i + 1][j]) / (2)
					+ (red[i][j + 1] / blue[i][j + 1] - red[i][j - 1] / blue[i][j - 1]) / (2);
			}
		}
	}

	*/
	/*for (int i = 0; i < width; i++) {
		ln[i] = new double[width];
		for (int j = 0; j < height; j++) {
			if (i == 0 && j == 0) {
				ln[i][j] = (-ytemp[i][j + 1]) / 2.0 + (xtemp[i + 1][j]) / 2.0;
				continue;
			}
			else if (i == width - 1 && j == 0) {
				ln[i][j] = (-ytemp[i][j + 1]) / 2.0 + (-xtemp[i - 1][j]) / 2.0;
				continue;
			}
			else if (i == 0 && j == height - 1) {
				ln[i][j] = ytemp[i][j - 1] / 2.0 + xtemp[i + 1][j] / 2.0;
				continue;
			}
			else if (i == width - 1 && j == height - 1) {
				ln[i][j] = ytemp[i][j - 1] / 2.0 + xtemp[i - 1][j] / 2.0;
				continue;
			}
			else if (i == 0) {
				ln[i][j] = (ytemp[i][j - 1] - ytemp[i][j + 1]) / 2.0
					+ (xtemp[i + 1][j]) / 2.0;
				continue;
			}
			else if (j == 0) {
				ln[i][j] = (-ytemp[i][j + 1]) / 2.0
					+ (xtemp[i + 1][j] - xtemp[i - 1][j]) / 2.0;
				continue;
			}
			else if (i == width - 1) {
				ln[i][j] = (ytemp[i][j - 1] - ytemp[i][j + 1]) / 2.0
					+ (-xtemp[i - 1][j]) / 2.0;
				continue;
			}
			else if (j == height - 1) {
				ln[i][j] = (ytemp[i][j - 1]) / 2.0
					+ (xtemp[i + 1][j] - xtemp[i - 1][j]) / 2.0;
				continue;
			}
			else {
				ln[i][j] = (ytemp[i][j - 1] - ytemp[i][j + 1]) / 2.0
					+ (xtemp[i + 1][j] - xtemp[i - 1][j]) / 2.0;
			}
		}
	}*/
	/*for (int i = 0; i < height; i++) {
		ln[i] = new double[width];
		for (int j = 0; j < width; j++) {
			ln[i][j] = xtemp[i][j] + ytemp[i][j];
		}
	}*/
	cout << "everything is ok before ln calculation" << endl;
	for (int i = 0; i < height; i++) {
		ln[i] = new double[width];
		for (int j = 0; j < width; j++) {
			if (alpha[i][j] > 0.0&&i>1 && j>1 && i<height - 2 && j<width - 2) {
				ln[i][j] = 0.0;
				for (int m = 0; m < 5; m++) {
					for (int n = 0; n < 5; n++) {
						ln[i][j] += (xtemp[i - 2 + m][j - 2 + n] + ytemp[i - 2 + m][j - 2 + n]) * weight[m][n];
						
					}
					
				}
			}
			else {
				ln[i][j] = 0.0;
			}
		}
	}
	for (int i = 0; i < height; i++) {
		for (int j = 0; j < width; j++) {
			
			ln[i][j] = ln[i][j]*canshu;
		}
	}

	minln = INFINITE;
	double maxln = -INFINITE;
	for (int i = 0; i < height; i++) {
		for (int j = 0; j < width; j++) {
			
			if (alpha[i][j]<1.0) {
				ln[i][j] = 0.0;
			}
			if (ln[i][j] > maxln) {
				maxln = ln[i][j];
			}
			if (ln[i][j] < minln) {
				minln = ln[i][j];
			}
		}
	}
	cout << "the min value of ln is: " << minln << endl;
	cout << "the max value of ln is: " << maxln << endl;

	original = maxln - minln;
	/*for (int i = 0; i < height; i++) {
		for (int j = 0; j < width; j++) {
			if (ln[i][j] < 0) {
				ln[i][j] = (ln[i][j]-minln)*(exp(((ln[i][j]-minln) -(maxln)) / maxln)) + minln;
			}
			
		}
	}*/
	
}
//figure out valid pixel in the whole map.
int calculateNumber(double **red, double **green, double **blue) {
	int count = 0;
	for (int i = 0; i < height; i++) {

		for (int j = 0; j < width; j++) {
			if (alpha[i][j]>0) {
				count++;
			}
		}
	}
	return count;
}

VectorXd heightCalculate(int validpixelNumber) {
	bool **validpixel = new bool*[height];
	int  **validOrder = new int*[height];
	int  *rowFromOrder;
	rowFromOrder = new int[validpixelNumber];
	int  *colFromOrder;
	colFromOrder = new int[validpixelNumber];
	int count = 0;

	for (int i = 0; i < height; i++) {
		validpixel[i] = new bool[width];
		validOrder[i] = new int[width];
		for (int j = 0; j < width; j++) {
			//if (alpha[i][j] > 0) {
			if (alpha[i][j]>0.0) {
				validpixel[i][j] = true;
				validOrder[i][j] = count;
				rowFromOrder[count] = i;
				colFromOrder[count] = j;
				count++;
			}
			else {
				validpixel[i][j] = false;
			}
			//gray[i][j] = gray[i][j];
		}
	}
	cout << "valid pixel is : " << count << endl;
	SpMat LMatrix(validpixelNumber, validpixelNumber);

	std::vector<T> coefficients;
	for (int i = 0; i < validpixelNumber; i++) {
		int r = rowFromOrder[i];
		int c = colFromOrder[i];
		for (int j = 0; j < validpixelNumber; j++) {
			if (i == j) {
				coefficients.push_back(T(i, j, 8 + u*u));
			}
			
			if (validpixel[r+1][c]&&j == validOrder[r+1][c]) {
				coefficients.push_back(T(i, j, -1));
			}
			if (validpixel[r + 1][c+1] && j == validOrder[r + 1][c+1]) {
				coefficients.push_back(T(i, j, -1));
			}
			if (validpixel[r - 1][c-1] && j == validOrder[r -1][c-1]) {
				coefficients.push_back(T(i, j, -1));
			}
			if (validpixel[r - 1][c]&&j == validOrder[r-1][c]) {
				coefficients.push_back(T(i, j, -1));
			}
			if (validpixel[r - 1][c+1] && j == validOrder[r - 1][c+1]) {
				coefficients.push_back(T(i, j, -1));
			}
			if (validpixel[r][c + 1]&&j == validOrder[r][c+1]) {
				coefficients.push_back(T(i, j, -1));
			}
			if (validpixel[r + 1][c-1] && j == validOrder[r + 1][c-1]) {
				coefficients.push_back(T(i, j, -1));
			}
			if (validpixel[r][c - 1]&&j == validOrder[r][c-1]) {
				coefficients.push_back(T(i, j, -1));
			}
		}
	}

	

	LMatrix.setFromTriplets(coefficients.begin(), coefficients.end());

	/*Eigen::BiCGSTAB<SparseMatrix<double> >  BCGST;
	BCGST.compute(sparMax);*/
	coefficients.clear();

	/*SpMat TLMatrix(validpixelNumber, validpixelNumber);
	TLMatrix = LMatrix.transpose();
	
	LMatrix = TLMatrix * LMatrix;*/

	Eigen::VectorXd matrixTheta(validpixelNumber);
	//SpMat MatrixC(validpixelNumber, validpixelNumber);
	//MatrixC.diagonal = u*u;
	/*std::vector<T> coefficients2;
	for (int i = 0; i < validpixelNumber; i++) {
		for (int j = 0; j < validpixelNumber; j++) {
			if (i == j) {
				coefficients2.push_back(T(i, j, u*u));
			}
		}
	}
	MatrixC.setFromTriplets(coefficients2.begin(), coefficients2.end());
	coefficients2.clear();*/
	for (int k = 0; k < validpixelNumber; k++) {
		int row = rowFromOrder[k];
		int column = colFromOrder[k];
		if (validpixel[row][column]) {
			matrixTheta(k) = theta;
		}
		else {
			
		}
	}
	

	Eigen::VectorXd matrixPhi(validpixelNumber);
	for (int i = 0; i < validpixelNumber; i++) {
		//matrixPhi(i) = ln[rowFromOrder[i]][colFromOrder[i]];
		matrixPhi(i) = ln[rowFromOrder[i]][colFromOrder[i]];
		/*if (alpha[rowFromOrder[i]][colFromOrder[i] + 1] != 0.0) {
			matrixPhi(i) += 0.5 * ln[rowFromOrder[i]][colFromOrder[i] + 1];
		}
		if (alpha[rowFromOrder[i]][colFromOrder[i] - 1] != 0.0) {
			matrixPhi(i) += 0.5 * ln[rowFromOrder[i]][colFromOrder[i] - 1];
		}
		if (alpha[rowFromOrder[i] + 1][colFromOrder[i]] != 0.0) {
			matrixPhi(i) += 0.5 * ln[rowFromOrder[i] + 1][colFromOrder[i] + 1];
		}
		if (alpha[rowFromOrder[i] - 1][colFromOrder[i]] != 0.0) {
			matrixPhi(i) += 0.5 * ln[rowFromOrder[i] - 1][colFromOrder[i]];
		}*/
	}

	Eigen::VectorXd MatrixB(validpixelNumber);
	//MatrixB = TLMatrix*matrixPhi;
	MatrixB = matrixPhi;
	for (int i = 0; i < validpixelNumber; i++) {
		MatrixB(i) = MatrixB(i) + matrixTheta(i)*u*u;
	}
	
	/*Eigen::VectorXd conb(validpixelNumber);*/
	Eigen::BiCGSTAB<SparseMatrix<double> >  BCGST;
	//LMatrix = LMatrix + u*u*MatrixC;
	BCGST.compute(LMatrix);
	Eigen::VectorXd validheight(validpixelNumber);

	validheight = BCGST.solve(MatrixB);
	std::cout << "everything is ok after solve." << endl;
	Eigen::VectorXd fheight(height*width);

	int index = 0;
	for (int m = 0; m < height; m++) {
		for (int n = 0; n < width; n++) {
			if (!validpixel[m][n]) {
				fheight(m*width + n) = 0.0;
			}
			else {
				fheight(m*width + n) = validheight(index++);
				//fheight(m*width + n) = gray[m][n]*20.0;
			}
		}
	}
	/*Eigen::VectorXd ffheight(height*width);
	for (int i = 0; i < height; i++) {
		
		for (int j = 0; j < width; j++) {
			ffheight(i*width + j) = 0.0;
			if (alpha[i][j] > 0.0&&i>1 && j>1 && i<height - 2 && j<width - 2) {
			
				for (int m = 0; m < 5; m++) {
					for (int n = 0; n < 5; n++) {
						ffheight(i*width + j) += fheight((i - 2 + m)*width + j - 2 + n) * weight[m][n];
					}
				}
			}
			
		}
	}*/

	/*for (int i = 0; i < height; i++) {

		for (int j = 0; j < width; j++) {
			fheight(i*width + j) = 0.0;
			if (alpha[i][j] > 0.0&&i>1 && j>1 && i<height - 2 && j<width - 2) {

				for (int m = 0; m < 5; m++) {
					for (int n = 0; n < 5; n++) {
						fheight(i*width + j) += ffheight((i - 2 + m)*width + j - 2 + n) * weight[m][n];
					}
				}
			}

		}
	}*/
	
	return fheight;
}

int main() {
	SurFacemesh *meshoff;
	meshoff = readPolygon();
	clock_t t1, t2, t3;
	t1 = clock();
	PIXEL *meshnormal;
	meshnormal = makeNormalMap(meshoff);
	treatOriginalData(meshnormal);
	free(meshoff);
	
	initDataTreat(red, green, blue);
	t2 = clock();
	cout << "time consuming of data collection and treatment is : " << (float)(t2 - t1) << endl;
	int validpixelNumber = calculateNumber(red, green, blue);
	Eigen::VectorXd result = heightCalculate(validpixelNumber);
	t3 = clock();
	cout << "time consuming of matrix calculation is : " << (float)(t3 - t2) << endl;
	ofstream myfile;
	std::ofstream ofs("./lionheadnormal.ppm", std::ios::out | std::ios::binary);
	ofs << "P6\n" << width << " " << height << "\n255\n";
	for (int i = 0; i < height; i++) {
		for (int j = 0; j < width; j++) {
			//red[i][j]+1.0)/2.0 * 255
			ofs << (unsigned char)((red[i][j]+1.0)/2.0*255) <<
				(unsigned char)((green[i][j]+1.0)/2.0*255) <<
				(unsigned char)((blue[i][j]+1.0)/2.0*255);
		}
	}
	ofs.close();
	std::ofstream ofs2("./lionhead2ln.ppm", std::ios::out | std::ios::binary);
	ofs2 << "P6\n" << width << " " << height << "\n255\n";
	for (int i = 0; i < height; i++) {
		for (int j = 0; j < width; j++) {
			ofs2 << (unsigned char)((ln[i][j]-minln)/original * 255) <<
				(unsigned char)((ln[i][j]-minln) /original * 255) <<
				(unsigned char)((ln[i][j]-minln) / original * 255);
		}
	}
	ofs2.close();
	myfile.open("elephant-basreliefconvex.off");
	myfile << "OFF" << " " << height*width << " " << (height - 1)*(width - 1) * 2 << " " << 0 << "\n";
	//int i = 0;

	for (double row = 0; row < height; row++) {
		for (double col = 0; col < width; col++) {
			myfile << row << " " << col << " " << result(row*width + col) << "\n";
			//cout << result(row*width + col) << endl;
		}
	}

	for (int seq = 0; seq < height*width; seq++) {
		int leftp = seq, rightp = seq + 1;
		if (rightp%width == 0) {
			continue;
		}
		if (leftp + width < height*width) {
			myfile << 3 << " " << leftp + width << " " << rightp << " " << leftp << "\n";
		}
		if (rightp - width > 0) {
			myfile << 3 << " " << leftp << " " << rightp << " " <<  rightp - width << "\n";
		}
	}
	myfile.close();
	system("pause");
	return 0;
}
