#include <GL/glut.h>
#include <iostream>
#include <fstream>
#include <vector> 
#include <cmath>
#include <string>
#include "Bresenham.h"

using namespace std;

////////////////////
//Global Variables//
////////////////////
float *PixelBuffer;
string inputFile;
char quit;
int windowSizeX = 1000, windowSizeY = 1000, style, mode, polyhedraCount, currentID, halfToning;
float viewXmin, viewXmax, viewYmin, viewYmax, deltaX, deltaY, deltaZ, Delta;
ifstream inFile;
int n;




//Vertex Struct
struct Vertex
{
	float x;
	float y;
	float z;
};

/*

struct RGB
{
	float R;
	float G;
	float B;
};
RGB ka, kd, ks;
*/
Vertex x, f;
double Ia, Il, k, ka, kd, ks;

//////////////////////////
//Function Defininitions//
//////////////////////////
void display();
void writeBack();
void getSettings(int, char*[]);
void getSettings2();
void boundBox();
void setScreen();
void setBoundaryBox();
Vertex toNDCtoPixel(float x, float y, float z, int mode);


//Define a Boundary Box struct
struct Boundary
{
	float Xmin, Ymin, Zmin, Xmax, Ymax, Zmax;
};

Boundary boundaryBox{ 0,0,0,0,0,0 };

//Find max function
float max(float a, float b) {
	if (a > b)
		return a;
	else
		return b;
}

//Find min function
float min(float a, float b) {
	if (a < b)
		return a;
	else
		return b;
}
//Bresenham line drawing
void drawBresenham(int x1, int y1, int x2, int y2, float light1, float light2) {
	Bresenham(x1, x2, y1, y2, PixelBuffer, windowSizeX, light1, light2);
}
///////////////////////////////////
//Polygon Object Class Definition//
///////////////////////////////////
class polyhedraObject {
public:
	int vertexCount;
	int triangleCount;

	struct triangle
	{
		int a;
		int b;
		int c;
		Vertex normals;
		vector<float> lightIntensity;
	};
	
	vector<triangle> triangles;
	vector<Vertex> vertices;
	//vector<RGB> lightIntensity;

	//Set count of vertices
	void setMatrix(int x) {
		vertexCount = x;
	}

	//Set count of edges
	void setTriangles(int x) {
		triangleCount = x;
	}

	//Add a new vertex to the vertice vertex
	void addVertex(float x, float y, float z) {
		Vertex newVertex = { x, y, z };
		vertices.push_back(newVertex);
	}

	//Set new edge connections
	void addTriangle(int a, int b, int c) {
		Vertex normal, vertex1, vertex2, vertex3, u, v; 

		//Get Normal
		vertex1 = vertices[a];
		vertex2 = vertices[b];
		vertex3 = vertices[c];
		u.x = vertex2.x - vertex1.x;
		u.y = vertex2.y - vertex1.y;
		u.z = vertex2.z - vertex1.z;
		v.x = vertex3.x - vertex1.x;
		v.y = vertex3.y - vertex1.y;
		v.z = vertex3.z - vertex1.z;
		normal.x = ((u.y*v.z) - (u.z*v.y))/(sqrt((u.x*u.x)+(u.y*u.y)+(u.z*u.z))+sqrt((v.x*v.x) + (v.y*v.y) + (v.z*v.z)));
		normal.y = ((u.z*v.x) - (u.x*v.z))/(sqrt((u.x*u.x) + (u.y*u.y) + (u.z*u.z)) + sqrt((v.x*v.x) + (v.y*v.y) + (v.z*v.z)));
		normal.z = ((u.x*v.y) - (u.y*v.x))/(sqrt((u.x*u.x) + (u.y*u.y) + (u.z*u.z)) + sqrt((v.x*v.x) + (v.y*v.y) + (v.z*v.z)));

		//Create the triangle
		triangle newTriangle = { a, b, c, normal };
		triangles.push_back(newTriangle);
	}

	///////////////
	//Phong Model//
	///////////////
	float phongModel(int point, Vertex normal) {
		float intensity, specular, ambient, diffuse, fpMag;		
		Vertex fp, p, l, v, r;
		p = vertices[point];
		//Calulate fp
		fp.x = f.x - p.x;
		fp.y = f.y - p.y;
		fp.z = f.z - p.z;
		fpMag = sqrt((fp.x*fp.x) + (fp.y*fp.y) + (fp.z*fp.z));
		
		//calculate light vector
		l.x = x.x - p.x;
		l.y = x.y - p.y;
		l.z = x.z - p.z;
		l.x = l.x / (sqrt((l.x*l.x) + (l.y*l.y) + (l.z*l.z)));
		l.y = l.y / (sqrt((l.x*l.x) + (l.y*l.y) + (l.z*l.z)));
		l.z = l.z / (sqrt((l.x*l.x) + (l.y*l.y) + (l.z*l.z)));
		
		//Caluclate Viewing Vector
		v.x = f.x - p.x;
		v.y = f.y - p.y;
		v.z = f.z - p.z;
		v.x = v.x / (sqrt((v.x*v.x) + (v.y*v.y) + (v.z*v.z)));
		v.y = v.y / (sqrt((v.x*v.x) + (v.y*v.y) + (v.z*v.z)));
		v.z = v.z / (sqrt((v.x*v.x) + (v.y*v.y) + (v.z*v.z)));

		//Calculate Reflection Vector
		float ndotl = (normal.x * l.x) + (normal.y * l.y) + (normal.z * l.z);
		r.x = -l.x + (2*ndotl)*normal.x;
		r.y = -l.y + (2*ndotl)*normal.y;
		r.z = -l.z + (2*ndotl)*normal.z;

		float rdotv = (r.x * v.x) + (r.y * v.y) + (r.z * v.z);
		float ndotv = (normal.x * v.x) + (normal.y * v.y) + (normal.z * v.z);

		ambient = Ia * ka;
		diffuse = kd * ndotl;
		specular = ks * pow(rdotv, n);
		if ((ndotl > 0 && ndotv < 0) || (ndotl < 0 && ndotv > 0)) {
			intensity = ambient;
		}
		else if (rdotv < 0) {
			intensity = ambient + ((Il / (fpMag + k)) * (diffuse));
		}
		else {
			intensity = ambient + ((Il / (fpMag + k)) * (diffuse + specular));
		}
		
		return intensity;
	}

	//////////////////////
	//Painters Algorithm//
	//////////////////////
	vector<int> paintersAlgorithm(char mode) {
		vector<int> layers;
		vector<double> max;
		for (int i = 0; i < triangleCount; i++) {
			layers.push_back(i);
			max.push_back(i);
		}
		
		if (mode == 'x') {//X
			int i, j, temp;
			//Calculate X Maxes
			for (int i = 0; i < triangleCount; i++) {
				if (vertices[triangles[i].a].x > vertices[triangles[i].b].x) {
					max[i] = vertices[triangles[i].a].x;
				}
				else {
					max[i] = vertices[triangles[i].b].x;
				}
				if (max[i] < vertices[triangles[i].c].x) {
					max[i] = vertices[triangles[i].c].x;
				}
			}
			//Calculate X layers
			for (i = 0; i < triangleCount - 1; i++) {
				for (j = 0; j < triangleCount - i - 1; j++) {
					if (max[j] < max[j + 1]) {
						temp = layers[j];
						layers[j] = layers[j + 1];
						layers[j + 1] = temp;
						temp = max[j];
						max[j] = max[j + 1];
						max[j + 1] = temp;
					}
				}		
			}
			return layers;
		}else if(mode == 'y') {//Y
			//Calculate Y Maxes
			for (int i = 0; i < triangleCount; i++) {
				if (vertices[triangles[i].a].y > vertices[triangles[i].b].y) {
					max[i] = vertices[triangles[i].a].y;
				}
				else {
					max[i] = vertices[triangles[i].b].y;
				}
				if (max[i] < vertices[triangles[i].c].y) {
					max[i] = vertices[triangles[i].c].y;
				}
			}
			//Calculate Y layers
			int i, j, temp;
			for (i = 0; i < triangleCount - 1; i++) {
				for (j = 0; j < triangleCount - i - 1; j++) {
					if (max[j] < max[j + 1]) {
						temp = layers[j];
						layers[j] = layers[j + 1];
						layers[j + 1] = temp;
						temp = max[j];
						max[j] = max[j + 1];
						max[j + 1] = temp;
					}
				}
			}
			return layers;
		}
		else {//Z
			//Calculate Z maxes
			for (int i = 0; i < triangleCount; i++) {
				if (vertices[triangles[i].a].z > vertices[triangles[i].b].z) {
					max[i] = vertices[triangles[i].a].z;
				}
				else {
					max[i] = vertices[triangles[i].b].z;
				}
				if (max[i] < vertices[triangles[i].c].z) {
					max[i] = vertices[triangles[i].c].z;
				}
			}
			//Calculate Z layers
			int i, j, temp;
			for (i = 0; i < triangleCount - 1; i++) {
				for (j = 0; j < triangleCount - i - 1; j++) {
					if (max[j] < max[j + 1]) {
						temp = layers[j];
						layers[j] = layers[j + 1];
						layers[j + 1] = temp;
						temp = max[j];
						max[j] = max[j + 1];
						max[j + 1] = temp;
					}
				}
			}
			return layers;
		}
	}

	/////////////////////////////////////
	//Write polygon to the Pixel Buffer//
	/////////////////////////////////////
	void drawPolyhedra() {
		Vertex temp1, temp2, temp3;
		int edge1, edge2, edge3;
		float tempLightA, tempLightB, tempLightC;

		//Phong Model
		for (int i = 0; i < triangleCount; i++) {
			triangles[i].lightIntensity.push_back(phongModel(triangles[i].a, triangles[i].normals));
			triangles[i].lightIntensity.push_back(phongModel(triangles[i].b, triangles[i].normals));
			triangles[i].lightIntensity.push_back(phongModel(triangles[i].c, triangles[i].normals));
			cout << triangles[i].lightIntensity[0] << " " << triangles[i].lightIntensity[1] << " " << triangles[i].lightIntensity[2] << endl;
		}

		//Painter's Algorithm
		vector<int> xLayers;
		vector<int> yLayers;
		vector<int> zLayers;
		xLayers = paintersAlgorithm('x');
		yLayers = paintersAlgorithm('y');
		zLayers = paintersAlgorithm('z');

		for (int i = 0; i < triangleCount; i++) {
			//Draw XY
			edge1 = triangles[zLayers[i]].a;
			edge2 = triangles[zLayers[i]].b;
			edge3 = triangles[zLayers[i]].c;
			temp1 = toNDCtoPixel(vertices[edge1].x, vertices[edge1].y, vertices[edge1].z, 1);
			temp2 = toNDCtoPixel(vertices[edge2].x, vertices[edge2].y, vertices[edge2].z, 1);
			temp3 = toNDCtoPixel(vertices[edge3].x, vertices[edge3].y, vertices[edge3].z, 1);
			tempLightA = triangles[zLayers[i]].lightIntensity[0];
			tempLightB = triangles[zLayers[i]].lightIntensity[1];
			tempLightC = triangles[zLayers[i]].lightIntensity[2];
			drawBresenham(temp1.x, temp1.y, temp2.x, temp2.y, tempLightA, tempLightB);
			drawBresenham(temp2.x, temp2.y, temp3.x, temp3.y, tempLightB, tempLightC);
			drawBresenham(temp3.x, temp3.y, temp1.x, temp1.y, tempLightC, tempLightA);

			//Draw XZ
			edge1 = triangles[yLayers[i]].a;
			edge2 = triangles[yLayers[i]].b;
			edge3 = triangles[yLayers[i]].c;
			temp1 = toNDCtoPixel(vertices[edge1].x, vertices[edge1].y, vertices[edge1].z, 2);
			temp2 = toNDCtoPixel(vertices[edge2].x, vertices[edge2].y, vertices[edge2].z, 2);
			temp3 = toNDCtoPixel(vertices[edge3].x, vertices[edge3].y, vertices[edge3].z, 2);
			tempLightA = triangles[yLayers[i]].lightIntensity[0];
			tempLightB = triangles[yLayers[i]].lightIntensity[1];
			tempLightC = triangles[yLayers[i]].lightIntensity[2];
			drawBresenham(temp1.x, temp1.z, temp2.x, temp2.z, tempLightA, tempLightB);
			drawBresenham(temp2.x, temp2.z, temp3.x, temp3.z, tempLightB, tempLightC);
			drawBresenham(temp3.x, temp3.z, temp1.x, temp1.z, tempLightC, tempLightA);

			//Draw YZ
			edge1 = triangles[xLayers[i]].a;
			edge2 = triangles[xLayers[i]].b;
			edge3 = triangles[xLayers[i]].c;
			temp1 = toNDCtoPixel(vertices[edge1].x, vertices[edge1].y, vertices[edge1].z, 3);
			temp2 = toNDCtoPixel(vertices[edge2].x, vertices[edge2].y, vertices[edge2].z, 3);
			temp3 = toNDCtoPixel(vertices[edge3].x, vertices[edge3].y, vertices[edge3].z, 3);
			tempLightA = triangles[xLayers[i]].lightIntensity[0];
			tempLightB = triangles[xLayers[i]].lightIntensity[1];
			tempLightC = triangles[xLayers[i]].lightIntensity[2];
			drawBresenham(temp1.y, temp1.z, temp2.y, temp2.z, tempLightA, tempLightB);
			drawBresenham(temp2.y, temp2.z, temp3.y, temp3.z, tempLightB, tempLightC);
			drawBresenham(temp3.y, temp3.z, temp1.y, temp1.z, tempLightC, tempLightA);

		}
	}

	//Function to find a new boundary for each polyhedra
	Boundary getBoundary() {
		Boundary instance{ 0,0,0,0,0,0 };
		for (int i = 0; i < vertexCount; i++) {
			instance.Xmax = max(instance.Xmax, vertices[i].x);
			instance.Ymax = max(instance.Ymax, vertices[i].y);
			instance.Zmax = max(instance.Zmax, vertices[i].z);
			instance.Xmin = min(instance.Xmin, vertices[i].x);
			instance.Ymin = min(instance.Ymin, vertices[i].y);
			instance.Zmin = min(instance.Zmin, vertices[i].z);
		}
		return instance;
	}

	

};

///////////////////////////
//Reset Pixel Buffer to 0//
///////////////////////////
void setPixelBuffer(float* PixelBuffer) {
	for (int i = 0; i < windowSizeX; i++) {
		for (int j = 0; j < windowSizeY; j++) {
			PixelBuffer[((windowSizeX * j) + i) * 3] = 0;
		}
	}
}

/////////////////////////////////
//Polygon Vector Initialization//
/////////////////////////////////
vector<polyhedraObject> polyhedras;

/////////////////
//Main Function//
/////////////////
int main(int argc, char *argv[])
{
	//allocate new pixel buffer, need initialization!!
	getSettings(argc, argv);
	setBoundaryBox();
	PixelBuffer = new float[windowSizeX * windowSizeY * 3];
	glutInit(&argc, argv);
	glutInitDisplayMode(GLUT_SINGLE);
	//set window size to windowSizeX by windowSizeX
	glutInitWindowSize(windowSizeX, windowSizeY);
	//set window position
	glutInitWindowPosition(0, 0);

	//create and set main window title
	int MainWindow = glutCreateWindow("Hello Graphics!!");
	glClearColor(0, 0, 0, 0); //clears the buffer of OpenGL
	//sets display function
	while (1) {
		display();
		cout << "Would you like to quit? (y/n)\nChooice: ";
		cin >> quit;
		if (quit == 'y') {
			exit(0);
		}
		getSettings2();
	}

	//glutDisplayFunc(display);

	glutMainLoop();//main display loop, will display until terminate
	return 0;
}



///////////////////////
//Make Pixel Function//
///////////////////////

void makePixel(int x, int y, float* PixelBuffer, float intensity)
{
	//Make sure it is within range
	if (x > 0 && x < 1000 && y > 0 && y < 1000) {
		PixelBuffer[((windowSizeX * y) + x) * 3] = intensity;
		PixelBuffer[(windowSizeX * y + x) * 3 + 1] = intensity;
		PixelBuffer[(windowSizeX * y + x) * 3 + 2] = intensity;
	}	
}

/////////////////////////////////////////////////////////////////////////////
//main display loop, this function will be called again and again by OpenGL//
/////////////////////////////////////////////////////////////////////////////
void display() {

	//Misc.
	glClear(GL_COLOR_BUFFER_BIT);
	glLoadIdentity();
	setPixelBuffer(PixelBuffer);

	setScreen();
	for (int i = 0; i < polyhedraCount; i++)
		polyhedras[i].drawPolyhedra();

	//draws pixel on screen, width and height must match pixel buffer dimension
	glDrawPixels(windowSizeX, windowSizeY, GL_RGB, GL_FLOAT, PixelBuffer);

	//window refresh
	glFlush();
}

//////////////////////////////
//Handles File Input//////////
//Handles Initial User Input//
//////////////////////////////
void getSettings(int argc, char* argv[]) {
	//Make sure that the right number of arguments are passed
	if (argc > 2) {
		cerr << "Too Many Arguments!\nStopping Execution";
		exit(1);
	}
	else if (argc == 1) {
		//If no input file specified, ask for one!
		//cout << "Specify Input File: ";
		//getline(cin, inputFile);
		inputFile = "input.txt";
		inFile.open(inputFile);
		if (!inFile) {
			cerr << "Unable to open input file \nStopping Execution";
			exit(1);
		}
	}
	else {
		//If failed to open, print error and exit
		inFile.open(argv[1]);
		if (!inFile) {
			cerr << "Unable to open input file \nStopping Execution";
			exit(1);
		}
	}
	string space;
	int vertexCount, TriangleCount, point1, point2, point3;
	float x, y, z;
	//Read input file
	inFile >> polyhedraCount;
	polyhedras.resize(polyhedraCount);
	for (int i = 0; i < polyhedraCount; i++) {
		getline(inFile, space);
		inFile >> vertexCount;
		polyhedras[i].setMatrix(vertexCount);
		for (int j = 0; j < vertexCount; j++) {
			inFile >> x;
			inFile >> y;
			inFile >> z;
			polyhedras[i].addVertex(x, y, z);
		}
		inFile >> TriangleCount;
		polyhedras[i].setTriangles(TriangleCount);
		for (int j = 0; j < TriangleCount; j++) {
			inFile >> point1;
			inFile >> point2;
			inFile >> point3;
			polyhedras[i].addTriangle(point1, point2, point3);
		}
	}
	//Go to polyhedra Menu
	getSettings2();

}

void getSettings2() {
	//Ask for the Phong Model Parameters
	cout << "Please enter the parameters of the Phong Model below:\nIa: ";
	cin >> Ia;
	while (Ia < 0 || Ia > 1) {
		cout << "Please Enter a value between 0 and 1\nIa: ";
		cin >> Ia;
	}
	cout << "IL: ";
	cin >> Il;
	while (Il < 0 || Il > 1) {
		cout << "Please Enter a value between 0 and 1\nIL: ";
		cin >> Il;
	}
	cout << "k: ";
	cin >> k;
	while (k < 0) {
		cout << "Please Enter a value that is positive\nk: ";
		cin >> k;
	}

	cout << "n: ";
	cin >> n;
	cout << "Ambient: ";
	cin >> ka;
	while (ka < 0 || ka > 1) {
		cout << "Please Enter a value between 0 and 1\nAmbient: ";
		cin >> ka;
	}
	cout << "Diffuse: ";
	cin >> kd;
	while (kd < 0 || kd > 1) {
		cout << "Please Enter a value between 0 and 1\nDiffuse: ";
		cin >> kd;
	}
	cout << "Specular: ";
	cin >> ks;
	while (ks < 0 || ks > 1) {
		cout << "Please Enter a value between 0 and 1\nSpecular: ";
		cin >> ks;
	}
	
	cout << "X\n  x: ";
	cin >> x.x;
	cout << "  y: ";
	cin >> x.y;
	cout << "  z: ";
	cin >> x.z;
	cout << "f\n  x: ";
	cin >> f.x;
	cout << "  y: ";
	cin >> f.y;
	cout << "  z: ";
	cin >> f.z;

	return;
}

/////////////////////////////////////////////////////////////////////////////////
//Sets Screen to divide it by 4 differnt quadrants, and adds which one is which//
/////////////////////////////////////////////////////////////////////////////////
void setScreen() {
	//Draw borders
	drawBresenham(0, 500, 1000, 500, 1, 1);
	drawBresenham(500, 0, 500, 1000, 1, 1);

	//Draw XY
	drawBresenham(0, 520, 10, 500, 1, 1);
	drawBresenham(10, 520, 0, 500, 1, 1);
	drawBresenham(20, 520, 25, 510, 1, 1);
	drawBresenham(30, 520, 20, 500, 1, 1);

	//Draw XZ
	drawBresenham(500, 520, 510, 500, 1, 1);
	drawBresenham(510, 520, 500, 500, 1, 1);
	drawBresenham(520, 520, 530, 520, 1, 1);
	drawBresenham(530, 520, 520, 503, 1, 1);
	drawBresenham(520, 503, 530, 503, 1, 1);

	//Draw YZ
	drawBresenham(0, 20, 5, 10, 1, 1);
	drawBresenham(10, 20, 0, 0, 1, 1);
	drawBresenham(20, 20, 30, 20, 1, 1);
	drawBresenham(30, 20, 20, 3, 1, 1);
	drawBresenham(20, 3, 30, 3, 1, 1);
}

//Set a new boundary box if values are outside the current range!
void setBoundaryBox() {
	Boundary temp;
	for (int i = 0; i < polyhedraCount; i++) {
		//Get the boundaries of the new values
		temp = polyhedras[i].getBoundary();
		//Check if the new values cause a new boundary, if so change boundary Box size
		//Add or subtract 10% to have some space between the screen edges
		if (max(boundaryBox.Xmax, temp.Xmax) > boundaryBox.Xmax) {
			boundaryBox.Xmax = temp.Xmax;
		}
		if (max(boundaryBox.Xmax, temp.Ymax) > boundaryBox.Ymax) {
			boundaryBox.Ymax = temp.Ymax;
		}
		if (max(boundaryBox.Zmax, temp.Zmax) > boundaryBox.Zmax) {
			boundaryBox.Ymax = temp.Ymax;
		}
		if (min(boundaryBox.Xmin, temp.Xmin) < boundaryBox.Xmin) {
			boundaryBox.Xmin = temp.Xmin;
		}
		if (min(boundaryBox.Ymin, temp.Ymin) < boundaryBox.Ymin) {
			boundaryBox.Ymin = temp.Ymin;
		}
		if (min(boundaryBox.Zmin, temp.Zmin) < boundaryBox.Zmin) {
			boundaryBox.Zmin = temp.Zmin;
		}
	}
	//Create new values for deltas
	deltaX = boundaryBox.Xmax - boundaryBox.Xmin;
	deltaY = boundaryBox.Ymax - boundaryBox.Ymin;
	deltaZ = boundaryBox.Zmax - boundaryBox.Zmin;
	//Take max of the deltas
	Delta = max(deltaX, deltaY);
	Delta = max(Delta, deltaZ);
}

//Turn world to NDC then to Pixel
Vertex toNDCtoPixel(float x, float y, float z, int mode) {
	Vertex point{ 0,0,0 };
	float xNDC, yNDC, zNDC;
	//If mode is 1, set for XY quadrant
	if (mode == 1) {
		xNDC = (x - boundaryBox.Xmin) / Delta;
		yNDC = (y - boundaryBox.Ymin) / Delta;
		point.x = int(xNDC * 300);
		point.y = int(yNDC * 300);
		point.y += 510;
		point.x += 10;
		return point;
	}//If mode is 2, set for XZ quadrant
	else if (mode == 2) {
		xNDC = (x - boundaryBox.Xmin) / Delta;
		zNDC = (z - boundaryBox.Zmin) / Delta;
		point.x = int(xNDC * 300);
		point.z = int(zNDC * 300);
		point.x += 510;
		point.z += 510;
		return point;
	}
	//Else it set for YZ quadrant
	yNDC = (y - boundaryBox.Ymin) / Delta;
	zNDC = (z - boundaryBox.Zmin) / Delta;
	point.y = int(yNDC * 300);
	point.z = int(zNDC * 300);
	point.y += 10;
	point.z += 10;
	return point;
}