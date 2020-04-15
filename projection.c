
#include <GL/glew.h>
#include <GL/freeglut.h>
#include <GL/freeglut_ext.h>
#include "matrix.h"
#include "initShader.h"
#include <math.h>
#include <time.h>
#include<stdio.h>

#define BUFFER_OFFSET( offset )   ((GLvoid*) (offset))
int num_vertices;
GLuint ctm_location;
GLuint proj_location;
vec4 twincube_location = { 0.0, 0.5, 0.0, 1.0 };
vec4 leftcube_location = { -0.5, -0.5, 0.0, 1.0 };
vec4 rightcube_location = { 0.5, -0.5, 0.0, 1.0 };
mat4 proj = { { 1.0, 0.0, 0.0, 0.0 },{ 0.0, 1.0, 0.0, 0.0 },{ 0.0, 0.0, 1.0, 0.0 },{ 0.0, 0.0, 0.0, 1.0 } };

GLfloat twincube_degree = 0.0, leftcube_degree = 0.0, rightcube_degree = 0.0;
mat4 twincube_ctm = { { 1.0, 0.0, 0.0, 0.0 },{ 0.0, 1.0, 0.0, 0.0 },{ 0.0, 0.0, 1.0, 0.0 },{ 0.0, 0.0, 0.0, 1.0 } };
mat4 leftcube_ctm = { { 1.0, 0.0, 0.0, 0.0 },{ 0.0, 1.0, 0.0, 0.0 },{ 0.0, 0.0, 1.0, 0.0 },{ 0.0, 0.0, 0.0, 1.0 } };
mat4 rightcube_ctm = { { 1.0, 0.0, 0.0, 0.0 },{ 0.0, 1.0, 0.0, 0.0 },{ 0.0, 0.0, 1.0, 0.0 },{ 0.0, 0.0, 0.0, 1.0 } };
int s;
int e;
vec4 mazeCentre;
typedef struct
{
	int left;
	int right;
	int top;
	int bottom;
}cell;

cell maze[8][8];

cell setCell(int l, int r, int t, int b)
{
	cell c;
	c.left = l;
	c.right = r;
	c.top = t;
	c.bottom = b;
	return c;
}


void initializemaze()
{
	srand(time(0));
	for (int i = 0; i < 8; i++)
	{
		maze[0][i].top =1; //top walls are all 1
		maze[7][i].bottom =1; //bottom walls are all 1
		maze[i][0].left =1; //left walls are all 1
		maze[i][7].right =1; //right walls are all 1
	}
	s = rand() % 8; //start column
	e = rand() % 8; //end column
	maze[7][s].bottom = 0; //start opening
	maze[0][e].top =0; //end opening
}

void createmaze(int b_up, int b_down, int b_left, int b_right)
{
	//randomly select a cell [i][j]
	//if single cell, return
	srand(time(0));
	if(( b_up == b_down) || (b_left == b_right))
	{
		return;
	}
	int r = rand() % (b_down-b_up); //random row
	r += b_up;
	//printf("%d\n", r);
	int c = rand() % (b_right-b_left); //random column
	c += b_left;
	//printf("%d\n", c);
	for (int i = b_left; i <= b_right; i++)
	{
		//all walls below the point are set to 1
		maze[r][i].bottom = 1;  
		maze[r + 1][i].top = 1;
	}
	for (int i = b_up; i <= b_down; i++)
	{
		//all walls after the point are set to 1
		maze[i][c].right = 1;
		maze[i][c+1].left = 1;
	}
	//take its bottom corner and add walls on top bottom left right dividing maze into 4
	//top vertical wall =0
	//left horizontal wall =1
	//bottom vertical wall =2
	//right horizontal wall =3
	
	int r_array[4] = { 0,1,2,3 };
	int size = 4;
	for (int i = 0; i < 3; i++)
	{
		int random = rand() % size;
		int random_wall = r_array[random];
		for (int k = random; k < (size-1); k++)
		{
			r_array[k] = r_array[k + 1];
		}
		size--;
		if (random_wall == 0)
		{
			//generate random number between b_up and r
			int random_cell = rand() % (r - b_up + 1);
			random_cell += b_up;
			maze[random_cell][c].right = 0;
			maze[random_cell][c + 1].left = 0;

		}
		else if (random_wall == 1)
		{
			//generate random number between b_left and c
			int random_cell = rand() % (c - b_left + 1);
			random_cell += b_left;
			maze[r][random_cell].bottom = 0;
			maze[r + 1][random_cell].top = 0;

		}
		else if (random_wall == 2)
		{
			//generate random number between r+1 and b_down
			int random_cell = rand() % (b_down - r);
			random_cell += r+1;
			maze[random_cell][c].right = 0;
			maze[random_cell][c + 1].left = 0;

		}
		else if (random_wall == 3)
		{
			//generate random number between c+1 and b_right
			int random_cell = rand() % (b_right - c);
			random_cell += c + 1;
			maze[r][random_cell].bottom = 0;
			maze[r + 1][random_cell].top = 0;

		}
	}
	//randomly choose 3 walls and remove walls from those 3
	//recursively do this for the 4 divisions untill each division is one cell
	createmaze(b_up, r, b_left, c); //top-left chamber
	createmaze(b_up, r, c+1, b_right);//top-right chamber
	createmaze(r+1, b_down, b_left, c); //bottom-left chamber
	createmaze(r+1, b_down, c+1, b_right);//bottom-right chamber


}
void printCell(cell c)
{
	if(c.top == 1) printf(" _ \n");
	else printf("   \n");

	if (c.left == 1) printf("|");
	else printf(" ");
	
	if (c.bottom == 1) printf("_");
	else printf(" ");

	if (c.right == 1) printf("|");
	else printf("");

}

void printMaze()
{
	for (int i = 0; i < 8; i++)
	{
		for (int j = 0; j < 8; j++)
		{
			//printCell(maze[i][j]);
			if(maze[i][j].top ==1) printf("*--*");
			else printf("    ");
		}
		printf("\n");
		for (int j = 0; j < 8; j++)
		{
			//printCell(maze[i][j]);
			if (maze[i][j].left == 1)  printf("|  ");
			else printf("   ");

			if (maze[i][j].right == 1) printf("|");
			else printf(" ");
			
			
 
		}
		printf("\n");
		for (int j = 0; j < 8; j++)
		{
			if (maze[i][j].bottom == 1) printf("*--*");
			else printf("    ");
		}
		
		printf("\n");
	}
}

vec4 v4(float x, float y, float z, float w)
{
	vec4 result = { x, y, z, w };
	return result;
}

mat4 matTrans(mat4 m)                        //Performs transformation of a matrix, rows become columns, columns become rows
{
	mat4 result = {
		{ m.x.x, m.y.x, m.z.x, m.w.x },
		{ m.x.y, m.y.y, m.z.y, m.w.y },
		{ m.x.z, m.y.z, m.z.z, m.w.z },
		{ m.x.w, m.y.w, m.z.w, m.w.w }
	};

	return result;
}

vec4 unit(vec4 v)
{
	GLfloat mag = magnitude(v);
	//printf("%f", mag);
	if (mag <= 0.0001)
	{
		mag = 0.0001;
	}
	vec4 result = scalarMult(1 / mag, v);
	return result;

}
float magnitude(vec4 u)
{
	GLfloat u_mag = sqrt(u.x*u.x + u.y*u.y + u.z*u.z + u.w*u.w);
	return u_mag;
}

vec4 scalarMult(float s, vec4 v)             //Multiplies scalar and 4x1 vector 
{
	vec4 result;
	result.x = s*v.x;
	result.y = s*v.y;
	result.z = s*v.z;
	result.w = s*v.w;
	return result;

}
vec4 vecAdd(vec4 u, vec4 v)                 //Adds 2 4x1 vectors
{
	vec4 result;
	result.x = u.x + v.x;
	result.y = u.y + v.y;
	result.z = u.z + v.z;
	result.w = u.w + v.w;
	return result;
}
vec4 vecSub(vec4 u, vec4 v)                // Subtracts 2 4x1 vectors
{
	vec4 result;
	result.x = u.x - v.x;
	result.y = u.y - v.y;
	result.z = u.z - v.z;
	result.w = u.w - v.w;
	return result;
}

vec4 matVecMult(mat4 m, vec4 v)            //multiplies 4x4 matrix with a 4x1 vector resulting in a 4x1 vector.
{
	vec4 result;
	result = vecAdd(vecAdd(scalarMult(v.x, m.x), scalarMult(v.y, m.y)), vecAdd(scalarMult(v.z, m.z), scalarMult(v.w, m.w)));
	return result;
}

mat4 matMult(mat4 m1, mat4 m2)              //Multiplies 2 4x4 matrices resulting in a 4x4 matrix
{
	mat4 result = { matVecMult(m1, m2.x), matVecMult(m1, m2.y), matVecMult(m1, m2.z), matVecMult(m1, m2.w) };
	return result;

}

GLfloat vecDot(vec4 u, vec4 v)               //Performs 4x1 vector dot multiplication , returns a  floating point number
{
	GLfloat d = (u.x*v.x) + (u.y*v.y) + (u.z*v.z) + (u.w*v.w);
	return d;
}

vec4 vecCross(vec4 u, vec4 v)             //Performs 4x1 vector cross multiplication, returns a 4x1 vector
{
	vec4 result;
	result.x = u.y*v.z - u.z*v.y;
	result.y = u.z*v.x - u.x*v.z;
	result.z = u.x*v.y - u.y*v.x;
	result.w = 0;
	return result;
}

mat4 translate(GLfloat x, GLfloat y, GLfloat z) //translation matrix
{
	mat4 result = { { 1.0, 0.0, 0.0, 0.0 },{ 0.0, 1.0, 0.0, 0.0 },{ 0.0, 0.0, 1.0, 0.0 },{ 0.0, 0.0, 0.0, 1.0 } };
	result.w.x = x;
	result.w.y = y;
	result.w.z = z;
	return result;
}


mat4 Scale(GLfloat x, GLfloat y, GLfloat z) //scale matrix
{
	mat4 result = { { 1.0, 0.0, 0.0, 0.0 },{ 0.0, 1.0, 0.0, 0.0 },{ 0.0, 0.0, 1.0, 0.0 },{ 0.0, 0.0, 0.0, 1.0 } };
	result.x.x = x;
	result.y.y = y;
	result.z.z = z;
	return result;
}

mat4 rotate_z(GLfloat theta)  //rotation matrix about z-axis
{
	mat4 result = { { 1.0, 0.0, 0.0, 0.0 },{ 0.0, 1.0, 0.0, 0.0 },{ 0.0, 0.0, 1.0, 0.0 },{ 0.0, 0.0, 0.0, 1.0 } };
	result.x.x = cos(theta);
	result.y.x = -sin(theta);
	result.x.y = sin(theta);
	result.y.y = cos(theta);
	return result;
}

mat4 rotate_x(GLfloat theta) //rotation matrix about x-axis
{
	mat4 result = { { 1.0, 0.0, 0.0, 0.0 },{ 0.0, 1.0, 0.0, 0.0 },{ 0.0, 0.0, 1.0, 0.0 },{ 0.0, 0.0, 0.0, 1.0 } };

	result.y.y = cos(theta);
	result.z.y = -sin(theta);
	result.y.z = sin(theta);
	result.z.z = cos(theta);

	return result;
}

mat4 rotate_y(GLfloat theta) //rotation matrix about y axis
{
	mat4 result = { { 1.0, 0.0, 0.0, 0.0 },{ 0.0, 1.0, 0.0, 0.0 },{ 0.0, 0.0, 1.0, 0.0 },{ 0.0, 0.0, 0.0, 1.0 } };
	result.x.x = cos(theta);
	result.z.x = sin(theta);
	result.x.z = -sin(theta);
	result.z.z = cos(theta);

	return result;
}

mat4 rotate_arb(GLfloat theta, vec4 v) //rotation about arbitrary axis
{
	mat4 result = { { 1.0, 0.0, 0.0, 0.0 },{ 0.0, 1.0, 0.0, 0.0 },{ 0.0, 0.0, 1.0, 0.0 },{ 0.0, 0.0, 0.0, 1.0 } };
	GLfloat d = sqrt((v.y*v.y) + (v.z*v.z));
	if (d == 0)                      // if y and z components are zero, rotate about x-axis
	{

		return rotate_x(theta);

	}
	if (d <= 0.0001)          //to avoid division by zero 
	{
		d = 0.0001;
	}
	mat4 rot_x = { { 1.0, 0.0, 0.0, 0.0 },{ 0.0, 1.0, 0.0, 0.0 },{ 0.0, 0.0, 1.0, 0.0 },{ 0.0, 0.0, 0.0, 1.0 } }; //rotation matrix to bring vector on x-z plane

	rot_x.y.y = v.z / d;
	rot_x.z.y = -v.y / d;
	rot_x.y.z = v.y / d;
	rot_x.z.z = v.z / d;


	mat4 rot_y = { { 1.0, 0.0, 0.0, 0.0 },{ 0.0, 1.0, 0.0, 0.0 },{ 0.0, 0.0, 1.0, 0.0 },{ 0.0, 0.0, 0.0, 1.0 } }; //rotation matrix to align vector with z-axis
	rot_y.x.x = d;
	rot_y.z.x = -v.x;
	rot_y.x.z = v.x;
	rot_y.z.z = d;


	result = matMult(matTrans(rot_x), matMult(matTrans(rot_y), matMult(rotate_z(theta), matMult(rot_y, rot_x)))); //final rotation matrix

	return result;
}

mat4 sMult4(float s, mat4 m)              //Multiplies a scalar and 4x4 matrix
{
	mat4 result;

	result.x.x = s*m.x.x;
	result.x.y = s*m.x.y;
	result.x.z = s*m.x.z;
	result.x.w = s*m.x.w;

	result.y.x = s*m.y.x;
	result.y.y = s*m.y.y;
	result.y.z = s*m.y.z;
	result.y.w = s*m.y.w;

	result.z.x = s*m.z.x;
	result.z.y = s*m.z.y;
	result.z.z = s*m.z.z;
	result.z.w = s*m.z.w;

	result.w.x = s*m.w.x;
	result.w.y = s*m.w.y;
	result.w.z = s*m.w.z;
	result.w.w = s*m.w.w;

	return result;

}

float matdet(float m[9])                   // Calculates determinant of a 3x3 matrix, returns a float
{
	return m[0] * m[4] * m[8] + m[1] * m[5] * m[6] + m[2] * m[3] * m[7] - m[6] * m[4] * m[2] - m[7] * m[5] * m[0] - m[8] * m[3] * m[1];
}

mat4 matMinor(mat4 m)                    // returns matrix minor of a 4x4 matrix. 
{

	mat3 a = { m.y.y, m.z.y, m.w.y, m.y.z, m.z.z, m.w.z, m.y.w, m.z.w, m.w.w };
	mat3 b = { m.y.x, m.z.x, m.w.x, m.y.z, m.z.z, m.w.z, m.y.w, m.z.w, m.w.w };
	mat3 c = { m.y.x, m.z.x, m.w.x, m.y.y, m.z.y, m.w.y, m.y.w, m.z.w, m.w.w };
	mat3 d = { m.y.x, m.z.x, m.w.x, m.y.y, m.z.y, m.w.y, m.y.z, m.z.z, m.w.z };

	mat3 e = { m.x.y, m.z.y, m.w.y, m.x.z, m.z.z, m.w.z, m.x.w, m.z.w, m.w.w };
	mat3 f = { m.x.x, m.z.x, m.w.x, m.x.z, m.z.z, m.w.z, m.x.w, m.z.w, m.w.w };
	mat3 g = { m.x.x, m.z.x, m.w.x, m.x.y, m.z.y, m.w.y, m.x.w, m.z.w, m.w.w };
	mat3 h = { m.x.x, m.z.x, m.w.x, m.x.y, m.z.y, m.w.y, m.x.z, m.z.z, m.w.z };

	mat3 i = { m.x.y, m.y.y, m.w.y, m.x.z, m.y.z, m.w.z, m.x.w, m.y.w, m.w.w };
	mat3 j = { m.x.x, m.y.x, m.w.x, m.x.z, m.y.z, m.w.z, m.x.w, m.y.w, m.w.w };
	mat3 k = { m.x.x, m.y.x, m.w.x, m.x.y, m.y.y, m.w.y, m.x.w, m.y.w, m.w.w };
	mat3 l = { m.x.x, m.y.x, m.w.x, m.x.y, m.y.y, m.w.y, m.x.z, m.y.z, m.w.z };

	mat3 n = { m.x.y, m.y.y, m.z.y, m.x.z, m.y.z, m.z.z, m.x.w, m.y.w, m.z.w };
	mat3 o = { m.x.x, m.y.x, m.z.x, m.x.z, m.y.z, m.z.z, m.x.w, m.y.w, m.z.w };
	mat3 p = { m.x.x, m.y.x, m.z.x, m.x.y, m.y.y, m.z.y, m.x.w, m.y.w, m.z.w };
	mat3 q = { m.x.x, m.y.x, m.z.x, m.x.y, m.y.y, m.z.y, m.x.z, m.y.z, m.z.z };

	mat4 r = { { matdet(a), matdet(b), matdet(c), matdet(d) }, //calculates 3x3 determinant for each entry in a 4x4 matrix.
	{ matdet(e), matdet(f), matdet(g), matdet(h) },
	{ matdet(i), matdet(j), matdet(k), matdet(l) },
	{ matdet(n), matdet(o), matdet(p), matdet(q) } };

	return r;

}
mat4 matCof(mat4 m)                     // turns a matrix minor into cofactor matrix
{
	mat4 r = { { m.x.x, -m.x.y, m.x.z, -m.x.w },
	{ -m.y.x, m.y.y, -m.y.z, m.y.w },
	{ m.z.x, -m.z.y, m.z.z, -m.z.w },
	{ -m.w.x, m.w.y, -m.w.z, m.w.w } };

	return r;
}

float matDet4(mat4 m)                  //Calculates determinant of a 4x4 matrix. Returns a float
{
	mat4 minor = matMinor(m);
	float d = m.x.x* minor.x.x - m.y.x* minor.y.x + m.z.x*minor.z.x - m.w.x*minor.w.x;
	return d;
}

mat4 matInv(mat4 m)                   //Returns inverse of a matrix if determinant is not zero.
{
	if (matDet4(m) != 0)
	{
		mat4 r = sMult4((1 / matDet4(m)), matTrans(matCof(matMinor(m))));
		return r;
	}
	else printf("Matrix is singular, Inverse not possible!");
	return;
}

mat4 look_at(GLfloat eyex, GLfloat eyey, GLfloat eyez, GLfloat atx, GLfloat aty, GLfloat atz, GLfloat upx, GLfloat upy, GLfloat upz)
{
	vec4 vpn = v4(eyex - atx, eyey - aty, eyez - atz, 0.0);
	vec4 n = unit(vpn);

	vec4 vup = v4(upx, upy, upz, 0.0);
	vec4 u = vecCross(vup, n);
	u = unit(u);

	vec4 v = vecCross(n, u);
	v = unit(v);

	vec4 e = v4(eyex, eyey, eyez, 1);

	mat4 modelView = { u,v,n,e };
	modelView.x.w = 0;
	modelView.y.w = 0;
	modelView.z.w = 0;

	modelView = matInv(modelView);


	return modelView;

}
mat4 projection(GLfloat Near, GLfloat Far, GLfloat top, GLfloat bottom, GLfloat left, GLfloat right)
{
	mat4 result = { { (-2 * Near) / (right - left), 0,0,0 },{ 0, (-2 * Near) / (top - bottom),0,0 },
	{ (right + left) / (right - left), (top + bottom) / (top - bottom), (Near + Far) / (Far - Near), -1 },{ 0, 0, (-2 * Far*Near) / (Far - Near),0 } };
	return result;

}

//single cube vertices
vec4 cubeVertices[36] = {
	//front
	{ -0.5,-0.5,0.5,1.0 },
	{ 0.5,-0.5, 0.5,1.0 },
	{ 0.5, 0.5, 0.5,1.0 },

	{ 0.5, 0.5, 0.5,1.0 },
	{ -0.5, 0.5, 0.5,1.0 },
	{ -0.5, -0.5, 0.5,1.0 },

	//right 
	{ 0.5, -0.5, 0.5,1 },
	{ 0.5, -0.5, -0.5,1 },
	{ 0.5, 0.5, -0.5,1 },

	{ 0.5, 0.5, -0.5,1 },
	{ 0.5, 0.5, 0.5,1 },
	{ 0.5, -0.5, 0.5,1 },

	//back
	{ -0.5, -0.5, -0.5,1 },
	{ 0.5, -0.5, -0.5,1 },
	{ 0.5, 0.5, -0.5,1 },

	{ 0.5, 0.5, -0.5,1 },
	{ -0.5, 0.5, -0.5,1 },
	{ -0.5, -0.5, -0.5,1 },

	//left
	{ -0.5, -0.5, 0.5,1 },
	{ -0.5, -0.5, -0.5,1 },
	{ -0.5, 0.5, -0.5,1 },

	{ -0.5, 0.5, -0.5,1 },
	{ -0.5, 0.5, 0.5,1 },
	{ -0.5, -0.5, 0.5,1 },

	//top
	{ -0.5, 0.5, 0.5,1 },
	{ 0.5, 0.5, 0.5,1 },
	{ 0.5, 0.5, -0.5,1 },

	{ 0.5, 0.5, -0.5,1 },
	{ -0.5, 0.5, -0.5,1 },
	{ -0.5, 0.5, 0.5,1 },

	//bottom
	{ -0.5, -0.5, 0.5,1 },
	{ 0.5, -0.5, 0.5,1 },
	{ 0.5, -0.5, -0.5,1 },

	{ 0.5, -0.5, -0.5,1 },
	{ -0.5, -0.5, -0.5,1 },
	{ -0.5, -0.5, 0.5,1 },
};
typedef struct
{
	GLfloat x;
	GLfloat y;
} vec2;

vec2 textures[36] = { { 0.5, 0.5 },{ 0.5, 0.0 },{ 0.0, 0.0 },{ 0.0, 0.0 },{ 0.0, 0.5 },{ 0.5, 0.5 },
{ 0.5, 0.5 },{ 0.5, 0.0 },{ 0.0, 0.0 },{ 0.0, 0.0 },{ 0.0, 0.5 },{ 0.5, 0.5 },
{ 0.5, 0.5 },{ 0.5, 0.0 },{ 0.0, 0.0 },{ 0.0, 0.0 },{ 0.0, 0.5 },{ 0.5, 0.5 },
{ 0.5, 0.5 },{ 0.5, 0.0 },{ 0.0, 0.0 },{ 0.0, 0.0 },{ 0.0, 0.5 },{ 0.5, 0.5 },
{ 0.5, 0.5 },{ 0.5, 0.0 },{ 0.0, 0.0 },{ 0.0, 0.0 },{ 0.0, 0.5 },{ 0.5, 0.5 },
{ 0.5, 0.5 },{ 0.5, 0.0 },{ 0.0, 0.0 },{ 0.0, 0.0 },{ 0.0, 0.5 },{ 0.5, 0.5 } };

vec2 stonetextures[36] = { { 0.5, 0.5 },{ 1, 0.5 },{ 1, 0.0 },{1, 0.0 },{ 0.5, 0.0 },{ 0.5, 0.5 },
{ 0.5, 0.5 },{ 1, 0.5 },{ 1, 0.25 },{ 1, 0.25 },{ 0.5, 0.25 },{ 0.5, 0.5 },
{ 0.5, 0.5 },{ 1, 0.5 },{ 1, 0.0 },{ 1, 0.0 },{ 0.5, 0.0 },{ 0.5, 0.5 },
{ 0.5, 0.5 },{ 1, 0.5 },{ 1, 0.25 },{ 1, 0.25 },{ 0.5, 0.25 },{ 0.5, 0.5 },
{ 0.5, 0.5 },{ 0.75, 0.5 },{ 0.75, 0.0 },{ 0.75, 0.0 },{ 0.5, 0.0 },{ 0.5, 0.5 },
{ 0.5, 0.5 },{ 0.75, 0.5 },{ 0.75, 0.0 },{ 0.75, 0.0 },{ 0.5, 0.0 },{ 0.5, 0.5 } };

vec2 grasstextures[36] = { { 0.0, 1.0 },{ 0.5, 1.0 },{ 0.5, 0.5 },{ 0.5, 0.5 },{ 0.0, 0.5 },{ 0.0, 1.0 },
{ 0.0, 1.0 },{ 0.5, 1.0 },{ 0.5, 0.5 },{ 0.5, 0.5 },{ 0.0, 0.5 },{ 0.0, 1.0 },
{ 0.0, 1.0 },{ 0.5, 1.0 },{ 0.5, 0.5 },{ 0.5, 0.5 },{ 0.0, 0.5 },{ 0.0, 1.0 },
{ 0.0, 1.0 },{ 0.5, 1.0 },{ 0.5, 0.5 },{ 0.5, 0.5 },{ 0.0, 0.5 },{ 0.0, 1.0 },
{ 0.0, 1.0 },{ 0.5, 1.0 },{ 0.5, 0.5 },{ 0.5, 0.5 },{ 0.0, 0.5 },{ 0.0, 1.0 },
{ 0.0, 1.0 },{ 0.5, 1.0 },{ 0.5, 0.5 },{ 0.5, 0.5 },{ 0.0, 0.5 },{ 0.0, 1.0 } };


vec4*vertices;
vec2*tex_coords;

void draw_maze()
{
	initializemaze();
	createmaze(0, 7, 0, 7);
	printMaze();
	int cube_count = 0;
	for (int i = 0; i < 8; i++)
	{
		for (int j = 0; j < 8; j++)
		{
			if (maze[i][j].top == 1)
			{
				cube_count++;
			}
			if (maze[i][j].bottom == 1)
			{
				cube_count++;
			}
			if (maze[i][j].right == 1)
			{
				cube_count++;
			}
			if (maze[i][j].left == 1)
			{
				cube_count++;
			}
		}
	}
	cube_count = cube_count * 3;
	cube_count += 64*4;
	printf("%d", cube_count);
	
	num_vertices = cube_count * 36;
	vertices = (vec4 *)malloc(sizeof(vec4) * (num_vertices));
	tex_coords = (vec2 *)malloc(sizeof(vec2) * (num_vertices));
	int index = 0;
	for (int i = 0; i < 8; i++)
	{
		for (int j = 0; j < 8; j++)
		{
			if (maze[i][j].left == 1)
			{
				for (int k = 0; k < 36; k++)
				{
					vertices[index] = matVecMult(matMult(translate(j, -i, 0.0), Scale(0.1, 1, 1)), cubeVertices[k]);
					tex_coords[index] = textures[k];
					index++;
				}
				for (int k = 0; k < 36; k++)
				{
					vertices[index] = matVecMult(matMult(translate(j, -i-0.5, 0.0), Scale(0.2, 0.2, 1.2)), cubeVertices[k]);
					tex_coords[index] = stonetextures[k];
					index++;
				}
				for (int k = 0; k < 36; k++)
				{
					vertices[index] = matVecMult(matMult(translate(j, -i + 0.5, 0.0), Scale(0.2, 0.2, 1.2)), cubeVertices[k]);
					tex_coords[index] = stonetextures[k];
					index++;
				}
			}
			if (maze[i][j].right == 1)
			{
				for (int k = 0; k < 36; k++)
				{
					vertices[index] = matVecMult(matMult(translate(j + 1, -i, 0.0), Scale(0.1, 1, 1)), cubeVertices[k]);
					tex_coords[index] = textures[k];
					index++;
				}
				for (int k = 0; k < 36; k++)
				{
					vertices[index] = matVecMult(matMult(translate(j + 1, -i-0.5, 0.0), Scale(0.2, 0.2, 1.2)), cubeVertices[k]);
					tex_coords[index] = stonetextures[k];
					index++;
				}
				for (int k = 0; k < 36; k++)
				{
					vertices[index] = matVecMult(matMult(translate(j + 1, -i + 0.5, 0.0), Scale(0.2, 0.2, 1.2)), cubeVertices[k]);
					tex_coords[index] = stonetextures[k];
					index++;
				}
			}
			if (maze[i][j].top == 1)
			{
				for (int k = 0; k < 36; k++)
				{
					vertices[index] = matVecMult(matMult(translate(0.5+j, 0.5-i, 0), matMult(rotate_z(3.14 / 2), Scale(0.1, 1, 1))), cubeVertices[k]);
					tex_coords[index] = textures[k];
					index++;
				}
				for (int k = 0; k < 36; k++)
				{
					vertices[index] = matVecMult(matMult(translate(j, 0.5 - i, 0), Scale(0.2, 0.2, 1.2)), cubeVertices[k]);
					tex_coords[index] = stonetextures[k];
					index++;
				}
				for (int k = 0; k < 36; k++)
				{
					vertices[index] = matVecMult(matMult(translate(j, -0.5 - i, 0), Scale(0.2, 0.2, 1.2)), cubeVertices[k]);
					tex_coords[index] = stonetextures[k];
					index++;
				}
			}
			if (maze[i][j].bottom == 1)
			{
				for (int k = 0; k < 36; k++)
				{
					vertices[index] = matVecMult(matMult(translate(0.5+j, -0.5-i, 0), matMult(rotate_z(3.14 / 2), Scale(0.1, 1, 1))), cubeVertices[k]);
					tex_coords[index] = textures[k];
					index++;
				}
				for (int k = 0; k < 36; k++)
				{
					vertices[index] = matVecMult(matMult(translate(j, -0.5 - i, 0), Scale(0.2, 0.2, 1.2)), cubeVertices[k]);
					tex_coords[index] = stonetextures[k];
					index++;
				}
				for (int k = 0; k < 36; k++)
				{
					vertices[index] = matVecMult(matMult(translate(j, 0.5 - i, 0), Scale(0.2, 0.2, 1.2)), cubeVertices[k]);
					tex_coords[index] = stonetextures[k];
					index++;
				}
			}
			
			for (int m = 0; m < 36; m++)
			{
				vertices[index] = matVecMult(matMult(translate(j+2, -i+2, -0.5), Scale(1.0, 1.0, 0.1)), cubeVertices[m]);
				tex_coords[index] = grasstextures[m];
				index++;

			}
			for (int m = 0; m < 36; m++)
			{
				vertices[index] = matVecMult(matMult(translate(j-1 , -i - 2, -0.5), Scale(1.0, 1.0, 0.1)), cubeVertices[m]);
				tex_coords[index] = grasstextures[m];
				index++;

			}
			for (int m = 0; m < 36; m++)
			{
				vertices[index] = matVecMult(matMult(translate(j-1 , -i + 2, -0.5), Scale(1.0, 1.0, 0.1)), cubeVertices[m]);
				tex_coords[index] = grasstextures[m];
				index++;

			}
			for (int m = 0; m < 36; m++)
			{
				vertices[index] = matVecMult(matMult(translate(j + 2, -i - 2, -0.5), Scale(1.0, 1.0, 0.1)), cubeVertices[m]);
				tex_coords[index] = grasstextures[m];
				index++;

			}
		}
		
	}
	mazeCentre.x = (8.1 - 0.1) / 2;
	mazeCentre.y = (-7.5 - 0.5) / 2;

	for (int i = 0; i < (num_vertices); i++)
	{
		vertices[i] = matVecMult(matMult(rotate_x(-3.14/2),translate(-mazeCentre.x, -mazeCentre.y, 0)), vertices[i]);
	}
	
	
	
}

void init(void)
{
	int width = 512;
	int height = 512;
	GLubyte my_texels[512][512][3];

	FILE *fp = fopen("p2texture04.raw", "r");
	fread(my_texels, width * height * 3, 1, fp);
	fclose(fp);
	
	
	GLuint program = initShader("vshader.glsl", "fshader.glsl");
	glUseProgram(program);

	GLuint mytex[1];
	glGenTextures(1, mytex);
	glActiveTexture(GL_TEXTURE0 + 0);
	glBindTexture(GL_TEXTURE_2D, mytex[0]);
	glTexImage2D(GL_TEXTURE_2D, 0, GL_RGB, width, height, 0, GL_RGB, GL_UNSIGNED_BYTE, my_texels);
	glTexParameterf(GL_TEXTURE_2D, GL_TEXTURE_WRAP_S, GL_REPEAT);
	glTexParameterf(GL_TEXTURE_2D, GL_TEXTURE_WRAP_T, GL_REPEAT);
	glTexParameterf(GL_TEXTURE_2D, GL_TEXTURE_MAG_FILTER, GL_NEAREST);
	glTexParameterf(GL_TEXTURE_2D, GL_TEXTURE_MIN_FILTER, GL_NEAREST);

	int param;
	glGetTexLevelParameteriv(GL_TEXTURE_2D, 0, GL_TEXTURE_WIDTH, &param);

	draw_maze();
	
	GLuint vao;
	glGenVertexArrays(1, &vao);
	glBindVertexArray(vao);

	GLuint buffer;
	glGenBuffers(1, &buffer);
	glBindBuffer(GL_ARRAY_BUFFER, buffer);
	glBufferData(GL_ARRAY_BUFFER, sizeof(vec4) * num_vertices + sizeof(vec2) * num_vertices, NULL, GL_STATIC_DRAW);
	glBufferSubData(GL_ARRAY_BUFFER, 0, sizeof(vec4) * num_vertices, vertices);
	glBufferSubData(GL_ARRAY_BUFFER, sizeof(vec4) * num_vertices, sizeof(vec2) * num_vertices, tex_coords);
	

	GLuint vPosition = glGetAttribLocation(program, "vPosition");
	glEnableVertexAttribArray(vPosition);
	glVertexAttribPointer(vPosition, 4, GL_FLOAT, GL_FALSE, 0, (GLvoid *)0);

	GLuint vTexCoord = glGetAttribLocation(program, "vTexCoord");
	glEnableVertexAttribArray(vTexCoord);
	glVertexAttribPointer(vTexCoord, 2, GL_FLOAT, GL_FALSE, 0,  (GLvoid *)(sizeof(vec4) * num_vertices));

	GLuint texture_location = glGetUniformLocation(program, "texture");
	glUniform1i(texture_location, 0);

	printf("texture_location: %i\n", texture_location);

	ctm_location = glGetUniformLocation(program, "ctm");
	proj_location = glGetUniformLocation(program, "projection");
	
	glEnable(GL_DEPTH_TEST);
	glClearColor(0.3, 0.8, 0.8,0.5);
	glDepthRange(1, 0);
}

void display(void)
{
	glClear(GL_COLOR_BUFFER_BIT | GL_DEPTH_BUFFER_BIT);

	glPolygonMode(GL_FRONT, GL_FILL);
	glPolygonMode(GL_BACK, GL_FILL);
	glUniformMatrix4fv(ctm_location, 1, GL_FALSE, (GLfloat *)&twincube_ctm);
	glUniformMatrix4fv(proj_location, 1, GL_FALSE, (GLfloat *)&proj);
	glDrawArrays(GL_TRIANGLES, 0, num_vertices);
	

	glutSwapBuffers();
}
float angle = 0;
vec4 eye = { 0,4,6,1 };
vec4 at = { 0,0,0,1 };
int flydown = 0;
int flyover = 1;
int walk = 0;
float alpha = 0;
float nearp=-1;
float farp=-12;
vec4 eyep;
vec4 atp;
int go_right = 0;
int go_left = 0;
int go_up= 1;
int go_down = 0;
int start_row;
int start_column;
int right;
int front;
float pro_left = -1;
float pro_right = 1;
float pro_top = 1;
float pro_bottom = -1;
int right_turn = 0;
int left_turn = 0;
int top_turn = 0;
int down_turn = 0;
void idle()
{
	
	if (flyover == 1)
	{
		eyep = matVecMult(rotate_y(angle*(3.14 / 180)), eye);
		atp = v4(0, 0, 0, 1);
		angle += 0.1;
		if (angle >= 360)
		{
			//flydown
			flyover = 0;
			flydown = 1;
		}
	}
	else if (flydown == 1)
	{
		vec4 entrance = { s + 0.5 , -7.6 -1 ,0,1 };
		entrance = matVecMult(matMult(rotate_x(-3.14 / 2), translate(-mazeCentre.x, -mazeCentre.y, 0)), entrance);
		entrance.y += 0.5;
		
		vec4 atnew = { s + 0.5 , -7.5 ,0,1 };
		//calculate vector(entrance -eye)

		atnew = matVecMult(matMult(rotate_x(-3.14 / 2), translate(-mazeCentre.x, -mazeCentre.y, 0)), atnew);
		atnew.y += 0.5;

		//slowly move towards entrance
		vec4 eyedirection = vecSub(entrance, eye);
		eyep = vecAdd(scalarMult(alpha, eyedirection), eye);

		//calculate vector(entrance -at)
		//slowly move towards entrance
		vec4 atdirection = vecSub(atnew, at);
		atp = vecAdd(scalarMult(alpha, atdirection), at);
		alpha += 0.0005;
		if (alpha >= 1)
		{
			flydown = 0;
			walk = 1;
			eye = entrance;
			at = atnew;
			alpha = 0;
			start_column = s;
			start_row = 7;
			nearp = -0.1;
			printf("(%f , %f , %f )", eyep.x, eyep.y, eyep.z);
			printf("(%f , %f , %f )", entrance.x, entrance.y, entrance.z);
		}



	}
	else if (walk == 1)
	{
		pro_left = -0.1;
		pro_right = 0.1;
		pro_top = 0.1;
		pro_bottom = -0.1;
		if (go_up == 1)
		{
			
			if(start_row < 0)
			{
				front = maze[start_row+1][start_column].top;
				right = -1;
			}
			else
			{
				right = maze[start_row][start_column].right;
				front = maze[start_row][start_column].bottom;
			}
			if (top_turn == 1)
			{
				vec4 atnew = { start_column + 1 , -start_row  ,0,1 }; //new at
				atnew = matVecMult(matMult(rotate_x(-3.14 / 2), translate(-mazeCentre.x, -mazeCentre.y, 0)), atnew);
				vec4 atdirection = vecSub(atnew, at);
				atp = vecAdd(scalarMult(alpha, atdirection), at);
	            /*vec4 next = { start_column + 1+0.5,-start_row- 0.7,0,1 }; //top wall position
				next = matVecMult(matMult(rotate_x(-3.14 / 2), translate(-mazeCentre.x, -mazeCentre.y, 0)), next);
				
				vec4 eyedirection = vecSub(next, eye);
				eyep = vecAdd(scalarMult(alpha, eyedirection), eye);*/
				alpha += 0.005;
				if (alpha >= 1)
				{
					go_left = 1;
					go_up = 0;
					alpha = 0;
					eye = eyep;
					at = atp;
					top_turn = 0;
					start_row += 1;
					

				}
			}
			else if (front == 0) //if current cell does not have a bottom wall, go ahead
			{
				vec4 next = { start_column + 0.5 , -start_row-0.4,0,1}; //top wall position
				next = matVecMult(matMult(rotate_x(-3.14 / 2), translate(-mazeCentre.x, -mazeCentre.y, 0)), next);
				vec4 eyedirection = vecSub(next, eye);
				eyep = vecAdd(scalarMult(alpha, eyedirection), eye);
				vec4 atnew = { start_column + 0.5  , -start_row  ,0,1 }; //new at
				atnew = matVecMult(matMult(rotate_x(-3.14 / 2), translate(-mazeCentre.x, -mazeCentre.y, 0)), atnew);
				vec4 atdirection = vecSub(atnew, at);
				atp = vecAdd(scalarMult(alpha, atdirection), at);
				alpha += 0.0005;
				if (alpha >= 1)
				{
					alpha = 0;
					eye = eyep;
					at = atp;
					if (right == -1)
					{//win condition
						printf("win!");
						walk = 0;
					}
					else if (right == 0)
					{//turn right
						go_right = 1;
						start_column += 1;
						go_up = 0;
					}
					else 
					{ //keep going up
						start_row -= 1; 
					}
					
				}
			}
			else
			{//go left
				start_column -= 1;
				top_turn = 1;
			}
		}
		
		else if (go_right == 1)
		{
			if (start_column>7)
			{
				front = maze[start_row][start_column - 1].right;
			}
			else
			{
				right = maze[start_row][start_column].bottom;
				front = maze[start_row][start_column].left;
			}
			if (right_turn == 1)
			{
				vec4 atnew = { start_column, -start_row -0.5  ,0,1 }; //new at
				atnew = matVecMult(matMult(rotate_x(-3.14 / 2), translate(-mazeCentre.x, -mazeCentre.y, 0)), atnew);
				vec4 atdirection = vecSub(atnew, at);
				atp = vecAdd(scalarMult(alpha, atdirection), at);
				/*
				vec4 next = { start_column - 0.7, -start_row -1,0,1 }; //top wall position
				next = matVecMult(matMult(rotate_x(-3.14 / 2), translate(-mazeCentre.x, -mazeCentre.y, 0)), next);
				vec4 eyedirection = vecSub(next, eye);
				eyep = vecAdd(scalarMult(alpha, eyedirection), eye);
				*/
				alpha += 0.005;
				if (alpha >= 1)
				{
					go_right = 0;
					go_up = 1;
					alpha = 0;
					eye = eyep;
					at = atp;
					right_turn = 0;
					start_column -= 1;
				}
			}
			else if (front == 0) //if current cell does not have a top wall, and
			{
				vec4 next = { start_column-0.4, -start_row, 0,1 }; //top wall position
				next = matVecMult(matMult(rotate_x(-3.14 / 2), translate(-mazeCentre.x, -mazeCentre.y, 0)), next);
				//calculate vector(atnew -eye)
				vec4 eyedirection = vecSub(next, eye);
				eyep = vecAdd(scalarMult(alpha, eyedirection), eye);

				vec4 atnew = { start_column  , -start_row ,0,1 };
				atnew = matVecMult(matMult(rotate_x(-3.14 / 2), translate(-mazeCentre.x, -mazeCentre.y, 0)), atnew);
				vec4 atdirection = vecSub(atnew, at);
				atp = vecAdd(scalarMult(alpha, atdirection), at);
				alpha += 0.0005;
				if (alpha >= 1)
				{
					alpha = 0;
					eye = eyep;
					at = atp;
					
					if (right == 0)
					{//turn down
						go_down = 1;
						start_row += 1;
						go_right = 0;
					}
					else
					{ //keep going right
						start_column += 1; 
					}
					
				}
			}
			
			else
			{ //turn up
				
				start_row -= 1;
				right_turn = 1;
				
				
			}


		}
		else if (go_left == 1)
		{
			if (start_column<0)
			{
				front = maze[start_row][start_column + 1].left;
			}
			
			else
			{
				right = maze[start_row][start_column].top;
				front = maze[start_row][start_column].right;
			}

			if (left_turn == 1)
			{
				
				vec4 atnew = { start_column + 0.7, -start_row + 0.5  ,0,1 }; //new at
				atnew = matVecMult(matMult(rotate_x(-3.14 / 2), translate(-mazeCentre.x, -mazeCentre.y, 0)), atnew);
				vec4 atdirection = vecSub(atnew, at);
				atp = vecAdd(scalarMult(alpha, atdirection), at);
				/*
				vec4 next = { start_column + 0.7, -start_row+1,0,1};//top wall position
				next = matVecMult(matMult(rotate_x(-3.14 / 2), translate(-mazeCentre.x, -mazeCentre.y, 0)), next);
				vec4 eyedirection = vecSub(next, eye);
				eyep = vecAdd(scalarMult(alpha, eyedirection), eye);
				*/
				alpha += 0.005;
				if (alpha >= 1)
				{
					go_down = 1;
					go_left = 0;
					alpha = 0;
					eye = eyep;
					at = atp;
					left_turn = 0;
					start_column += 1;

				}
			}
			else if (front == 0) //if current cell does not have a top wall, and
			{
				vec4 next = { start_column + 0.4,-start_row,0,1 };//top wall position
				next = matVecMult(matMult(rotate_x(-3.14 / 2), translate(-mazeCentre.x, -mazeCentre.y, 0)), next);
				//calculate vector(atnew -eye)
				vec4 eyedirection = vecSub(next, eye);
				eyep = vecAdd(scalarMult(alpha, eyedirection), eye);
				vec4 atnew = { start_column , -start_row ,0,1 };
				atnew = matVecMult(matMult(rotate_x(-3.14 / 2), translate(-mazeCentre.x, -mazeCentre.y, 0)), atnew);
				vec4 atdirection = vecSub(atnew, at);
				atp = vecAdd(scalarMult(alpha, atdirection), at);
				alpha += 0.0005;
				if (alpha >= 1)
				{
					alpha = 0;
					eye = eyep;
					at = atp;
					
					if (right == 0)
					{ //turn up
						go_up = 1;
						start_row -= 1;
						go_left = 0;
					}
					else
					{ //keep going left
						start_column -= 1;
					}
				}
			}
			else
			{
				
				start_row += 1;
				left_turn = 1;
				
			}
		}
		else if (go_down == 1)
		{
			
			if (start_row > 7)
			{
				front = maze[start_row - 1][start_column].bottom;
				right = -1;
			}
			else
			{
				right = maze[start_row][start_column].left;
				front = maze[start_row][start_column].top;
			}
			if (down_turn == 1)
			{
				vec4 atnew = { start_column  - 0.5, -start_row  ,0,1 }; //new at
				atnew = matVecMult(matMult(rotate_x(-3.14 / 2), translate(-mazeCentre.x, -mazeCentre.y, 0)), atnew);
				vec4 atdirection = vecSub(atnew, at);
				atp = vecAdd(scalarMult(alpha, atdirection), at);
				/*
				vec4 next = { start_column+0.5-1,-start_row + 0.7 ,0,1 }; //top wall position
				next = matVecMult(matMult(rotate_x(-3.14 / 2), translate(-mazeCentre.x, -mazeCentre.y, 0)), next);
				vec4 eyedirection = vecSub(next, eye);
				eyep = vecAdd(scalarMult(alpha, eyedirection), eye);
				*/
				alpha += 0.005;
				if (alpha >= 1)
				{
					go_right = 1;
					go_down = 0;
					alpha = 0;
					eye = eyep;
					at = atp;
					down_turn = 0;
					start_row -= 1;

				}
			}
			else if (front == 0) //if current cell does not have a top wall, and
			{
				vec4 next = { start_column + 0.5 , -start_row+0.5,0,1}; //top wall position
				next = matVecMult(matMult(rotate_x(-3.14 / 2), translate(-mazeCentre.x, -mazeCentre.y, 0)), next);
				//calculate vector(atnew -eye)
				vec4 eyedirection = vecSub(next, eye);
				eyep = vecAdd(scalarMult(alpha, eyedirection), eye);
				vec4 atnew = { start_column + 0.5 , -start_row ,0,1 };
				atnew = matVecMult(matMult(rotate_x(-3.14 / 2), translate(-mazeCentre.x, -mazeCentre.y, 0)), atnew);
				vec4 atdirection = vecSub(atnew, at);
				atp = vecAdd(scalarMult(alpha, atdirection), at);
				alpha += 0.0005;
				if (alpha >= 1)
				{
					alpha = 0;
					eye = eyep;
					at = atp;
					if(right == -1)
					{
						printf("wrong exit");
					}
					else if (right == 0)
					{
						go_left = 1;
						start_column -= 1;
						go_down = 0;
					}
					else{ start_row += 1; }
				}
			}
			else
			{
				
				start_column += 1;
				down_turn = 1;
				
				
			}
		}
		
		
	}
	
		
		twincube_ctm = look_at(eyep.x, eyep.y, eyep.z, atp.x, atp.y, atp.z, 0, 1, 0);
		proj = projection(nearp, farp, pro_top, pro_bottom, pro_left, pro_right);
		glutPostRedisplay();
	
}



int main(int argc, char **argv)
{
	glutInit(&argc, argv);
	glutInitDisplayMode(GLUT_RGBA | GLUT_DOUBLE | GLUT_DEPTH);
	glEnable(GL_DEPTH_TEST);
	glutInitWindowSize(512, 512);
	glutInitWindowPosition(100, 100);
	glutCreateWindow("cubes");
	glewInit();
	init();
	glutDisplayFunc(display);
	glutIdleFunc(idle);
	glutMainLoop();
	return 0;
}
