#include <GL/glew.h>
#include <GL/freeglut.h>
#include <GL/freeglut_ext.h>
typedef struct              //4x1 Vector 
{
	float x;
	float y;
	float z;
	float w;
}vec4;

vec4 v4(float, float, float, float);
typedef struct               //4x4 matrix
{
	vec4 x;
	vec4 y;
	vec4 z;
	vec4 w;
}mat4;

typedef float mat3[9];      //3x3 matrix 1-D array 

void printvect(vec4);
vec4 scalarMult(GLfloat, vec4);
vec4 vecAdd(vec4, vec4);
vec4 vecSub(vec4, vec4);
GLfloat vecDot(vec4, vec4);
vec4 vecCross(vec4, vec4);
void printmat(mat4);
mat4 sMult4(GLfloat, mat4);
mat4 matAdd(mat4, mat4);
mat4 matSub(mat4, mat4);
mat4 matTrans(mat4);
vec4 matVecMult(mat4, vec4);
mat4 matMult(mat4, mat4);
mat4 matMinor(mat4);
GLfloat matdet(mat3);
mat4 matCof(mat4);
GLfloat matDet4(mat4);
mat4 matInv(mat4);
mat4 translate(GLfloat,GLfloat, GLfloat);
mat4 Scale(GLfloat, GLfloat, GLfloat);
mat4 rotate_z(GLfloat);
mat4 rotate_x(GLfloat);
mat4 rotate_y(GLfloat);
mat4 rotate_arb(GLfloat, vec4 v);
vec4 unit(vec4 );
GLfloat magnitude(vec4 );

