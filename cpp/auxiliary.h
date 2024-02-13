#ifndef _AUXILIARY_H_
#define _AUXILIARY_H_

#include <math.h>
#include <stdio.h>
#include <stdlib.h>
#include <stdarg.h>
#include <string.h>

#ifdef WIN32
#include <windows.h>
#endif

#if defined(__APPLE__)&& defined(__MACH__)
#include <GLUT/glut.h>
#include <OpenGL/gl.h>
#include <OpenGL/glu.h>
#else
#include <GL/gl.h>			/* OpenGL header file */
#include <GL/glu.h>			/* OpenGL utilities header file */
#include <GL/glut.h>	
#endif

#define samesign(i,j) \
( ((i >= 0 && j >=0) || (i < 0 && j < 0)) ? (1) : (0) )

static int font_style;

/* normalize a 3D vector from Nate Robbins' tutorials*/
double normalize(double v[3]);

/* 
 * calculate surface normal for a triangle 
 * Normals are returned in n[3]
 */
void calcNormal(const double pt1[3],const double pt2[3],const double pt3[3],double n[3]);


double degreesToRadians(double degrees);
double radiansToDegrees(double radians);

/* dot produnt of 2 3D vectors */
double dot_product( const double v1[3], const double v2[3] );
double magnitude( const double v1[3]);

/* cross product of 2 3D vectors */
void cross_product( const double a[3], const double b[3], double cross[3] );

/* transform a vector of length n using an n x n matrix */
void transform_vector(double * start_vector, double * matrix, double * end_vector, int n); 

/* transpose an n x n matrix sored in array length (n x n) */
void transpose_matrix( double before[], double after[], int n );


void setfont(char* name, int size);
void drawstr(GLuint x, GLuint y, char* format, ...);


#endif
