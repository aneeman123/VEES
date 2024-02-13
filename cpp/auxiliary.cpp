/* Utilities for normal calculation 
 * I think I may have gotten these from nate Robbins' online
 * tutorial
 * Alisa Neeman
 * <aneeman@cse.ucsc.edu>
 */

#include <auxiliary.h>

double magnitude( const double v1[3]) {
  double m;
  m = v1[0]*v1[0] +  v1[1]*v1[1] +  v1[2]*v1[2];
  return sqrt(m);
}

/* normalize a 3D vector from Nate Robbins' tutorials*/
double normalize(double v[3])
{
    double length;

    length = sqrt(v[0]*v[0] + v[1]*v[1] + v[2]*v[2]);
    v[0] /= length;
    v[1] /= length;
    v[2] /= length;

    return length;
}

/* 
 * calculate surface normal for a triangle 
 * Normals are returned in n[3]
 */
void calcNormal(
   const double pt1[3],const double pt2[3],const double pt3[3], double n[3]) {
	double a[3], b[3];
	a[0] = pt2[0] - pt1[0];
	a[1] = pt2[1] - pt1[1];
	a[2] = pt2[2] - pt1[2];

	b[0] = pt3[0] - pt1[0];
	b[1] = pt3[1] - pt1[1];
	b[2] = pt3[2] - pt1[2];

	//cross product of a and b
	n[0] = a[1]*b[2] - a[2]*b[1];
	n[1] = a[2]*b[0] - a[0]*b[2];
	n[2] = a[0]*b[1] - a[1]*b[0];

}



/* To convert degrees into radians you 
   should multiply by pi/180. */
double degreesToRadians(double degrees) {
  return ((degrees*3.14159)/180);
}

double radiansToDegrees(double radians) {

  return ((radians*180)/3.14159);
}

/* dot product of two 3D vectors */
double dot_product( const double v1[3], const double v2[3] ) {
  double zero_prod,one_prod,two_prod, sum;
  zero_prod = v1[0]*v2[0];
  one_prod = v1[1]*v2[1];
  two_prod = v1[2]*v2[2];
  
  sum = two_prod + one_prod + zero_prod;
  return sum;
  
}

/* cross product of two 3D vectors */
void cross_product( const double a[3],const double b[3], double cross[3] ) 
{
  cross[0] = a[1]*b[2] - a[2]*b[1];
  cross[1] = a[2]*b[0] - a[0]*b[2];
  cross[2] = a[0]*b[1] - a[1]*b[0];

}

/* transform a vector of length n by multiplying with an n x n matrix 
 * pre! array bounds are ok. conventional matrix[row][column] 
 *
 */
void transform_vector(double * start_vector, double * matrix, double * end_vector, int n) 
{  
  int i,j;
  double result;
  for( i = 0; i < n; i++ ) { /* for each column in matrix */
    result = 0;
    for (j = 0; j < n; j++) {  /* for each start_vector entry */
	result = result + (start_vector[j] * matrix[(i*n) + j]);
    }
    end_vector[i] = result;
  }
}

/*
 * transpose an n x n matrix sorted in array length (n x n)
 */
void transpose_matrix( double before[], double after[], int n )
{
  int i, j;
  for ( i = 0; i < n; i++ ) {
    for( j = 0; j < n; j++)
    after[i*n +j] = before[j*n + i];
  }
}




/** Set glut font  */
//int font_style;
void setfont(char* name, int size)
{
  font_style = *((unsigned int *)GLUT_BITMAP_HELVETICA_10);
    if (strcmp(name, "helvetica") == 0) {
	if (size == 12) 
	  font_style = *( (unsigned int *)GLUT_BITMAP_HELVETICA_12);
	else if (size == 18)
	  font_style =  *((unsigned int *)GLUT_BITMAP_HELVETICA_18);
    } else if (strcmp(name, "times roman") == 0) {
      font_style =  *((unsigned int * )GLUT_BITMAP_TIMES_ROMAN_10);
	if (size == 24)
	  font_style =  *( (unsigned int *)GLUT_BITMAP_TIMES_ROMAN_24);
    } else if (strcmp(name, "8x13") == 0) {
      font_style =  *((unsigned int *)GLUT_BITMAP_8_BY_13);
    } else if (strcmp(name, "9x15") == 0) {
      font_style =  *((unsigned int *)GLUT_BITMAP_9_BY_15);
    }
}

/** 
 * drawing a string on the lower_left: position the raster
 * and call glutBitmapCharacter, character by character,
 * on the string
 */

void drawstr(GLuint x, GLuint y, char* format, ...)
{
    va_list args;
    char buffer[255], *s;
 
   
    va_start(args, format);
    vsprintf(buffer, format, args);
    va_end(args);
    
    glRasterPos2i(x, y);
    for (s = buffer; *s; s++)
	glutBitmapCharacter((void *)font_style, *s);
   
}


