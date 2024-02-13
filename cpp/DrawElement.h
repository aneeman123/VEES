/**
** (C) Copyright 2006, The Regents of the University of California    **
** All Rights Reserved.                                               **
**                                                                    **
** Commercial use of this program without express permission of the   **
** University of California, Berkeley, is strictly prohibited.  See   **
** file 'COPYRIGHT'  in main directory for information on usage and   **
** redistribution,  and for a DISCLAIMER OF ALL WARRANTIES.           **
**
** Developed by:
** Alisa Neeman (aneeman@cse.ucsc.edu) 
 **
 ** DrawElement is an abstract class for a drawing method.
 ** The idea is to have lots of ways to draw an element.
 **/

#ifndef _DRAW_ELEMENT_
#define _DRAW_ELEMENT_
// includes for the domain classes
#include <Domain.h>
#include <Node.h>

/* helper classes for draw */
#include <SingleDomEleIter.h>
#include <SingleDomNodIter.h>
#include <NodeIter.h>

#include <math.h>
#include <stdio.h>
#include <stdio.h>

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
#include "materials.h"
#include "ColorRange.h"
/**
 **
 ** DrawElement.cpp is an abstract class for elements
 ** with different kinds of geometry. It carries the scalars
 ** for color along with draw methods
 **/



class DrawElement {

 public:
  DrawElement();
  virtual ~DrawElement();
  virtual void drawElement(Element * el, double * colors, 
		      	   ColorRange range, int colorStyleFlag );
  virtual bool canDraw( Element * el );
  
  virtual int getClassTag();

};
/* For 20 node brick, 8 nodes make a loop (top & bottom faces) 
 * below is the drawing order (counterclockwise for each face) */

static int brick20Loops[2][8] ={ {5,13,6,14,7,15,8,16},//bottom face
                           {1,9,2,10,3,11,4,12}};//top face
static int brick20Verts[4][3]= { {5,17,1}, {6,18,2}, {7,19,3},{8,20,4}};

/* for brick8n, parallelpiped sides */
static int brick8Loops[4][4]= {
  {4,5,1,0},{3,7,4,0},{7,6,2,3},{1,5,6,2}
  };

static int oppose8Loops[4][4]= { 
{2,3,6,7}, {1,2,5,6}, {0,1,4,5}, {0,3,4,7} }; // nodes in face opposite   
                                            

enum DrawElementType{ DRAW_WIRE, DRAW_WIRE_DISPLACED, 
		      DRAW_SOLID,DRAW_SOLID_DISPLACED,DRAW_GAUSS_PTS,
		      DRAW_STIFFNESS,
			    HIDE_ELEMENT,  DRAW_PLANE_IN_A_BOX
		     };

enum colorStyle{ BY_NODE, BY_ELEMENT,BY_GAUSS_PT, SINGLE_COLOR };

#endif
