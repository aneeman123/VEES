/*****************************************************************
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
** Filename: VisWin.h OpenGL window for 3D viewing
**  GUI API independent.
******************************************************************/


#ifndef VISWIN_H
#define VISWIN_H


#include <stdio.h>
#include "Camera.h"              //camera abstraction


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

#define X 0
#define Y 1
#define Z 2

#define HIT_BUFSIZE 512
enum MOUSE{MOUSE_PRESS,MOUSE_RELEASE,MOUSE_LEFT_BUTTON,MOUSE_MIDDLE_BUTTON,
		   MOUSE_RIGHT_BUTTON};

#include "Viewer.h"
#include "auxiliary.h"
#include "materials.h"
class VisWin {
  


 public:
  VisWin ( Viewer *domv, int width, int height );
  void init( int width, int height );
  void draw();
  void drawIncrNumber( char * format, int num );
  void drawTwoNumbers( char * format, int num1, int num2 );
  int handle(int button, int state, int x, int y); // handling mouse events
  int pick(int x, int y);
  void rotateVolume();
  void centerVolume();
  void drawStaticScene();
  void reshape( int width, int height );

  /*********************/
  void initLight();
 void lightMinusX(); 
 void lightMinusY();  
 void lightMinusZ();
 void lightAbove();
 
 void initCamera();
  
  void zoomIn();
  void zoomOut();
  void spinX();
  void spinY();
  void spinZ();
  void resetView();
  void setOrtho();
   /* set viewing parameters */
   void setZoom(double zoom);   
   void setDrawMethod( int drawMethod );

 private:

   Camera cam; 
   int oldX;
   int oldY; // follow mouse
 
   int followMouse; // 0 don't foolow, 1 follow (left mouse
   int xSpin, ySpin, zSpin; 
   double center[3]; // center of volume
   int width, height;
   double scale;
   double maxDim;
   Viewer * theViewer;
};

#endif	/* VISWIN_H */


