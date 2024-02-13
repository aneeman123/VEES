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

** camera.h: General header for three dimensional camera object
******************************************************************/


#ifndef CAMERA_H
#define CAMERA_H


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

#define X 0
#define Y 1
#define Z 2



class  Camera {
 

 public:
  Camera();
  void initCamera( int z_length );
  void lookAt();
  void calculateCameraPosition();
  void setZoom( double zoom );
  void zoomIn();
  void zoomOut();
  double getZoom();

  private:
  /* gluLookAt */
  GLdouble eye[3]; //camera position
  GLdouble at[3];  //position we're looking at
  GLdouble up[3];  //looking down "z", which direction is up

 
  GLdouble NO_ZOOM;
  GLdouble curZoom; // distance (R) to object we're looking at
  GLdouble curThetaX;   //how much azimuth
  GLdouble curThetaY;    //how much elevation


};

#endif	/* CAMERA_H */


