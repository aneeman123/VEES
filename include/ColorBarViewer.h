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
** ColorBarViewer draws a single color bar. 
**/

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

#include <materials.h>
#include <auxiliary.h>
#include <ColorRange.h>
#include <Camera.h>
#include <Viewer.h>

#ifndef _COLORBAR_VIEWER_
#define _COLORBAR_VIEWER_


class  ColorBarViewer:public Viewer {

 public:
  ColorBarViewer();
  void setRange( ColorRange * range );
  void draw(); 
  double getCenterAndMaxDim( double center[3] );
  int getNumElements();

 private:
  ColorRange * theRange;
};

#define BAR_LEN 200
#define BAR_HEIGHT 10

#endif
