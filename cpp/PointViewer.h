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
** PointViewer draws a single gauss point or fiber.
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

// includes for the domain classes
#include <Domain.h>
#include <Node.h>

#include "Camera.h"
#include "materials.h"
#include "auxiliary.h"
#include "ColorRange.h"

#include "DrawElement.h"


#include "Viewer.h"

#ifndef _POINT_VIEWER_
#define _POINT_VIEWER_

class  PointViewer:public Viewer {

 public:
  PointViewer();
  void setDomain( Domain * d );
  void draw(); 
  double getCenterAndMaxDim( double center[3] );
  int getNumElements();
  void setElement( int num );


 private:
  Domain * theDomain;
  Element * theElement;
  int eleNum; //the element to be drawn
  int gaussNum;
  int fiberLocation;

  double max[3]; // element bounding box
  double min[3];
  double maxDim;
 
};

#endif
