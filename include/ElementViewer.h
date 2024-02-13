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
** ElementViewer draws a single element. Integration points or individual
** fibers can be selected
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
#include "DrawWireFrame.h"
#include "DrawDisplacedWireFrame.h"
#include "DrawSolid.h"
#include "DrawSolidDisplaced.h"
//#include "DrawPlaneInABox.h"
#include "HideElement.h"
#include "DrawGaussPoints.h"

#include "Viewer.h"

#ifndef _ELEMENT_VIEWER_
#define _ELEMENT_VIEWER_

class  ElementViewer:public Viewer {

 public:
  ElementViewer();
  void setDomain( Domain * d );
  void draw(); 
  double getCenterAndMaxDim( double center[3] );
  int getNumElements();
  void setElement( int num );
  int getElement();
  void setSelected( int i ); // should this be in superclass?
  void  drawCompass();

 private:
  Domain * theDomain;
  Element * theElement;
  int eleNum; //the element to be drawn

  double max[3]; // element bounding box
  double min[3];
  double maxDim;

  DrawElement * drawElementTypes[10]; //these are singleton instances of each 
                                   //type of DrawElement object
                                    //NOT REQUIRING STATE

};

#endif
