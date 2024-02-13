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
** DomainViewer.h: Header for
** Wrapper class for OpenSees Domain to provide management for visualization
**/

#include <stdlib.h>
#include <OPS_Globals.h>
#include <StandardStream.h>
#include <ArrayOfTaggedObjects.h>

// includes for the domain classes
#include <Domain.h>
#include <Node.h>


/* helper classes for draw */
#include <SingleDomEleIter.h>
#include <SingleDomNodIter.h>
#include <NodeIter.h>

#include <math.h>
#include <stdio.h>

#ifdef WIN32
#include <windows.h>
#include <map>
#else
//#include <map.h>
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

#include <File_Channel.h> // alisa's class

#include "Camera.h"
#include "materials.h"
#include "auxiliary.h"
#include "ColorRange.h"

#include "DrawElement.h"
#include "DrawWireFrame.h"
#include "DrawDisplacedWireFrame.h"
#include "DrawSolid.h"
#include "DrawSolidDisplaced.h"
#include "HideElement.h"
#include "DrawGaussPoints.h"
#include "DrawStiffness.h"

#include "ElementStrings.h"



#ifndef _DOMAIN_VIEWER_H
#define _DOMAIN_VIEWER_H

/* 
 * DomainViewer.h
 * @author Alisa Neeman
 * Wrapper class for OpenSees Domain to provide management for visualization
 */



/* in DrawElement:

enum DrawElementType{ DRAW_WIRE, DRAW_WIRE_DISPLACED, 
		      DRAW_SOLID,DRAW_SOLID_DISPLACED,DRAW_GAUSS_PTS,
			    HIDE_ELEMENT,  DRAW_PLANE_IN_A_BOX};
*/


#include "ColorSchemeNames.h"

#include "Viewer.h"

class  DomainViewer:public Viewer {
 

 public:
  DomainViewer();
  void setDomain( Domain * d );
  void initRanges();
  void draw(); // draw elements
  void drawNodes();

  // methods for reading in state from file
  int updateState(const char * filename, int commitStep);
  void updateRanges();
  void updateStiffness();

  void setDrawMethodForAll( int drawMethod );
  void setDrawMethodForOne( int eleNum, int drawMethod );
  bool setScalarColorMethodForAll( int drawMethod );
  void setColorRange( double min, double max, int invert );



  void initGauss(); // allocate and set locations
  bool setMeanStress();
  void setDeviatoricStress();
  bool setDisplacement( double * theData );


  void populateElementShowList( int *names, int & numTypes );
  void toggleShowElement( int classTag, int hideElement );


  void   setBasicBounds();
  double getCenterAndMaxDim( double center[3] );

  const ColorRange& getColorRange(int colorScheme);
  void updateColorMap( bool isLogScale, double dataRange); //update currRange
  
  void   drawCompass(); //legend for axes

  int getNumElements();
  void saveMeanStress();
  void saveDeviatoricStress();
  void saveStiffness( int commitStep );
  void saveVTKdata( int commitStep );
  void saveStressTensor( int commitStep );
  void saveRotationAccuracy( int commitStep );
  void saveRegularVolumeStiffness( int commitStep );
  void setMinStiffnessEigenvalue( );
  void saveEigenTensor( int commitStep);
  void saveStressLodeAngle( int commitStep);
  void saveEigenLodeAngle( int commitStep);
  double cubeRoot(double val);
  int findElement(double loc[3] );
  double minDistBetweenPts();

  int Matrix2TensorSysR4(const Matrix& M, double T[4][4][4][4]);

  private:
  Domain * theDomain;
 
  int maxNodeIndex;
  int maxEleIndex;
  int maxEleTypeIndex;
 

  int drawStyle; //BY_NODE, BY_ELEMENT, BY_GAUSS_PT
  int drawMethod; //DRAW_WIREFRAME, etc 


  ColorRange ranges[NUM_COLOR_SCHEMES];
  
  double * data0; // color maps
  double * data1;
  double * data2;
  double * data3;

  double ** gaussColor; // each element receives an array of gauss points
  double * gaussColorPtr; // for colorRange, treat as single array
  double ** devGaussColor; // each element receives an array of gauss points
  double * devGaussColorPtr; // for colorRange, treat as single array
  double ** stiffEigenColor;
  double * stiffEigenColorPtr;

  int numGaussPts;

  DrawGaussPoints * gsdraw; // DrawGaussPoints hold state, thus 1 per element
  DrawStiffness * stfdraw; // same

  int currRange; //current color map
  double colorMax, colorMin;
  bool invert;

  double center[3];  //bounding box center
  double compassSize; // X,Y,Z legend proportional to volume
  double compassLoc[3]; //origin of compass, relative to 0,0,0

  DrawElement * drawElementTypes[10]; //these are singleton instances of each 
                                   //type of DrawElement object
                                    //NOT REQUIRING STATE

  DrawElement ** drawMethods; //one drawMethod per element
  File_Channel * theFileChannel; //update state

  ReynoldsGlyph * rg;
  
};
static char * drawEleNames[] = { "Mesh",  "Displaced Mesh",
				  "Solid","Displaced Solid", 
                                "Gauss Points","Stiffness" };
#define numDrawEleNames 6


#endif
