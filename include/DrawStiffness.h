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
** DrawStiffness.cpp is a class to draw a finite element's
** stiffness at integration points. .
**/

#include "DrawElement.h"
#include <OPS_Stream.h>
#include <StandardStream.h>
#include <Recorder.h>
#include <ElementRecorder.h>
#include <Information.h>
#include <Response.h>
#include <Matrix.h>

#ifdef WIN32
#include <windows.h>
#endif

#ifndef _DRAW_STIFFNESS_
#define _DRAW_STIFFNESS_

#include "ReynoldsGlyph.h" // static arrays of triangles
#include "Camera.h"
#include "auxiliary.h"
#include "ColorSchemeNames.h"
//#define NUM_TRI_SMALL 320
//#define NUM_PTS_SMALL 312


class DrawStiffness:public DrawElement{
 public:
  DrawStiffness();
  ~DrawStiffness();

  void drawElement(Element * el, double * colors, 
		   ColorRange range, int colorStyleFlag );
  bool canDraw( Element * el );
  int getClassTag();
  double getMinDistanceBetweenPts();
  int setElement( Element * ele );
  void setPoints( const Matrix & other, double radius );
  
  int getNumPts();
  const Matrix &  getMatrix();
  double getSphereSize();
  
  void setToStress();    // draw stress 
  void setToStiffness(); // draw low stiffness eigentensor 

  // associated data
  bool setMeanStress( double * data );
  void setMinStiffness( double * data );

  bool insideElement( double thePoint[3] );
  void setFacetSign();
  double  orient3D( double a[3], double b[3], double c[3], double d[3] );

  // print to File
  void printStiffness( FILE * theFile ); // matrices
  void printPointsVTK( FILE * theFile ); //gauss points
  void printMinStiffVTK( FILE * theFile ); 
  void printEigenmodeVTK( FILE * theFile );// scalar, 0 to 1 derived from RGB
  void printRotationTest(FILE * theFile);
  double matrixDifference(Matrix m_A, Matrix m_B);
  void normalizeMatrix( Matrix & theMatrix );
  void interpolateStiffness(Matrix & m_stiff, double loc[3]);
  void interpolateStress(Matrix & m_stiff, double loc[3]);
  void interpolateEigenTensor(Matrix & m_eigen, double loc[3]);

 private:
 
  ReynoldsGlyph* rg;
  // how about eight GlyphViewers? The problem is the resolution is too high.
  // I could slash and burn them.

  double  Mdistance( Matrix pts, int j, int k );

  GLUquadricObj* sphere; // one sphere, drawn many times. 
  int numPts;
  Matrix points; //numGaussPoints * 3 matrix
  double sphereSize;
  double minDistBetweenPts;
  Element * el;
  bool drawByMatrixAlone;
  bool selectable; // element selected for viewing in another viewer

  //---- inside/outside for interpolation ---- 
  double facetSign[6]; // normal points in or out of polytope
 
};


static int node2Gauss[8] = // lookup table for closest node to gauss location
  {6,2,5,1,7,3,4,0}; 

static int gauss2Node[8] = // get gauss point in node order
  {7,3,1,5,6,2,0,4};       

#endif
