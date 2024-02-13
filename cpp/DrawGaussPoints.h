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
** DrawGaussPoints.cpp is a class to draw a finite element's
** integration points. .
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

#ifndef _DRAW_GAUSS_PTS_
#define _DRAW_GAUSS_PTS_


class DrawGaussPoints:public DrawElement{
 public:
  DrawGaussPoints();
  ~DrawGaussPoints();

  void drawElement(Element * el, double * colors, 
		   ColorRange range, int colorStyleFlag );
  bool canDraw( Element * el );
  int getClassTag();
  int setElement( Element * ele );
  void setPoints( const Matrix & other, double radius );
  
  int getNumPts();
  const Matrix &  getMatrix();
  double getSphereSize();
  
  void setSelected(int gaussNum); // highlight chosen gauss point

  // associated data
  bool setMeanStress( double * data );
  void setDeviatoricStress( double * data );
 private:
  
  double  Mdistance( Matrix pts, int j, int k );

  GLUquadricObj* sphere; // one sphere, drawn many times. 
  int numPts;
  Matrix points; //numGaussPoints * 3 matrix
  double sphereSize;
  Element * el;
  bool drawByMatrixAlone;
  int selected; // gauss pt selected for viewing in another viewer
};

#endif
