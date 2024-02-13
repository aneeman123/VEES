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
** Viewer is Superclass to DomainViewer, ElementViewer and GaussViewer.
** This class is passed to a VisWin who will need to ask it
** to draw itself, and reveal its bounding box and center point
** needed to set camera & projection.
*/
#ifndef _VIEWER_
#define _VIEWER_

#include <stdio.h>
class Viewer {
 public:
  Viewer();
  virtual void draw(); 
  virtual double getCenterAndMaxDim( double center[3] );
  virtual int getNumElements();
  virtual void setSelected(int i);

};

#endif
