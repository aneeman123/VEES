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

#include "Viewer.h"

Viewer::Viewer() {}

void Viewer::draw(){
  printf("Error in code: Viewer.draw() should be implemented in subclass\n");
}

double Viewer::getCenterAndMaxDim( double center[3] )
{
  printf("Error in code: Viewer.getCenterAndMaxDim should be implemented in subclass\n");
  return 0;
}

int Viewer::getNumElements()
{
  printf("Error in code: Viewer.getNumElements() should be implemented in subclass\n");
  return 0;
}

void  Viewer::setSelected(int i) {
  printf("Error in code: Viewer.setSelected(int i) should be implemented in subclass\n");
}
