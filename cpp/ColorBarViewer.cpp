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
#include <ColorBarViewer.h>
ColorBarViewer::ColorBarViewer() {
  theRange = NULL;
}
void ColorBarViewer::setRange( ColorRange * range ){
  theRange = range;
}

void ColorBarViewer::draw() {
  int i;

  glBegin(GL_QUAD_STRIP);
  for( i = 0; i < BAR_LEN; i++ ) {
	
	set_color_by_param( (float) i/BAR_LEN , 1.0 );
	glNormal3f(0,0,1);
	glVertex3f( i*3 - BAR_LEN, -7,0); // x range: 10 to -20 
	glVertex3f( i*3 - BAR_LEN, 7,0); // y range 0 to 5
	
  }
  glEnd();

}

/** bounding box is 0 to BAR_LEN along x, 
 * 0 to BAR_HEIGHT on y, 0 on z */ 
double ColorBarViewer::getCenterAndMaxDim( double center[3] ) {
  center[X] = BAR_LEN/2;
  center[Y] = BAR_HEIGHT/2;
  center[Z] = 0;

  return BAR_LEN;
}


int ColorBarViewer::getNumElements() {
  return 1;
}
