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

#include "ElementViewer.h"

/** default constructor */
ElementViewer::ElementViewer(){
  eleNum = -1;
  theDomain = NULL;
  theElement = NULL;
}

void ElementViewer::setDomain( Domain * d ) {
  theDomain = d;

  drawElementTypes[DRAW_WIRE] = new DrawWireFrame();
  drawElementTypes[DRAW_WIRE_DISPLACED] = new DrawDisplacedWireFrame();
  drawElementTypes[DRAW_SOLID] = new DrawSolid();
  drawElementTypes[DRAW_SOLID_DISPLACED] = new DrawSolidDisplaced();
  drawElementTypes[DRAW_GAUSS_PTS] = new DrawGaussPoints();
  /* we don't need HideElement */
}

int ElementViewer::getNumElements() {
  if (theElement == NULL )
	return 0;
  return 1;
}

void ElementViewer::setSelected( int i ) {

  ((DrawGaussPoints *)drawElementTypes[DRAW_GAUSS_PTS])->setSelected(i);
}

void ElementViewer::draw() {
  int pick;
  double color[1];
  color[0]= 0.0;

  ColorRange range;
  // save pick state, only enable during gauss point draw 
  glGetIntegerv(GL_RENDER_MODE,&pick); 

  if (theElement != NULL ) {
	//	DrawPlaneInABox * dpb = new	DrawPlaneInABox( theElement );  //test
	if( pick != GL_SELECT )  // only render wireframe for non-pick
	  drawElementTypes[DRAW_WIRE]->drawElement( 
				 theElement,color, range, SINGLE_COLOR );
	

	drawElementTypes[DRAW_GAUSS_PTS]->drawElement( 
	       	 theElement,color, range, SINGLE_COLOR );  
	//dpb->drawElement( theElement,color, range, SINGLE_COLOR );
  }
  drawCompass();
} 

/**
 * get center and max dimensions of a single finite element.
 * @pre! element number has been set
 */
double ElementViewer::getCenterAndMaxDim( double center[3] ) {

  Node ** theNodes; // nodes for a single element
  Node * theNode;
  Vector v;
  int i,j,numNodes, numCrds;

  if ( eleNum == -1 ) {
	center[0] = center[1] = center[2] = 0;
	return 0;
  }
  theElement = theDomain->getElement( eleNum );

  if (theElement == NULL ) {
	printf("Error, no element with number %d\n", eleNum);
	center[0] = center[1] = center[2] = 0;
	return 0.0;
  }
  
  theNodes = theElement->getNodePtrs();
  numNodes = theElement->getNumExternalNodes();

  //make the first node both min and max.
  //go through the rest of the nodes and get bigger or smaller x, y, z values
  // if no z values, just use zero for z. 
  if ( theNodes == NULL ) {
	printf("Error fetching element nodes\n");
	return 0.0;
  }
  theNode = theNodes[0];
  v = theNode->getCrds(); 
  numCrds = v.Size();

  // baseline from first node
  for(i = 0; i < numCrds; i++ )
	max[i] = min[i] = v[i];

  if ( numCrds < 3 )
	max[2] = min[2] = 0;

  center[0] = center[1] = center[2] = 0;

  for (j = 0; j < numNodes; j++ ){	  

	v = theNodes[j]->getCrds(); 

	for( i = 0; i < numCrds; i++ ) {

	  center[i] = center[i] + v[i];
	
	  if ( max[i] < v[i] )
		max[i] = v[i];
	  if ( min[i] > v[i] )
		min[i] = v[i];
	} 

  } // while more nodes

  //find average of all nodes (unweighted centroid :)
  for ( i = 0; i < numCrds; i++ ) {
	center[i] = center[i]/numNodes;
  }
  // biggest absolute diff in bounding box
  maxDim = 0;
  for (i = 0; i < numCrds; i++ ) {
	if ( fabs(max[i] - min[i]) > maxDim )
	  maxDim =  fabs(max[i] - min[i]);
  }

  return maxDim; 

}

/**
 * Set the element to be viewed */
void ElementViewer::setElement( int num ) {
 
  if ( theDomain == NULL ) {
	printf( "Error : domain has not been set so element cannot be shown\n");
	return;
  }
  eleNum = num;

  if ( eleNum <= -1 )
	theElement = NULL;

  else {
	theElement = theDomain->getElement( eleNum );
	
	if (theElement == NULL ) {
	  printf("Error (setElement), no element with number %d\n", eleNum);
	  return;
	}
	((DrawGaussPoints *)drawElementTypes[DRAW_GAUSS_PTS])->setElement( theElement );
   

  }// actual 
}

int  ElementViewer::getElement() {
  return eleNum;
}

/**
 * Draw X, Y, Z axes in red, green blue at lowest x,y,z, corner
 */
void  ElementViewer::drawCompass() {
  double base,height, len;
  double compassSize = 1.3;

  GLUquadric * cylinder;
  cylinder = gluNewQuadric();
  glLineWidth(2);

  height = compassSize;
  base = compassSize/8;

  len = height/3;
  set_color_by_param_shiny( 1.0, 1.0 ); //red
  // set_saturation_by_param(0.0, 1.0); // black
  Vector v = theDomain->getPhysicalBounds();
  glPolygonMode(GL_FRONT, GL_FILL); 

  //------- Z axis ---------
   set_color_by_param_shiny( 0, 1.0 ); //blue
  glPushMatrix();
  
  //  glTranslatef( compassLoc[X], compassLoc[Y], compassLoc[Z]);
   gluCylinder( cylinder , base/7 , base/7, 
	       height*.85, 8, 1 );
   
 glPushMatrix();
  glTranslatef(0,0, height*.85);
  glutSolidCone(base/3, base/2,6,1);
  glPopMatrix();

 
  glPopMatrix();
  
  //---------- X axis ---------
  set_color_by_param_shiny( 1.0, 1.0 ); //red
  glPushMatrix();
  //  glTranslatef( compassLoc[X], compassLoc[Y], compassLoc[Z]);
  
  glRotatef(90, 0.0, 1.0, 0.0); //rotate around Y axis, swap x & z
  gluCylinder( cylinder , base/7 , base/7, 
	       height*.85, 8 , 1);
  glPushMatrix();
  glTranslatef(0,0, height*.85);
  glutSolidCone(base/3, base/2,6,1);
  glPopMatrix();
 

  glPopMatrix();

  //------- Y axis --------
  
   set_color_by_param_shiny( 0.5, 1.0 ); //green

  glPushMatrix();
  //  glTranslatef( compassLoc[X], compassLoc[Y], compassLoc[Z]);
  
  glRotatef(270, 1.0, 0.0, 0.0); //rotate around X axis
  gluCylinder( cylinder , base/7 , base/7, 
	       height*.85 , 8 , 1 );

  glPushMatrix();
  glTranslatef(0, 0,height*.85 );
  glutSolidCone(base/3, base/2,6,1);
  glPopMatrix();

 
  glPopMatrix();
  set_color_by_param_shiny(1.0,1.0); // change params re lighting
  
  delete cylinder; // clean up
  glLineWidth(1);
}
