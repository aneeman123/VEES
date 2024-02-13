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
** PointViewer draws a single element. Integration points or individual
** fibers can be selected
**/

#include "PointViewer.h"

/** default constructor */
PointViewer::PointViewer(){
  eleNum = -1;
  theDomain = NULL;
  theElement = NULL;
}

void PointViewer::setDomain( Domain * d ) {
  theDomain = d;

  drawElementTypes[DRAW_WIRE] = new DrawWireFrame();
  drawElementTypes[DRAW_WIRE_DISPLACED] = new DrawDisplacedWireFrame();
  drawElementTypes[DRAW_SOLID] = new DrawSolid();
  drawElementTypes[DRAW_SOLID_DISPLACED] = new DrawSolidDisplaced();
  drawElementTypes[DRAW_GAUSS_PTS] = new DrawGaussPoints();
  /* we don't need HideElement */
}

int PointViewer::getNumElements() {
  if (theElement == NULL )
	return 0;
  return 1;
}


void PointViewer::draw() {
  int pick;
  double color[1];
  color[0]= 0.0;

  ColorRange range;
  // save pick state, only enable during gauss point draw 
  glGetIntegerv(GL_RENDER_MODE,&pick); 

  if (theElement != NULL ) {

	if( pick != GL_SELECT )  // only render wireframe for non-pick
		drawElementTypes[DRAW_WIRE]->drawElement( 
								 theElement,color, range, SINGLE_COLOR );

 	  drawElementTypes[DRAW_GAUSS_PTS]->drawElement( 
								 theElement,color, range, SINGLE_COLOR );  

  }

  //   else {
		  //	setElement( 1 );
	  //drawElementTypes[DRAW_WIRE]->drawElement( 
	  //					 theElement,color, range, SINGLE_COLOR ); 
  // }

} 

/**
 * get center and max dimensions of a single finite element.
 * @pre! element number has been set
 */
double PointViewer::getCenterAndMaxDim( double center[3] ) {

  Node ** theNodes; // nodes for a single element
  Node * theNode;
  Vector v;
  int i,j,numNodes, numCrds;

  if ( eleNum == -1 ) {
	printf( "Error getting center and bounding box; no elemnt selected\n");
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

  for (j = 1; j < numNodes; j++ ){	  
	v = theNodes[j]->getCrds(); 
	for( i = 0; i < numCrds; i++ ) {
	  if ( max[i] < v[i] )
		max[i] = v[i];
	  if ( min[i] > v[i] )
		min[i] = v[i];
	} 
  } // while more nodes
  
	// center is average of min and max (bounding box)
  for (i = 0; i < 3; i++ )
	center[i] = (max[i] + min[i])/2;
  
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
void PointViewer::setElement( int num ) {
  printf("rec'd %d\n", num);
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
  }// actual 
}
