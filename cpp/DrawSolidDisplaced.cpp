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
** DrawSolidDisplaced.cpp : draw a finite element as a solid (except for 2Node 
** element, these will simply be lines) with nodes moved by committed
** displacement. 
**/

#include "materials.h"
#include "auxiliary.h"
#include "DrawSolidDisplaced.h"
#include "Camera.h"
DrawSolidDisplaced::DrawSolidDisplaced() {}
 
void  DrawSolidDisplaced::drawElement(Element * theElement, double * theColors, 
						   ColorRange range, int flag )
{
  int   numNodes = theElement->getNumExternalNodes(); 

  //check for within range
  if ( flag == BY_NODE ) { 
	if ( !range.nodesInRange( theColors, theElement ) )
	  return;
  }
 
	  
  else if ( flag == BY_ELEMENT ) {
	if ( ! range.inRange(theColors[theElement->getTag()]) )
	  return;
  }

  if ( flag == SINGLE_COLOR )
	set_color_by_param( theColors[0], 1.0);

  if ( numNodes == 2 ) 
	drawTwoNode( theElement, theColors, flag );

  else if( numNodes == 4 ) 
	drawFourNode( theElement, theColors, flag );

  else if( numNodes == 8 ) 
	drawEightNode( theElement, theColors,  flag );

}
  
void DrawSolidDisplaced::drawTwoNode( Element * theElement, double * theColors, int colorflag ) {

  Node ** theNodes; // nodes for a single element
  int numNodes,i,pick;
  Vector v;
  Vector d;

  theNodes = theElement->getNodePtrs();
  numNodes = theElement->getNumExternalNodes(); 

  glGetIntegerv(GL_RENDER_MODE,&pick); 
  if( pick == GL_SELECT )
	glPushName( theElement->getTag() ); //make element pickable


  /* for now this is just a line: should it be a tube? */
 
   glBegin( GL_LINES);
	  for ( i = 0; i < numNodes; i++ ) {
		 v = theNodes[i]->getCrds();
		 d = theNodes[i]->getDisp();
		if ( colorflag == BY_NODE ) {
		  set_color_by_param( theColors[theNodes[i]->getTag()], 1.0 );
		}
		else if ( colorflag == BY_ELEMENT ) {
		  set_color_by_param( theColors[theElement->getTag()], 1.0 );
		}	 
		if (v.Size() == 3 )
		  glVertex3f(v[X]+d[X], v[Y]+d[Y],v[Z]+d[Z]);
		else if (v.Size() == 2 )
		  glVertex3f(v[X]+d[X], v[Y]+d[Y],0 );
		
	  }
   glEnd();



  glGetIntegerv(GL_RENDER_MODE,&pick); 
  if( pick == GL_SELECT )
	glPopName(); 
}

void DrawSolidDisplaced::drawFourNode( Element * theElement, double * theColors, int colorflag ) {
  Node ** theNodes; // nodes for a single element
  int numNodes,i,pick;
  Vector v,v2,d,d2;
  double facetNorm[3];  /* normal calculation */
  double pts[3][3];   /* vector to array transform */

   theNodes = theElement->getNodePtrs();
   numNodes = theElement->getNumExternalNodes(); 


 
  //calculate normal for quad
  for ( i = 0; i < 3; i++ ) {
	v = theNodes[i]->getCrds();
	d = theNodes[i]->getDisp();
	pts[i][X] = v[X] + d[X];
	pts[i][Y] = v[Y] + d[Y];
	if ( v.Size() == 3 )
	  pts[i][Z] = v[Z]+d[Z];
  }
  calcNormal(pts[X], pts[Y], pts[Z], facetNorm);

  
 glGetIntegerv(GL_RENDER_MODE,&pick); 
 if( pick == GL_SELECT ) {
	glPushName( theElement->getTag() ); //make element pickable
	
 }
  glBegin( GL_QUADS );
	  for ( i = 0; i < numNodes; i++ ) {
		v = theNodes[i]->getCrds();
		d = theNodes[i]->getDisp();		
		if ( colorflag == BY_NODE ) {
		  set_color_by_param( theColors[theNodes[i]->getTag()], 1.0 );
		}
		else if ( colorflag == BY_ELEMENT ) {		
		  set_color_by_param( theColors[theElement->getTag()], 1.0 );
		}
		
		if (v.Size() == 3 ) {
		  glNormal3f(facetNorm[X], facetNorm[Y],facetNorm[Z]);
		  glVertex3f(v[X] + d[X], v[Y] + d[Y],v[Z]+d[Z]);

		}
		else if (v.Size() == 2 ) {
		  glVertex3f(v[X] + d[X], v[Y] + d[Y],0 );
		  
		}
		else  glVertex3f(v[X] + d[X],0,0); 
	
	  }
  glEnd();


  glGetIntegerv(GL_RENDER_MODE,&pick); 
  if( pick == GL_SELECT )
	glPopName(); 
} // 4 nodes such as a quad element


void DrawSolidDisplaced::drawEightNode( Element * theElement, double * theColors, int colorflag ) {	
  Node ** theNodes; // nodes for a single element
  int numNodes,i,j,pick;
  Vector v,v2,d,d2;
  double facetNorm[3];  /* normal calculation */
  double pts[3][3];   /* vector to array transform */

  theNodes = theElement->getNodePtrs();
   numNodes = theElement->getNumExternalNodes(); 

  for (i = 0; i < 3; i++)
	facetNorm[i] = 0;
 
  glGetIntegerv(GL_RENDER_MODE,&pick); 
  if( pick == GL_SELECT )
	glPushName( theElement->getTag() ); //make element pickable



  // draw 4 sides of parallelpiped
  // this should work for 8n brick
  //  set_color_by_param(0.5, 1.0);
	
	  glBegin(GL_QUADS);
	  for ( i = 0; i < 4; i++ ) {
	
		// entire quad will have same normals 
		//if (i != 0)
		for ( j = 0; j < 3; j++ ) {
		  v = theNodes[brick8Loops[i][j]]->getCrds();
		  d = theNodes[brick8Loops[i][j]]->getDisp();
		  pts[j][X] = v[X] + d[X];
		  pts[j][Y] = v[Y] + d[Y];
		  pts[j][Z] = v[Z] + d[Z];
		}

		calcNormal( pts[X], pts[Y], pts[Z], facetNorm);

		for( j = 0; j < 4; j++ ) {
			
		    v = theNodes[brick8Loops[i][j]]->getCrds();
		    d = theNodes[brick8Loops[i][j]]->getDisp();	
			if ( colorflag == BY_NODE ) {
			  set_color_by_param( theColors[theNodes[brick8Loops[i][j]]->getTag()], 1.0 );
			}
			else if ( colorflag == BY_ELEMENT ) {
			  set_color_by_param( theColors[theElement->getTag()], 1.0 );
			}

			glNormal3f( facetNorm[X], facetNorm[Y], facetNorm[Z] );
			glVertex3f( v[X] + d[X], v[Y]+ d[Y], v[Z] +d[Z]);
		}
	  } // 4 sides of parallelpiped
	  
	  glEnd();

	  
	  //draw top quad.

		for ( j = 0; j < 3; j++) {
		  v = theNodes[j]->getCrds();
		  d = theNodes[j]->getDisp();
		  pts[j][X] = v[X] + d[X];
		  pts[j][Y] = v[Y] + d[Y];
		  pts[j][Z] = v[Z] + d[Z];
		}
		calcNormal( pts[X], pts[Y], pts[Z], facetNorm);

   glBegin( GL_QUADS );	// not sure why top & bottom have same normal
	   for( i = 0; i < 4; i++ ){ 
		 v = theNodes[i]->getCrds();
		 d = theNodes[i]->getDisp();

		 if ( colorflag == BY_NODE ) {
		   set_color_by_param( theColors[theNodes[i]->getTag()], 1.0 );
		 }
		 else if ( colorflag == BY_ELEMENT ) {
		   set_color_by_param( theColors[theElement->getTag()], 1.0 );
		 }
	
		 glNormal3f( facetNorm[X], facetNorm[Y], facetNorm[Z] );
		 glVertex3f( v[X] + d[X], v[Y]+ d[Y], v[Z]+d[Z] );
	   }
   glEnd();


   //draw bottom quad
   		for ( j = 4; j < 7; j++) {
		  v = theNodes[j]->getCrds();
		  d = theNodes[j]->getDisp();
		  pts[j-4][X] = v[X] + d[X];
		  pts[j-4][Y] = v[Y] + d[Y];
		  pts[j-4][Z] = v[Z] + d[Z];
		}
		calcNormal( pts[X], pts[Y], pts[Z], facetNorm);

	 glBegin( GL_QUADS );	
	   for( i = 4; i < 8; i++ ){ 
		 v = theNodes[i]->getCrds();
		 d = theNodes[i]->getDisp();

		 if ( colorflag == BY_NODE ) {
		   set_color_by_param( theColors[theNodes[i]->getTag()], 1.0 );
		 }
		 else if ( colorflag == BY_ELEMENT ) {
		   set_color_by_param( theColors[theElement->getTag()], 1.0 );
		 }
	
		 glNormal3f( facetNorm[X], facetNorm[Y], facetNorm[Z] );
		 glVertex3f( v[X] + d[X], v[Y]+ d[Y], v[Z]+d[Z] );
	   }
   glEnd();

	glGetIntegerv(GL_RENDER_MODE,&pick); 

  if( pick == GL_SELECT )
	glPopName(); 
 
} // 8 node element - some type of brick

  



/** We can draw any type of element as a solid */
bool  DrawSolidDisplaced::canDraw( Element * theElement ){
  int   numNodes = theElement->getNumExternalNodes(); 
  if ( numNodes == 2 )
	return true;
  if ( numNodes == 4 )
	return true;
  if ( numNodes == 8 )
	return true;
  
  return false;
}

int DrawSolidDisplaced::getClassTag() {
  return DRAW_SOLID_DISPLACED;
}
