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
** DrawSolid: draw an element as a solid (unless it is 2node, then all you get
** is a line). 
**/

#include "materials.h"
#include "auxiliary.h"
#include "DrawSolid.h"
#include "Camera.h"
DrawSolid::DrawSolid() {}
 
void  DrawSolid::drawElement(Element * theElement, double * theColors, 
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
  
void DrawSolid::drawTwoNode( Element * theElement, double * theColors, int colorflag ) {

   Node ** theNodes; // nodes for a single element
  int numNodes,i,pick;
  Vector v;
  Vector v1;
 
 
  theNodes = theElement->getNodePtrs();
  numNodes = theElement->getNumExternalNodes(); 

  glGetIntegerv(GL_RENDER_MODE,&pick); 
  if( pick == GL_SELECT )
	glPushName( theElement->getTag() ); //make element pickable


  /* for now this is just a line: should it be a tube? */
 
   glBegin( GL_LINES);
	  for ( i = 0; i < numNodes; i++ ) {
		 v = theNodes[i]->getCrds();

		if ( colorflag == BY_NODE ) {
		  set_color_by_param( theColors[theNodes[i]->getTag()], 1.0 );
		}
		else if ( colorflag == BY_ELEMENT ) {
		  set_color_by_param( theColors[theElement->getTag()], 1.0 );
		}	 
		if (v.Size() == 3 ) {
		  glVertex3f(v[0], v[1],v[2]);
		  glNormal3f(v[0], -v[2], v[1]);
		}
		else if (v.Size() == 2 )
		  glVertex3f(v[0], v[1],0 );
		
	  }
   glEnd();



  glGetIntegerv(GL_RENDER_MODE,&pick); 
  if( pick == GL_SELECT )
	glPopName(); 




}

void DrawSolid::drawFourNode( Element * theElement, double * theColors, int colorflag ) {
  Node ** theNodes; // nodes for a single element
  int numNodes,i,pick;
  Vector v,v2,d,d2;
  double norm[3];  /* normal calculation */
  double pts[3][3];   /* vector to array transform */

   theNodes = theElement->getNodePtrs();
   numNodes = theElement->getNumExternalNodes(); 


 
  //calculate normal for quad
  for ( i = 0; i < 3; i++ ) {
	v = theNodes[i]->getCrds();
	pts[i][0] = v[0];
	pts[i][1] = v[1];
	if ( v.Size() == 3 )
	  pts[i][2] = v[2];
  }
  calcNormal(pts[0], pts[1], pts[2], norm);

  
 glGetIntegerv(GL_RENDER_MODE,&pick); 
 if( pick == GL_SELECT ) {
	glPushName( theElement->getTag() ); //make element pickable
	printf("push %d\n",theElement->getTag());
 }
  glBegin( GL_QUADS );
	  for ( i = 0; i < numNodes; i++ ) {
		v = theNodes[i]->getCrds();
	 		
		if ( colorflag == BY_NODE ) {
		  set_color_by_param( theColors[theNodes[i]->getTag()], 1.0 );
		}
		else if ( colorflag == BY_ELEMENT ) {		
		  set_color_by_param( theColors[theElement->getTag()], 1.0 );
		}
		
		if (v.Size() == 3 ) {
		  glNormal3f(norm[0], norm[1],norm[2]);
		  glVertex3f(v[0], v[1],v[2]);
		 
		}
		else if (v.Size() == 2 ) {
		  glVertex3f(v[0], v[1],0 );
		  
		}
		else  glVertex3f(v[0],0,0); 
		
	  }
  glEnd();


  glGetIntegerv(GL_RENDER_MODE,&pick); 
  if( pick == GL_SELECT )
	glPopName(); 
} // 4 nodes such as a quad element


void DrawSolid::drawEightNode( Element * theElement, double * theColors, int colorflag ) {	
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
		  pts[j][X] = v[X] ;
		  pts[j][Y] = v[Y] ;
		  pts[j][Z] = v[Z] ;
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
			glVertex3f( v[X] , v[Y], v[Z] );
		}
	  } // 4 sides of parallelpiped
	  
	  glEnd();

	  
	  //draw top quad.

		for ( j = 0; j < 3; j++) {
		  v = theNodes[j]->getCrds();
		  d = theNodes[j]->getDisp();
		  pts[j][X] = v[X] ;
		  pts[j][Y] = v[Y] ;
		  pts[j][Z] = v[Z] ;
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
		 glVertex3f( v[X] , v[Y], v[Z] );
	   }
   glEnd();


   //draw bottom quad
   		for ( j = 4; j < 7; j++) {
		  v = theNodes[j]->getCrds();
		  d = theNodes[j]->getDisp();
		  pts[j-4][X] = v[X] ;
		  pts[j-4][Y] = v[Y] ;
		  pts[j-4][Z] = v[Z] ;
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
		 glVertex3f( v[X] , v[Y], v[Z] );
	   }
   glEnd();

	glGetIntegerv(GL_RENDER_MODE,&pick); 

  if( pick == GL_SELECT )
	glPopName(); 
} // 8 node element - some type of brick

  



/** We can draw any type of element as a solid */
bool  DrawSolid::canDraw( Element * theElement ){
  int   numNodes = theElement->getNumExternalNodes(); 
  if ( numNodes == 2 )
	return true;
  if ( numNodes == 4 )
	return true;
  if ( numNodes == 8 )
	return true;
  
  return false;
}

int DrawSolid::getClassTag() {
  return DRAW_SOLID;
}
