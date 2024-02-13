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
** DrawWireFrame.cpp : draw a finite element as a wireframe. Note: the element
** will still be pickable as if it is a solid.
**/

#include "DrawWireFrame.h"
#include "materials.h"
#include "DrawSolid.h"
DrawWireFrame::DrawWireFrame() {}

void DrawWireFrame::drawElement(Element * theElement, double * theColors, 
						   ColorRange range, int flag )
{
  int   numNodes = theElement->getNumExternalNodes(); 
  int count = 0;
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
	set_saturation_by_param( theColors[0], 1.0);
		

  if ( numNodes == 2 ) 
	drawTwoNode( theElement, theColors, flag );

  else if( numNodes == 4 ) 
	drawFourNode( theElement, theColors, range, flag );

  else if( numNodes == 8 ) 
	drawEightNode( theElement, theColors,   range,flag );
  else if (numNodes == 20)
	drawTwentyNode( theElement, theColors,   range,flag );
}

/* We can draw any type of element as wireframe */
bool  DrawWireFrame::canDraw( Element * theElement )
{
 int   numNodes = theElement->getNumExternalNodes(); 
 if ( numNodes == 2 )
  return true;
 if ( numNodes == 4 )
  return true;
 if ( numNodes == 8 )
  return true;
 if ( numNodes == 20 )
  return true;

 return false;
}

void DrawWireFrame::drawTwentyNode( Element * theElement, double * theColors, ColorRange range, int colorflag )
{
  Node ** theNodes; // nodes for a single element
  int numNodes,i;
  DrawSolid dts;
  Vector v,v2;

  int lines[6][2]={{1,2},{3,4}, {1,2},{3,4},{1,2},{3,4}};
  
  theNodes = theElement->getNodePtrs();
  numNodes = theElement->getNumExternalNodes(); 

	  //draw top and bottom quads.
  glBegin( GL_LINES );	
    for( i = 0; i < 6; i++ ){
	  v = theNodes[lines[i][0]]->getCrds();  
	
 	  glVertex3f(v[0], v[1],v[2]);
	  v = theNodes[lines[i][1]]->getCrds();  
	
 	  glVertex3f(v[0], v[1],v[2]);
  }
  glEnd();
  /*
  glBegin( GL_LINE_LOOP );	
	 for( i = 0; i < 8; i++ ){ 
	   v = theNodes[zFaceBack[i]]->getCrds();
	   
	   if ( colorflag == BY_NODE )
		 set_color_by_param( theColors[theNodes[zFaceBack[i]]->getTag()], 1.0 );
	   else if ( colorflag == BY_ELEMENT )
		 set_color_by_param( theColors[theElement->getTag()], 1.0 );
	   
	   glVertex3f(v[0], v[1],v[2]);
	  }
  glEnd();
  */
  /*
	for(i = 0; i < 4; i++) {
	glBegin( GL_LINE_STRIP );
	   for(j = 0; j < 3; j++) {
	   // color
	   //line
	   v = theNodes[sides[i][j]]->getCrds();

	   glVertex3f(v[0], v[1],v[2]);
	   }
	   glEnd();
	}
  */

}

void DrawWireFrame::drawEightNode( Element * theElement, double * theColors, ColorRange range, int colorflag )
{
  Node ** theNodes; // nodes for a single element
  int numNodes,i, pick;
  DrawSolid dts;
  Vector v,v2;

  theNodes = theElement->getNodePtrs();
  numNodes = theElement->getNumExternalNodes(); 

  glGetIntegerv(GL_RENDER_MODE,&pick); 
  if( pick == GL_SELECT ) {
	
	dts.drawElement(theElement, theColors, range, colorflag);

  }

  else {

		 glBegin( GL_LINES);	
		  for( i = 0; i < 4; i++ ){ // draw vertical bars
			v = theNodes[i]->getCrds();
	  		v2 = theNodes[i+4]->getCrds();

			if ( colorflag == BY_NODE ) {
				set_color_by_param( theColors[theNodes[i]->getTag()], 1.0 );
			}
			else if ( colorflag == BY_ELEMENT ) {
			set_color_by_param( theColors[theElement->getTag()], 1.0 );
			}

			glVertex3f( v[0], v[1], v[2] );

			if ( colorflag == BY_NODE )
			set_color_by_param( theColors[theNodes[i+4]->getTag()], 1.0 );
			else if ( colorflag == BY_ELEMENT )
			set_color_by_param( theColors[theElement->getTag()], 1.0 );
		

	   		glVertex3f( v2[0], v2[1], v2[2] );	  
		  }
    	glEnd();
	  //draw top and bottom quads.
		glBegin( GL_LINE_LOOP );	
		for( i = 0; i < 4; i++ ){ 
			v = theNodes[i]->getCrds();
	
		if ( colorflag == BY_NODE )
			  set_color_by_param( theColors[theNodes[i]->getTag()], 1.0 );
		else if ( colorflag == BY_ELEMENT )
			  set_color_by_param( theColors[theElement->getTag()], 1.0 );
	
			glVertex3f(v[0], v[1],v[2]);
		  }	
		glEnd();
	  
	glBegin( GL_LINE_LOOP );	
	  for( i = 4; i < 8; i++ ){ 
		v = theNodes[i]->getCrds();

		if ( colorflag == BY_NODE )
		  set_color_by_param( theColors[theNodes[i]->getTag()], 1.0 );
		else if ( colorflag == BY_ELEMENT )
		  set_color_by_param( theColors[theElement->getTag()], 1.0 );
	
		glVertex3f(v[0], v[1],v[2]);
	  }
  glEnd();
  }
}

void DrawWireFrame::drawFourNode( Element * theElement, double * theColors, ColorRange range, int colorflag )
{
  Node ** theNodes; // nodes for a single element
  int numNodes,i,pick;
  DrawSolid dts;
  Vector v,v2;

  theNodes = theElement->getNodePtrs();
  numNodes = theElement->getNumExternalNodes(); 

  glGetIntegerv(GL_RENDER_MODE,&pick); 
  if( pick == GL_SELECT ) {
	dts.drawElement(theElement, theColors, range, colorflag);
	return;
  }

  glBegin( GL_LINE_LOOP );
	  for ( i = 0; i < numNodes; i++ ) {
		v = theNodes[i]->getCrds();

		if ( colorflag == BY_NODE )
		  set_color_by_param( theColors[theNodes[i]->getTag()], 1.0 );
		else if ( colorflag == BY_ELEMENT )
		  set_color_by_param( theColors[theElement->getTag()], 1.0 );

		

		if (v.Size() == 3 ) {
		  glVertex3f(v[0], v[1],v[2]);
		}
		else if (v.Size() == 2 ) {
		  glVertex3f(v[0], v[1],0 );
		
		}
		else  glVertex3f(v[0],0,0); 
	  }
  glEnd(); 

}

void DrawWireFrame::drawTwoNode( Element * theElement, double * theColors, int colorflag )
{  Node ** theNodes; // nodes for a single element
  int numNodes,i,pick;
  DrawSolid dts;
  Vector v,v2;

  theNodes = theElement->getNodePtrs();
  numNodes = theElement->getNumExternalNodes(); 

 glGetIntegerv(GL_RENDER_MODE,&pick); 
  if( pick == GL_SELECT )
	glPushName( theElement->getTag() ); //make element pickable

  glBegin( GL_LINES);
	  for ( i = 0; i < numNodes; i++ ) {

		if ( colorflag == BY_NODE )
		  set_color_by_param( theColors[theNodes[i]->getTag()], 1.0 );
		else if ( colorflag == BY_ELEMENT )
		  set_color_by_param( theColors[theElement->getTag()], 1.0 );
		
		Vector v = theNodes[i]->getCrds();
		
		if (v.Size() == 3 )
		  glVertex3f(v[0], v[1],v[2]);
		else if (v.Size() == 2 )
		  glVertex3f(v[0], v[1],0 );
		
	  }
  glEnd();
}  //2 nodes such as trusses


int DrawWireFrame::getClassTag(){ 
   return DRAW_WIRE; 
}



