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
** DrawWireFrame.h : draw a finite element as a wireframe. Note: the element
** will still be pickable as if it is a solid.
**/

#ifndef _DRAW_WIRE_
#define _DRAW_WIRE_

#include "DrawElement.h"

class DrawWireFrame:public DrawElement{
 public:
  DrawWireFrame();
 void drawElement(Element * el, double * colors, 
						   ColorRange range, int colorStyleFlag );
  bool canDraw( Element * el );
  int getClassTag();

 private:
  void drawTwoNode( Element * theElement, double * theColors, int colorFlag );
  void drawFourNode( Element * theElement, double * theColors,  
					 ColorRange range, int colorFlag );
  void drawEightNode( Element * theElement, double * theColors, 
					  ColorRange range, int colorFlag );

  void drawTwentyNode( Element * theElement, double * theColors, ColorRange range, int colorflag );


};

#endif



