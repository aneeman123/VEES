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
** DrawSolid.cpp : draw a finite element as a solid (except for 2Node 
** element, these will somply be lines.
**/

#ifndef _DRAW_SOLID_
#define _DRAW_SOLID_

#include "DrawElement.h"


class DrawSolid:public DrawElement{
 public:
  DrawSolid();

  void drawElement(Element * el, double * colors, 
						   ColorRange range, int colorStyleFlag );

  bool canDraw( Element * el );
  int getClassTag();
 private:
  void drawTwoNode( Element * theElement, double * theColors, int colorFlag );
  void drawFourNode( Element * theElement, double * theColors, int colorFlag );
  void drawEightNode( Element * theElement, double * theColors, int colorFlag );

};

#endif



