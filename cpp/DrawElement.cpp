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
** DrawElement.cpp is an abstract class for a drawing method.
** The idea is to have lots of ways to draw an element.
**/
#include "DrawElement.h"

DrawElement::DrawElement() {}
DrawElement::~DrawElement() {}


void DrawElement::drawElement(Element * el, double * colors, 
							  ColorRange range, int colorStyleFlag )
{
  printf("Error: drawElement must be implemented in subclass\n");
}

/** return if this method is capable of drawing 
 * this element type */
bool  DrawElement::canDraw( Element * el )
{ 
  return false; 
}

/**
 * Find out what DrawMethod this is */
int DrawElement::getClassTag(){
  printf("Error: getClassTag must be implemented in subclass\n");
  return -1;
}
