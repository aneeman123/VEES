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
** HideElement.h : the invisible element. We keep a momento of
** how it was drawn when last seen, so that if someone wants to 
** see it again, it looks the same as before. Note that these must be
** individually created & destroyed because of element-owned state.
**/

#include "HideElement.h"

HideElement::HideElement(){}

/** Nothing to draw for the invisible element */
 void HideElement::drawElement(Element * el, double * colors, 
							   ColorRange range, int colorStyleFlag ){
   return;
}

/** there's no type of element that can't be invisible */  
bool HideElement::canDraw( Element * el ){ 
   return true; 
}

void HideElement::setLastDrawnAs( DrawElement * d) {
  lastDrawnAs = d;
}

/**
 * Reverting to visible form means returning the old draw method
 */ 
DrawElement * HideElement::revert() {
  return lastDrawnAs;
} 
 
int HideElement::getClassTag() {
  return HIDE_ELEMENT;
}
