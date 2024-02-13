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

#ifndef _HIDE_ELEMENT_
#define _HIDE_ELEMENT_

#include "DrawElement.h"


class HideElement:public DrawElement{
 public:
  HideElement();

 void drawElement(Element * el, double * colors, 
				  ColorRange range, int colorStyleFlag );
 int getClassTag();

 bool canDraw( Element * el );

 void setLastDrawnAs( DrawElement * d);
 
 DrawElement * revert(); 
 

 private:
 DrawElement * lastDrawnAs;
};

#endif



