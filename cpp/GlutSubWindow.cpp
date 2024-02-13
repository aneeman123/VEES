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
** Carries the state for a Glut subwindow. Calculates the new size
** and position on a resize. 
**/

#include "GlutSubWindow.h"

GlutSubWindow::GlutSubWindow() {}

void GlutSubWindow::init(  int superWidth, int superHeight, int placeX, 
			   int placeY, int width, int height, int num )
{
  w =(float) width;
  h = (float) height;
  posX = (float) placeX; // upper left
  posY = (float) placeY; // upper left
  minWidth =(float) superWidth;
  minHeight = (float) superHeight;
  winNum = num;
}
 
int GlutSubWindow::getWinNum() {
  return winNum;
}

/**
 * Calculate new size and position of subwindow on resize.
 * @pre: caller checks that new size is not below minimum
 */
void GlutSubWindow::resize(  int & width, int & height, int superWidth, 
   int superHeight,int & placeX,  int & placeY) 
{
  float increaseW, increaseH;
  
  increaseW = ( (float) superWidth )/minWidth;
  increaseH = ( (float) superHeight )/minHeight;
 
  width = (int)( increaseW*w );
  height = (int)( increaseH*h );

  placeX = (int)( posX*  increaseW );
  placeY = (int) ( posY*increaseH );
}

int  GlutSubWindow::getWidth() {
  return ((int)w);
}


int  GlutSubWindow::getHeight() {
  return ((int)h);
} 


int  GlutSubWindow::getX() {//upper left corner relative to full GUI
  return ((int)posX);
}

int  GlutSubWindow::getY()  {//upper left corner relative to full GUI
  return ((int)posY);
}

/**
 * given a 2D array representing the entire GUI window, 
 * write the pixels for this subwindow
 * @pre: subwindow has the focus
 */
void GlutSubWindow::writeSubwindow( unsigned char ** array) {
  // get the width and height of the subwindow
  // get the offset of the subwindow (x,y)
  // note this is upper left

  int heightOffset, leftOffset, x, y;
  unsigned char* rgb_data; 

  rgb_data =  new unsigned char[3*getWidth()];  
  heightOffset =  getY();
  leftOffset = getX()*3;

  
   for( y = 0; y < getHeight(); y++ ) {// for each row
	glReadPixels(0, getHeight() - 1 - y, getWidth()*3, 
				 1, GL_RGB, GL_UNSIGNED_BYTE, rgb_data);
	
	// calculate local window position
	for(x =0;  x < getWidth()*3; x++) {
	  array[y + heightOffset][x + leftOffset] = rgb_data[x];
	}
  }
   
   delete[] rgb_data;
}
