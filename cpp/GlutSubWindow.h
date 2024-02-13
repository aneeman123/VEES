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

#ifdef WIN32
#include <windows.h>
#endif

#if defined(__APPLE__)&& defined(__MACH__)
#include <GLUT/glut.h>
#include <OpenGL/gl.h>
#include <OpenGL/glu.h>
#else
#include <GL/gl.h>			/* OpenGL header file */
#include <GL/glu.h>			/* OpenGL utilities header file */
#include <GL/glut.h>	
#endif

class GlutSubWindow{
 public:
  GlutSubWindow();

  void init(  int superWidth, int superHeight, int placeX, 
	      int placeY, int width, int height, int num );

  void resize(  int & width, int & height, int superWidth, int superHeight,
		int & placeX, int & placeY );
  int getWinNum();
  
  int getWidth();
  int getHeight();
  int getX(); //upper left corner relative to full GUI
  int getY();
  void writeSubwindow( unsigned char ** array);

 private:
  float w;
  float h;

  float posX;
  float posY;

  float minWidth; //base container window size for scaling factor
  float minHeight; 

  int winNum;

};
