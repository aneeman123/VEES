/*****************************************************************
**
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
** Filename: VisWin.cpp OpenGL window for 3D viewing
**  GUI API independent.
******************************************************************/

#include "VisWin.h"
#include <stdio.h>
#include <math.h>

/*
 *	Setup OpenGL initialization parameters
 * pre! DomainViewer has a Domain instance inside it. 
 */
VisWin::VisWin ( Viewer *domv, int width, int height )
{ 
  int numElements;
  double  max,winSize, aspect;
  

  scale = 1.0;
  

  init ( width, height );  

  xSpin = ySpin = zSpin = 0;
  cam.setZoom( 0 );


  theViewer = domv;

  max = theViewer->getCenterAndMaxDim( center ); // center the volume over 0,0,0
  

  aspect = ((double)width)/((double)height);

  // works for aspect ratio close to 1
  if (aspect < 1.5 && aspect > 0.5 ) {
    if ( width < height)
      winSize = width;
    else winSize = height;

  }
  
  //else assume a long thin object in a long thin window
  else {
    if (width < height)
      winSize = height;
    else 
      winSize = width;
  }

 
  numElements = theViewer->getNumElements();

  // some reasonable defaults
  if (max > winSize && numElements > 1 ) {
	scale = winSize/(max*1.2);
  }

  else if (max > winSize && numElements == 1 ) {
   scale = winSize/max;
  }

  else if ( max < winSize && numElements > 1)
    scale = (winSize - max)/30;
 
  else if ( max < winSize && numElements == 1)
    //scale =  (winSize - max)/5; works for glyphs
    scale = winSize/max;

  else scale = 10; //default
  

   maxDim = max;
   // printf("width %d height %d maxDim %f winSize %f scale %f\n", width, height, max, winSize, scale);


 
}

 /** recenter the volume based on current state of  viewer.*/
void VisWin::centerVolume() {
  double max,winSize;
  int numElements;
  max = theViewer->getCenterAndMaxDim( center );
  if ( width < height)
	winSize = width;
  else winSize = height;
  
  numElements = theViewer->getNumElements();

  if ( max < winSize && numElements > 1)
    scale = (winSize - max)/30;
  else if ( max < winSize && numElements == 1)
  
   scale =  (winSize - max)/5;
 
  else scale = 1.0;

  maxDim = max;

}


void VisWin::resetView()
{
  xSpin = ySpin = zSpin = 0;
  cam.setZoom( 0 );
}
void VisWin::spinX() {
  xSpin = (xSpin + 4) % 360; 
}

void VisWin::spinY() {
  ySpin = (ySpin + 4) % 360; 
}

void VisWin::spinZ() {
  zSpin = (zSpin + 4) % 360; 
}

/*
 * Setup OpenGL initialization parameters.
 * Precondition: domain handle already set 
 * @param width: width of openGL window
 * @param height: height of openGL window 
 */
void VisWin::init( int wWidth, int wHeight )
{
  
  glMatrixMode(GL_PROJECTION);
  glLoadIdentity();
  gluPerspective(60.0, (((double) wWidth)/ ((double) wHeight)), 0.1, 1000);
  glEnable(GL_DEPTH_TEST); //allow us to reset z buffer as we draw things

  glShadeModel(GL_SMOOTH);
  glPolygonMode(GL_FRONT_AND_BACK, GL_FILL);

  glEnable(GL_POLYGON_SMOOTH);
  glEnable(GL_LINE_SMOOTH);
  
  glHint( GL_POLYGON_SMOOTH_HINT,  GL_NICEST );
  glHint( GL_LINE_SMOOTH_HINT, GL_NICEST );
  
  glEnable(GL_DEPTH_TEST);
  glEnable(GL_NORMALIZE);

  initLight();
  width = wWidth;
  height = wHeight;
  //  cam.initCamera( 100 );

}

/**
 * Sets view as orthographic projection, with box based on 
 * dimensions & scale of object to be drawn. Affects only
 * this VisWin.
 * @pre! Domain has been loaded and camera initialized */
void VisWin::setOrtho()
{

  // box with aspect ratio of 1.0 so as not to
  // morph object

  GLdouble left    = -maxDim*scale*.4; //x
  GLdouble right   =  maxDim*scale*.4;
  
  GLdouble bottom  = -maxDim*scale*.4; //y
  GLdouble top     =  maxDim*scale*.4; 
 
  GLdouble nearDim   = -2*maxDim*scale; //z
  GLdouble farDim    =  2*maxDim*scale;

  glMatrixMode(GL_PROJECTION); 
  glLoadIdentity(); 
  glOrtho(left, right, bottom, top, nearDim, farDim); 
  
  glMatrixMode(GL_MODELVIEW); 
  glLoadIdentity();  
}

/* Pick up mouse pick events on GL_window.
 * Note Y position is reversed. L->R = 0x to MAXx
 * Bottom to top = maxY to 0Y
 * negative values are erroneous. keyboard events case insensitive
 * @return the selected drawn object number, or -1 for no pick
 */
int VisWin::handle(int button, int state, int x, int y)
{
  switch( state ) {
  case MOUSE_RELEASE:
    if ( button == MOUSE_LEFT_BUTTON ) {
	 
      // Now try to pick at the current mouse position
	  return  (pick( x, y ));//y reversal done in viswin
    }
	break;

  } //switch
  return -1;
}


/*
 * InitLight give us ambient , diffuse and a 
 * position light on the scene .
 * so we can see color. these are all GL_LIGHT0
 * */
void VisWin::initLight(){

  GLfloat white_light[] = {1.0, 1.0, 1.0,1.0};
  GLfloat lmodel_ambient[] = {0.4,0.4,0.4,1.0};
 
  glClearColor(0.0, 0.0, 0.0, 0.0);
  glShadeModel(GL_SMOOTH);

  // glLightfv(GL_LIGHT0,GL_AMBIENT,white_light);
  glLightModelfv(GL_LIGHT_MODEL_AMBIENT,lmodel_ambient);
 
  //general characteristics
  // glLightModeli(GL_LIGHT_MODEL_TWO_SIDE,GL_TRUE);
  glEnable(GL_LIGHTING);
  //  glEnable(GL_LIGHT0);
  glEnable(GL_DEPTH_TEST);


 }



/** Handle a pick: decide which numbered geometric 
 * shape was selected.
 * @returns: selection closest to camera eye (front) or
 * -1 if no hits
*/
int VisWin::pick(int x, int y )
{
  
  GLint		vp[4];			// viewport parameter
  GLint		hits;			// number of hits
  GLuint  hitbuf[512] = {0};
  int i, j;
  
  //  get viewport info
  glGetIntegerv (GL_VIEWPORT, vp);
  
  //  setup for pick/name buffer
 
  glSelectBuffer(HIT_BUFSIZE, hitbuf);
  glInitNames();
 // glPushName(0);
  glRenderMode(GL_SELECT);
  
  //  about to re-render everything using 5x5 pick-window,
  //  better save current transform first
 
  glMatrixMode(GL_PROJECTION);
  glPushMatrix();
  glLoadIdentity();
  
  // Pick just in a 3*3 window around the mouse position
  // vp[3]-y because the window coordinates and opengl coordinates
  // are reversed in y direction
  
  gluPickMatrix( (GLdouble) x, (GLdouble) (vp[3]-y), 3.0, 3.0, vp );
  gluPerspective(60.0, (double) width/ (double) height, 0.1, 1000);
  
  draw();
 
  glMatrixMode(GL_PROJECTION);
  glPopMatrix();
  
  //  get back into render mode
  hits = glRenderMode(GL_RENDER);
  
  // PROF. WILHELMS' CODE 
 unsigned int *ptr,numnames,zmin,zmax,lastname;
 unsigned int tri = 0;
 int done = 0;
 ptr = hitbuf;

 for (i=0; i < hits; i++)
   { if (done) break;
     numnames = *ptr;
     ptr++;
     zmin = *ptr;
     ptr++;
     zmax = *ptr;
     ptr++;
     // take the first name from the first hit that's
     // not empty 
     for (j = 0; j < numnames; j++)
       {
		 lastname = *ptr;
         tri = *ptr;
         done = 1;
		 printf("picked element %d\n", tri);
		 ptr++;
       }
   }

 if (hits == 0 ) 
   return -1;
 return tri;
}



/*
 * This method is called whenever there is a need to refresh the 
 * OpenGL window. 
 * Precondition: camera position set
 * Postcondition: Buffers must be swapped afterward by caller 
 * (platform/GUI dependent)
 */
void VisWin::draw()
{
  glEnable(GL_NORMALIZE);
  glDisable(GL_LIGHT0);
  glDisable(GL_LIGHT1);
  glDisable(GL_LIGHT2);
  glDisable(GL_LIGHT3);
  glDisable(GL_LIGHT4);
  glDisable(GL_LIGHT5);
  glDisable(GL_LIGHT6);

  glMatrixMode(GL_MODELVIEW);
  glClearColor(1.0, 1.0, 1.0, 1.0); //white background
  glClear( GL_COLOR_BUFFER_BIT | GL_DEPTH_BUFFER_BIT );
  
  initLight(); /* ambient */ 

  glLoadIdentity ();          //clear the matrix  

  lightMinusX();
  cam.calculateCameraPosition();
  cam.lookAt();

  drawStaticScene();

//	glutSwapBuffers(); caller should swap - GUI dependent
}


/* respond to window resize 
 * Postcondition: Buffers must be swapped afterward by caller \
 * (platform/GUI dependent)
 */
void VisWin::reshape( int width, int height ) {

  glViewport( 0, 0, width, height ); //set viewport to window new size 
  glMatrixMode( GL_PROJECTION );
  glLoadIdentity();
  gluPerspective( 60.0, (double) width/ (double) height, 0.1, 1000 );

  glMatrixMode (GL_MODELVIEW); //need modelview for gluLookAt call
  cam.lookAt();

}

void VisWin::lightMinusX() {
  GLfloat light_pos[4];// ={ -2.0, 2.0, 2.0, 1.0 };
  GLfloat size;

  GLfloat light_diff[4] = {1.0,1.0,1.0,1.0 };
 
 
  if ( maxDim > 0.0 && scale > 1.0 ) {	
	size = maxDim*1.2;
  }
  
  else
	size = 3.0;
  
  light_pos[0] = size*2;
  light_pos[1] = size*2;//size;
  light_pos[2] = size*2;//size;
  light_pos[3] = 1.0;  // carries diffuse, specular and attenuation
  
   glLightfv( GL_LIGHT0, GL_POSITION, light_pos );
   glEnable( GL_LIGHT0 );
}

void VisWin::lightMinusZ() {
  GLfloat light_pos[4];// ={ 1.90, 1.70, -2.0, 1.0 }; 
  GLfloat light_diff[] = {1.0,1.0,1.0,1.0 };

  light_pos[0] = 1.9*scale;
  light_pos[1] = 1.7*scale;
  light_pos[2] = -2.0*scale;
  light_pos[3] = 1.0;

  glLightfv( GL_LIGHT1, GL_POSITION, light_pos );
  //glLightfv( GL_LIGHT1, GL_SPECULAR, light_diff );
  glEnable( GL_LIGHT1 );
}

void VisWin::lightMinusY() {
  GLfloat light_pos[4];// ={ 2.0, -2.0, 2.0, 1.0 }; 
  GLfloat light_diff[] = {1.0,1.0,1.0,1.0 };

  light_pos[0] = 2.0*scale;
  light_pos[1] = -2.0*scale;
  light_pos[2] = 2.0*scale;
  light_pos[3] = 1.0;

  glLightfv( GL_LIGHT2, GL_POSITION, light_pos );
  glEnable( GL_LIGHT2 );
}

void VisWin::lightAbove() {
  GLfloat light_pos[] ={ 0.0, 4.0, 0.0, 1.0 }; 
  GLfloat light_diff[] = {1.0,1.0,1.0,1.0 };
  glLightfv( GL_LIGHT3, GL_POSITION, light_pos );
  glEnable( GL_LIGHT3 );
}

/**
 * Let the user see the time/increment step in the openGL window. 
 * hardcoded location for small window.
*/
void VisWin::drawIncrNumber(  char * format, int commitStep ) {

  //just like draw commands
  glEnable(GL_NORMALIZE);
  glMatrixMode(GL_MODELVIEW);
  initLight(); // ambient  
  glLoadIdentity ();          //clear the matrix  
  cam.calculateCameraPosition();
  cam.lookAt();  

 
  glScalef(scale,scale,scale); // like drawStaticScene without rotate
  glTranslatef(-center[X],-center[Y],-center[Z]);
  lightMinusX(); //good 

  gluOrtho2D(-100, 100, -100, 100); //flatten to draw text
 
  // move the raster position and write text
  set_saturation_by_param(0.1, 1.0);
  setfont("helvetica", 16);
  drawstr(-30, -113, format, commitStep ); //y direction is reverse
}

void VisWin::drawTwoNumbers( char * format, int num1, int num2 ) {

  //just like draw commands
  glEnable(GL_NORMALIZE);
  glMatrixMode(GL_MODELVIEW);
  initLight(); // ambient  
  glLoadIdentity ();          //clear the matrix  
  cam.calculateCameraPosition();
  cam.lookAt();  

 
  glScalef(scale,scale,scale); // like drawStaticScene without rotate
  glTranslatef(-center[X],-center[Y],-center[Z]);
  lightMinusX(); //good 

  gluOrtho2D(-100, 100, -100, 100); //flatten to draw text
 
  // move the raster position and write text
  set_saturation_by_param(0.1, 1.0);
  setfont("helvetica", 16);
  drawstr(-90, -113, format, num1, num2 ); //y direction is reverse
}




/**
 * rotate the volume and then draw the objects
 */
void VisWin::drawStaticScene() {
  // printf("starting VisWin::drawStaticScene\n");
  //lightAbove(); nada
  glScalef(scale,scale,scale);

  rotateVolume();
  glTranslatef(-center[X],-center[Y],-center[Z]);
  lightMinusX(); //good 
 
  theViewer->draw();
}

/* rotate volume (not camera) */
void VisWin::rotateVolume() {

 glRotatef(xSpin, 0.0, 1.0, 0.0);
 glRotatef(ySpin, 1.0, 0.0, 0.0);
 glRotatef(zSpin, 0.0, 0.0, 1.0);
  
}

/* set how much to move camera in from initial position */
void VisWin::setZoom(double zoom) {
  cam.setZoom( zoom );
}


void VisWin::zoomIn() {
  cam.zoomIn();
}

void VisWin::zoomOut() {
  cam.zoomOut();
}
