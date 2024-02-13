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
** glutWin.cpp
** Class file for main glut window which holds openGL window.
** Acts as a wrapper for openGL calls in VisWin object.
** It also creates all the little glui objects for interaction
** It is a mix of C and C++ because glut uses pointers to functions
**/

#ifndef _GLUTWIN_H_
#define _GLUT_WIN_H_

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



#include <stdlib.h>
#include <Domain.h>
#include <OPS_Globals.h>
#include <OPS_Stream.h>
#include <StandardStream.h>
#include <GL/gl.h>			/* OpenGL header file */
#include <GL/glu.h>			/* OpenGL utilities header file */
#include <GL/glut.h>
#include <GL/glui.h>

#include "GlutSubWindow.h"

#include "VisWin.h"
#include "DomainViewer.h"
#include "ElementStrings.h"
#include "ColorSchemeNames.h"
#include <ElementViewer.h>
#include <GlyphViewer.h>
#include <ColorBarViewer.h>
#include <auxiliary.h>

//#include <DomainBuilder.h>



class GlutWin {

	GlutWin(); // making constructor private, access through "instance" method
public:

 static GlutWin & Instance();
  
  void createGUI(int argc, char **argv , Domain * dm);
  void createGUI(Domain * dm); 
  void setMainWindow(int i);
  int getMainWindow();
  int getWidth();
  int getHeight();
  void setCurrWidth( int w );
  void setCurrHeight( int h );


  /* VisWin manipulation methods */
  void drawDomain();
  void drawElement();
  void drawGlyph();
  void drawColorBar();

  void reshape( int width, int height );
  void zoomIn( int winNum ); //element viewer, domain viwer,..
  void zoomOut(  int winNum );
  void spinX( int winNum );
  void spinY( int winNum );
  void spinZ( int winNum );
  void resetView( int winNum );
  void recvMouseAction(int button, int state, int x, int y, int winNum);
 

  /* DomainViewer manipulation methods */
  void changeDrawMethod( int method );
  void changeColorMethod( int method );
  void setColorRange( double max, double min, int invert );
  void showHideElement(int index);
  void updateState(); //new
  void updateColorMap();

  /* GlyphViewer draw methods */
  void setGlyphDrawMethod( unsigned char method );
  
  /* element viewer */
  void setElement( int elementNum);

  void makeMovie(int forwards );
  void writePPM();
  void writePPM2();
  void writeSubwindow(int subwin, unsigned char ** array);
  int snapshot();
  
  GlutSubWindow getSubWindow (int i);

  /*** widget panels ***/
  GLUI *glui; // render constrols, main window
  GLUI *glui2; // load timestep panel
  GLUI * glui3; // adjust colormap
  GLUI * glui4; // drawing method for integration point

  /*** widgets ********/
  GLUI_Listbox *drawTypeList;
  GLUI_Listbox * scalarColorList;
  GLUI_Listbox * gaussDrawTypeList;
  GLUI_Spinner * colorMax;
  GLUI_Spinner * colorMin;
  GLUI_Checkbox * invertColorRange;
  GLUI_Rollout * eleChooserRollout;
  GLUI_Checkbox * *eleChoices; // the names for these will come from Domain
  GLUI_EditText * colorRangeNum;
  GLUI_EditText * commitNum;
  GLUI_EditText *selectElement;

  float minColor, maxColor;
  int colorInvert, selectedElement;

  int commitStep;//new


  int names[numEleTypes]; //name list for show/hide element type
  int nameBools[numEleTypes];




  private:

  VisWin * domainWin;
  VisWin * eleWin;
  VisWin * intgrPtWin;
  VisWin * colorBarWin;

  DomainViewer *dv;
  ElementViewer *ev;
  GlyphViewer  *gv;
  ColorBarViewer *cv;

  
  int mainWindow;
  int drawMethod;
  int mainWidth, mainHeight;
  int currWidth,currHeight;

  /**** widget data *******/
  int drawType, colorType, gaussDrawType; 
  int numTypes; //number of finite element types in Domain
  char filename[sizeof(GLUI_String)];//new
  double colorMapRange;
  int logColorScale;
  int endStep; // last commit step
  int saveImage;
  int saveVTK;
  int saveRotation;
  int saveStress;
  int saveEigenTensor;
  int saveEigenLode;
  GlutSubWindow subWins[4];



};

enum widgets { NONE, DRAW_TYPE, COLOR_TYPE, COLOR_MAX, COLOR_MIN, 
			   INVERT_COLOR_RANGE, 
			   SET_FILE_NAME, SET_COMMIT_STEP,LOAD_STEP,
	       INCREMENT_STEP,DECREMENT_STEP,COLOR_MAP, MOVIE, BACKWARDS_MOVIE,
	       SELECT_ELEMENT, GAUSS_DRAW_TYPE, SHOW_HIDE_ELEMENT_TYPE=1000 };

enum viewers { DOMAIN_VIEWER, ELE_VIEWER, INTGR_PT_VIEWER, COLOR_BAR_VIEWER };



/** element show/hide checkboxes start at 1000 */

#define GAP 25
#define DOMAIN_WIN_SIZE 350 //350 x 350
#define ELEMENT_WIN_SIZE 200
#define COLOR_BAR_SIZE 50



#endif
