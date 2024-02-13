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

#include "GlutWin.h"
#include "auxiliary.h"
/**********************************************************
 * GlutWin class methods. Basically wrap the VisWin methods 
 * and call GUI API dependent swap buffer.
 *********************************************************/

/** make handle to global so free floating functions can use */
GlutWin::GlutWin() {
  // viswin = new VisWin( & dv );
  maxColor = 1.0;
  minColor = 0.0;
  commitStep = 0;
  colorMapRange = 0.0;
  logColorScale = 0;
  selectedElement = 0;

  //time step iterator
  strcpy(filename, "domain");
  commitStep = endStep = 1;
  saveImage = saveVTK = saveRotation = saveStress =  saveEigenLode = 0;
  saveEigenTensor = 0;
}

/** save the handle to the main window (little gluis need it) */
void GlutWin::setMainWindow( int winNum ) {
  mainWindow = winNum;
}

/** return the handle to the main (container) opengl window */
int GlutWin::getMainWindow() {
  return mainWindow;
}

/** ask wrapped openGl window to draw itself
 * large left side window showing entire domain
 */
void GlutWin::drawDomain() {
  domainWin->draw();
  glutSwapBuffers(); 
}

/** redraw upper right single element window */ 
void GlutWin::drawElement() {
 char com[] = "ELEMENT %d";
 eleWin->draw();

 // label with lement number
 eleWin->drawIncrNumber( com, ev->getElement() );
 glutSwapBuffers();
}

/** redraw lower right single glyph window */
void GlutWin::drawGlyph() {
  char com[] = "STEP %d Integr. Point %d";
  intgrPtWin->setOrtho();
  intgrPtWin->draw();

  // label with commit step and integration pt num
  intgrPtWin->drawTwoNumbers( com, commitStep,
			      gv->getIntgrPoint() );
  glutSwapBuffers();
}

/** redraw lower LEFT color bar window */
void GlutWin::drawColorBar() {
  colorBarWin->draw();
  glutSwapBuffers();
}


/** ask wrapped openGL subwindows to respond to resizing */
void GlutWin::reshape( int width, int height ) {

  domainWin->reshape( width, height );
  eleWin->reshape( width, height );
  intgrPtWin->reshape( width, height );
  colorBarWin->reshape( width, height );

}

GlutWin & GlutWin::Instance()
{
  static GlutWin glutwin; // a global, but 'static' keeps it private
  return glutwin;
}

/** ask subwindow to rotate about x axis.
 * @param winNum: which subwindow to spin 
 */
void  GlutWin::spinX( int winNum ) {
  if ( winNum == DOMAIN_VIEWER )
	domainWin->spinX();
  else if ( winNum == ELE_VIEWER )
	eleWin->spinX();
  else if ( winNum == INTGR_PT_VIEWER )
	intgrPtWin->spinX();
}

/** ask subwindow to rotate about y axis.
 * @param winNum: which subwindow to spin 
 */
void  GlutWin::spinY( int winNum ) {
  if ( winNum == DOMAIN_VIEWER )
	domainWin->spinY();
  else if ( winNum == ELE_VIEWER )
	eleWin->spinY();
  else if ( winNum == INTGR_PT_VIEWER )
	intgrPtWin->spinY();
}

/** ask subwindow to rotate about z axis.
 * @param winNum: which subwindow to spin 
 */
void  GlutWin::spinZ( int winNum ) {
  if ( winNum == DOMAIN_VIEWER )
	domainWin->spinZ();
  else if ( winNum == ELE_VIEWER )
	eleWin->spinZ();
  else if ( winNum == INTGR_PT_VIEWER )
	intgrPtWin->spinZ();
}

/**
 * ask subwindow to undo all zooms and rotations
 * @param winNum: which subwindow to reset
 */
void GlutWin::resetView( int winNum ) {
  if ( winNum == DOMAIN_VIEWER )
	domainWin->resetView();
  else if ( winNum == ELE_VIEWER )
	eleWin->resetView();
  else if ( winNum == INTGR_PT_VIEWER )
	intgrPtWin->resetView();
}

void GlutWin::setGlyphDrawMethod( unsigned char method ) {
  gv->setDrawMethod(method);
}

void GlutWin::setElement( int elementNum ) {
  ev->setElement( elementNum );
  eleWin->centerVolume();
  glutSetWindow( subWins[ELE_VIEWER].getWinNum() );//switch to elewin
  drawElement();
}

void GlutWin::recvMouseAction(int button, int state, int x, int y, int winNum)
{
  int i; 
  // handle main window pick by setting element in
  // element viewer window
  if ( winNum == DOMAIN_VIEWER ) { 
	i = domainWin->handle( button, state, x, y);

	if (i != -1 ) {
	  ev->setElement( i );
	  eleWin->centerVolume();
	  glutSetWindow( subWins[ELE_VIEWER].getWinNum() );//switch to elewin
	  drawElement();
	}
  }
  else if ( winNum == ELE_VIEWER ) {
	i = eleWin->handle( button, state, x, y );
	ev->setSelected( i );
	if (i != -1 ) {
	  gv->setElement( ev->getElement() );
	  gv->setIntegrPoint( i );

	  glutSetWindow( subWins[ELE_VIEWER].getWinNum() );//switch to elewin
	  drawElement();

	  glutSetWindow(subWins[INTGR_PT_VIEWER].getWinNum() );//switch to elem
	  drawGlyph();
	}
  }
  else if (winNum == INTGR_PT_VIEWER ) {
	; // NOT SURE WHAT HAPPENS HERE
  }

}

GlutSubWindow GlutWin::getSubWindow (int i) {
  if ( i  > 3 || i < 0 ) {
	printf("Unable to retrieve subwindow %d\n", i);
	return subWins[0];
  }
  return subWins[i];
}

void GlutWin::changeDrawMethod( int method )
{
  dv->setDrawMethodForAll( method );

}

void GlutWin::changeColorMethod( int method )
{
  if ( dv->setScalarColorMethodForAll( method ) ) {

  //update GUI to reflect range
    const ColorRange & range = dv->getColorRange( method );
    colorRangeNum->set_float_val(range.getBaseValue());
    
    if( range.isLogScaled() )
      logColorScale = 1;
    else
      logColorScale = 0;
    
    glui3->sync_live();
    
    dv->draw();
    glutSwapBuffers();
  }
}

void GlutWin::setColorRange( double min, double max, int invert ) { 
  dv->setColorRange(min, max, invert );
  dv->draw();
  glutSwapBuffers(); 
} 

void GlutWin::showHideElement(int index) {
  dv->toggleShowElement( names[index],
                         nameBools[index]);

  fflush(stdout);
}

/**
 * Load a time step */
void GlutWin::updateState(){
  printf("\n\n GlutWin::update State!!!!!!!! step %d\n\n", commitStep);
  fflush(stdout);

  if( dv->updateState(filename ,commitStep) != 0 ) {
    printf("updateState failed\n");
    fflush(stdout);
  }
    
  //ev->updateState();
   gv->updateState(commitStep);

   glutSetWindow(subWins[INTGR_PT_VIEWER].getWinNum() );//switch to elem
   drawGlyph();
   
}

void GlutWin::updateColorMap() {

  printf("GlutWin::updateColorMap() %f\n", colorRangeNum->get_float_val());
  if ( logColorScale ) 
    dv->updateColorMap( true,colorRangeNum->get_float_val() );
  else
    dv->updateColorMap( false, colorRangeNum->get_float_val() );
}

/**
 * Load timesteps 1 thru commit step */
void GlutWin::makeMovie(int widget) {
  int i;
  // This should really happen in a differnt class, but lets try 
  int stop;
  if (widget == MOVIE ){ // from commit to up to end
    stop = endStep;

    for ( i = commitStep; i <= endStep; i++ ) {
      printf("\n\n GlutWin::makeMovie!!!!!!!! step %d\n\n", i);
      fflush(stdout);
   
      dv->updateState(filename , i);
      glutSetWindow(subWins[ DOMAIN_VIEWER].getWinNum() );
      drawDomain();
      
      //ev->updateState();
      
      gv->updateState(i);
      glutSetWindow(subWins[INTGR_PT_VIEWER].getWinNum() );//switch to elem
      drawGlyph();

      commitStep = i;
      glui2->sync_live();
      if (saveImage) {
	writePPM2();
      }
      if(saveVTK) {
	dv->saveEigenTensor(i);
	dv->saveRegularVolumeStiffness(i);
      }
      if (saveRotation) {
	dv->saveRotationAccuracy(i);
      }
      if (saveStress) {
      	dv->saveStressTensor(i);
	dv->saveStressLodeAngle( i);
      }
      if(saveEigenTensor) {
	dv->saveEigenTensor(i);
      }
     if(saveEigenLode) {
	dv->saveEigenLodeAngle(i);
      }

    }

  }
  else {
    stop = 1;

    for( i = commitStep; i >= stop; i--) {// from commit to 1
      printf("\n\n GlutWin::makeMovie, no file save!! step %d\n\n", i);
      fflush(stdout);
   
      dv->updateState(filename , i);
      glutSetWindow(subWins[ DOMAIN_VIEWER].getWinNum() );
      drawDomain();
      
      //ev->updateState();
      
      gv->updateState(i);
      glutSetWindow(subWins[INTGR_PT_VIEWER].getWinNum() );//switch to elem
      drawGlyph();

      commitStep = i;
      glui2->sync_live();

      if (saveImage) {
	writePPM2();

      }
      if (saveVTK) {
	dv->saveRegularVolumeStiffness(i);
      }
      if (saveRotation) {
	dv->saveRotationAccuracy(i);
      }
      if (saveStress) {
      	dv->saveStressTensor(i);
	dv->saveStressLodeAngle( i);
      }
      if(saveEigenTensor) {
	dv->saveEigenTensor(i);
      }
    if(saveEigenLode) {
	dv->saveEigenLodeAngle(i);
      }
    }
  } 

}


/************************************************************
 * Free functions since function pointers must be passed
 * to GLUT. (It cannot be done with class member functions) 
 *************************************************************/
static void domainDisplay(void)
{
  	GlutWin::Instance().drawDomain();
}

/* resize main window  */
static void domainReshape(int width, int height)
{
  //recalculate width & hieight since this 
  // is a subwindow 
  GlutWin::Instance().reshape( width, height ); 

}

static void elementDisplay(void)
{
  GlutWin::Instance().drawElement();
}

static void glyphDisplay(void)
{
	GlutWin::Instance().drawGlyph();
}

static void colorDisplay(void)
{
	GlutWin::Instance().drawColorBar();
}


static void mainMouseMotion(int x, int y)
{

  
}



/* button: GLUT_LEFT_BUTTON (0), GLUT_MIDDLE_BUTTON(1), or GLUT_RIGHT_BUTTON(2)
 * state: GLUT_UP (1=release) or  GLUT_DOWN (0=press)
 * x,y: window relative coords when mouse button state changed
 */
static void mainMouse(int button, int state, int x, int y)
{
  int mouseButton,mouseState;

  if (button == GLUT_LEFT_BUTTON)
	mouseButton = MOUSE_LEFT_BUTTON; 
  else if (button == GLUT_RIGHT_BUTTON )
	mouseButton = MOUSE_RIGHT_BUTTON;
  else 
	mouseButton = MOUSE_MIDDLE_BUTTON;

  if ( state == GLUT_DOWN )
	mouseState = MOUSE_PRESS;
  else 
	mouseState = MOUSE_RELEASE;
  
  	GlutWin::Instance().recvMouseAction( 
					   mouseButton, mouseState, x, y, DOMAIN_VIEWER);

}

static void elementMouse(int button, int state, int x, int y)
{
  int mouseButton,mouseState;

  if (button == GLUT_LEFT_BUTTON)
	mouseButton = MOUSE_LEFT_BUTTON; 
  else if (button == GLUT_RIGHT_BUTTON )
	mouseButton = MOUSE_RIGHT_BUTTON;
  else 
	mouseButton = MOUSE_MIDDLE_BUTTON;

  if ( state == GLUT_DOWN )
	mouseState = MOUSE_PRESS;
  else 
	mouseState = MOUSE_RELEASE;
  
  	GlutWin::Instance().recvMouseAction( 
					   mouseButton, mouseState, x, y, ELE_VIEWER);

}

/**
 * handle all the different kinds of callbacks from the glui
 */
static void gluiCallback( int widget ) {
 
  // +++++++++ Named Widget +++++++++++
  switch( widget ) {
    case( DRAW_TYPE ):
	  GlutWin::Instance().changeDrawMethod(  
 	       	 GlutWin::Instance().drawTypeList->get_int_val());
	break;

    case( COLOR_TYPE ):
	  GlutWin::Instance().changeColorMethod( 
			  GlutWin::Instance().scalarColorList->get_int_val());
	break;
    case( COLOR_MAX ):
	// don't let min and max cross each other
	  if(  GlutWin::Instance().maxColor <  GlutWin::Instance().minColor )
	    GlutWin::Instance().colorMax->set_float_val(
									   GlutWin::Instance().minColor);
	    GlutWin::Instance().setColorRange( GlutWin::Instance().minColor, 
                       GlutWin::Instance().maxColor,
                  GlutWin::Instance().colorInvert );
	  break;

    case(COLOR_MIN):
	  if(  GlutWin::Instance().minColor > GlutWin::Instance().maxColor )
	    GlutWin::Instance().colorMin->set_float_val(
									   GlutWin::Instance().maxColor);
 	    GlutWin::Instance().setColorRange( GlutWin::Instance().minColor, 
                        GlutWin::Instance().maxColor,
                  GlutWin::Instance().colorInvert );

	  break;
    case( INVERT_COLOR_RANGE ):
	  GlutWin::Instance().setColorRange( GlutWin::Instance().minColor, 
                       GlutWin::Instance().maxColor,
                  GlutWin::Instance().colorInvert );
	 break;
 
    case ( INCREMENT_STEP ): {
	  
	  GlutWin::Instance().commitStep++;
	  GlutWin::Instance().updateState();  // now load the step
	  GlutWin::Instance().glui2->sync_live();
	  if ( GlutWin::Instance().snapshot())
	     GlutWin::Instance().writePPM2();
	   break;
     // now load the step
     }

  case ( DECREMENT_STEP ): {
	if (GlutWin::Instance().commitStep > 0 ) {
	  GlutWin::Instance().commitStep--;
	  GlutWin::Instance().updateState();  // now load the step
   
	  GlutWin::Instance().glui2->sync_live();
	  if ( GlutWin::Instance().snapshot())
	     GlutWin::Instance().writePPM2();
	}
	  break;
    
     }
    case( LOAD_STEP ) : {
	   GlutWin::Instance().updateState();
	   GlutWin::Instance().glui2->sync_live(); 
	   if ( GlutWin::Instance().snapshot())
	     GlutWin::Instance().writePPM2();
	   break;
    } 
   case( COLOR_MAP ) : {
        GlutWin::Instance().updateColorMap();
         break; 
   }

  case ( MOVIE ) : {
    // start at 1, run update until selected commit step
    GlutWin::Instance().makeMovie( MOVIE );
    break;
  }

    case( BACKWARDS_MOVIE )  : {
    // start at curr, run down to 1
    GlutWin::Instance().makeMovie( BACKWARDS_MOVIE );
    break;
  }
  case ( SELECT_ELEMENT ) : {
    GlutWin::Instance().setElement( 
		    GlutWin::Instance().selectElement->get_int_val()  );

    break;
  }
   
  case (GAUSS_DRAW_TYPE) : {
    printf("GAUSS_DRAW_TYPE");
    break;
  }

 }
  

  if( widget >= SHOW_HIDE_ELEMENT_TYPE ) {
	GlutWin::Instance().showHideElement(widget-SHOW_HIDE_ELEMENT_TYPE);
	
  }
  
}


/* keyboard for element viewer */
static void eleKeyboard(unsigned char key, int x, int y)
{
  switch (key) {
      case 27:	exit(0);
	   
      case('+'):
		GlutWin::Instance().zoomIn( ELE_VIEWER );
		break;
      case('-'):
		GlutWin::Instance().zoomOut( ELE_VIEWER );
		break;
      case 'y':
		GlutWin::Instance().spinY( ELE_VIEWER ); 
 	    break;
      case 'x':
		GlutWin::Instance().spinX( ELE_VIEWER ); 
		break;
      case 'z':
		GlutWin::Instance().spinZ( ELE_VIEWER ); 
		break;
     case 'r':
	   GlutWin::Instance().resetView( ELE_VIEWER ); 
	   break;
   }
   GlutWin::Instance().drawElement();
}

/* keyboard for domain viewer */
static void domKeyboard(unsigned char key, int x, int y)
{
  switch (key) {
      case 27:	exit(0);
	   
      case('+'):
		GlutWin::Instance().zoomIn( DOMAIN_VIEWER );
		break;
      case('-'):
		GlutWin::Instance().zoomOut( DOMAIN_VIEWER );
		break;
      case 'y':
		GlutWin::Instance().spinY( DOMAIN_VIEWER ); 
 	    break;
      case 'x':
		GlutWin::Instance().spinX( DOMAIN_VIEWER ); 
		break;
      case 'z':
		GlutWin::Instance().spinZ( DOMAIN_VIEWER ); 
		break;
     case 'r':
	   GlutWin::Instance().resetView( DOMAIN_VIEWER ); 
	   break;
   }
   GlutWin::Instance().drawDomain();
}

/* keyboard for domain viewer */
static void intgrKeyboard(unsigned char key, int x, int y)
{
  switch (key) {
      case 27:	exit(0);
	   
      case('+'):
		GlutWin::Instance().zoomIn( INTGR_PT_VIEWER );
		break;
      case('-'):
		GlutWin::Instance().zoomOut( INTGR_PT_VIEWER );
		break;
      case 'y':
		GlutWin::Instance().spinY( INTGR_PT_VIEWER ); 
 	    break;
      case 'x':
		GlutWin::Instance().spinX( INTGR_PT_VIEWER ); 
		break;
      case 'z':
		GlutWin::Instance().spinZ( INTGR_PT_VIEWER ); 
		break;
     case 'r':
	   GlutWin::Instance().resetView( INTGR_PT_VIEWER ); 
	   break;
   
     case '0':
     case '1':
     case '2':
      case'3':
     case '4':
     case '5':
     case '6':
     case '7':
     case '8':
  case '9':
  case 'S':
  case 'G':
	  GlutWin::Instance().setGlyphDrawMethod( key );

	  break;
   }
   GlutWin::Instance().drawGlyph();
}

/* no string arguments */
void GlutWin::createGUI(Domain * dm) { 
  createGUI( 0, NULL, dm );
}

void GlutWin::zoomIn(  int winNum ) {
  if ( winNum == DOMAIN_VIEWER )
	domainWin->zoomIn();
  else if ( winNum == ELE_VIEWER )
	eleWin->zoomIn();
  else if ( winNum == INTGR_PT_VIEWER )
	intgrPtWin->zoomIn();
}

void GlutWin::zoomOut(  int winNum ) {
  if ( winNum == DOMAIN_VIEWER )
	domainWin->zoomOut();
  else if ( winNum == ELE_VIEWER )
	eleWin->zoomOut();
  else if ( winNum == INTGR_PT_VIEWER )
	intgrPtWin->zoomOut();
}

/** The glui idle callback must ask the main window to post 
 * a redisplay 
 */
void myGlutIdle() {
  glutSetWindow(GlutWin::Instance().getMainWindow());
  glutPostRedisplay();
}



/** create window holding subwindows */
void mainDisplay(void)
{
    glClearColor(0.8, 0.8, 0.8, 0.0);
    glClear(GL_COLOR_BUFFER_BIT | GL_DEPTH_BUFFER_BIT);
    
	/* black text */
    glColor3ub(0, 0, 0);

    /* label above subwindow inside gap contains menu */
    setfont("helvetica", 12);
    // x, y, format
    drawstr(GAP,396, "MENU: Zoom In [+] Out [-]   Rotate [x] [y] [z]  Reset [r]");
 
	drawstr( GAP+DOMAIN_WIN_SIZE+GAP, GAP+ELEMENT_WIN_SIZE+(GAP/2),"  Element  " );

	drawstr( GAP+DOMAIN_WIN_SIZE+GAP, 
			 GAP+ELEMENT_WIN_SIZE+ GAP + ELEMENT_WIN_SIZE+(GAP/2),
			 "  Integration Point  " );
    glutSwapBuffers();
}

/**
 * when window gets resized, reset subwindow sizes, projection,etc.
 */
void mainReshape(int width,  int height) 
{
  int w, h, x, y;
  GlutSubWindow subWin;

  if ( width < GlutWin::Instance().getWidth() ||
       height < GlutWin::Instance().getHeight() )
    {
      width = GlutWin::Instance().getWidth();
      height = GlutWin::Instance().getHeight();
    }

  GlutWin::Instance().setCurrWidth( width );
  GlutWin::Instance().setCurrHeight( height );

  glViewport( 0, 0, width, height );
  glMatrixMode( GL_PROJECTION );
  glLoadIdentity();
  gluOrtho2D( 0, width, height, 0 );
  glMatrixMode( GL_MODELVIEW );
  glLoadIdentity();
  
  subWin = GlutWin::Instance().getSubWindow( DOMAIN_VIEWER);
  glutSetWindow( subWin.getWinNum() );
  subWin.resize( w, h, width, height, x, y );
  glutReshapeWindow( w, h );
  glutPositionWindow( x, y);
  glutSetWindow( subWin.getWinNum() );

  subWin = GlutWin::Instance().getSubWindow( ELE_VIEWER);
  glutSetWindow( subWin.getWinNum() );
  subWin.resize( w, h, width, height, x, y );
  glutReshapeWindow( w, h );
  glutPositionWindow( x, y);
  glutSetWindow( subWin.getWinNum() );

  subWin = GlutWin::Instance().getSubWindow(INTGR_PT_VIEWER);
  glutSetWindow( subWin.getWinNum() );
  subWin.resize( w, h, width, height, x, y );
  glutReshapeWindow( w, h );
  glutPositionWindow( x, y);
  glutSetWindow( subWin.getWinNum() );


  subWin = GlutWin::Instance().getSubWindow(COLOR_BAR_VIEWER);
  glutSetWindow( subWin.getWinNum() );
  subWin.resize( w, h, width, height, x, y );
  glutReshapeWindow( w, h );
  glutPositionWindow( x, y);
  glutSetWindow( subWin.getWinNum() );

}

  /** get base dimension of container window */
int GlutWin::getWidth(){
    return mainWidth;
}

  /** get base dimension of container window */
int GlutWin::getHeight(){
    return mainHeight;
}

/** set current dimension of container window */
void GlutWin::setCurrWidth( int w ){
    currWidth = w;
}

  /** set current dimension of container window */
void GlutWin::setCurrHeight( int h ){
    currHeight = h;
}

/**
 * create main window, subwindows, and some little glui widgets 
 * and add callbacks for widgets */
void GlutWin::createGUI(int argc, char **argv, Domain * dm ) {
  int i;
  currWidth = mainWidth =   DOMAIN_WIN_SIZE + ELEMENT_WIN_SIZE + 3*GAP;
  currHeight = mainHeight =  DOMAIN_WIN_SIZE + COLOR_BAR_SIZE + 3*GAP + 55;
 
  dv = new DomainViewer();
  ev = new ElementViewer();
  gv = new GlyphViewer();
  cv = new ColorBarViewer();

  dv->setDomain( dm );

  ev->setDomain( dm );
 
  gv->setDomain( dm );

  domainWin   = new VisWin( dv, DOMAIN_WIN_SIZE , DOMAIN_WIN_SIZE );
  eleWin      = new VisWin( ev, ELEMENT_WIN_SIZE, ELEMENT_WIN_SIZE );
  intgrPtWin  = new VisWin( gv,  ELEMENT_WIN_SIZE, ELEMENT_WIN_SIZE );
  colorBarWin = new VisWin( cv,  DOMAIN_WIN_SIZE, COLOR_BAR_SIZE );

  // *************** PARENT CONTAINER WINDOW *************
 

  glutInit(&argc, argv);
  glutInitDisplayMode (GLUT_DEPTH | GLUT_DOUBLE | GLUT_RGB);
  glutInitWindowSize (mainWidth, mainHeight); 
  glutInitWindowPosition (100, 100);
  mainWindow = glutCreateWindow ("vees");
  
  glutDisplayFunc(mainDisplay); 
  glutReshapeFunc(mainReshape);

  // ************ DOMAINVIEWER WINDOW *****************
 
  // big sub window working fine. 
   i = glutCreateSubWindow(mainWindow, GAP, GAP, 
				   DOMAIN_WIN_SIZE, DOMAIN_WIN_SIZE);
   // void GlutSubWindow::init(  int superWidth, int superHeight, int placeX, 
   //  int placeY, int width, int height, int num )
   subWins[DOMAIN_VIEWER].init(  mainWidth, mainHeight, GAP, GAP,
		     DOMAIN_WIN_SIZE, DOMAIN_WIN_SIZE, i ); 
  
   glutReshapeFunc(domainReshape);
   glutDisplayFunc(domainDisplay);
   glutKeyboardFunc(domKeyboard);

   glutMotionFunc( mainMouseMotion );
   glutMouseFunc( mainMouse );
   

   glui = GLUI_Master.create_glui( "Rendering Controls" );

   //glui= GLUI_Master.create_glui_subwindow (mainWindow, 
   //									GLUI_SUBWINDOW_LEFT);

  // **** draw methods *****************


  drawTypeList = glui->add_listbox( "Draw Function:", &drawType, DRAW_TYPE,
						 	gluiCallback );
  for( i=0; i < numDrawEleNames; i++ ) {
  
     drawTypeList->add_item( i, drawEleNames[i] );
  } 

  // ********* scalar color range ************
  colorMax = glui->add_spinner("Max Color", GLUI_SPINNER_FLOAT,&maxColor,
						      COLOR_MAX,gluiCallback );
  colorMin = glui->add_spinner("Min Color", GLUI_SPINNER_FLOAT,&minColor,
					      	   COLOR_MIN, gluiCallback);
  colorMax->set_float_limits( 0.0, 1.0, GLUI_LIMIT_CLAMP);
  colorMin->set_float_limits( 0.0, 1.0, GLUI_LIMIT_CLAMP);

  

  invertColorRange = glui->add_checkbox("Invert Color Range", &colorInvert,
			      		INVERT_COLOR_RANGE,  gluiCallback );


  // ****** hide/show elements ****
 
  //domain viewer finds element type names
  dv->populateElementShowList( names, numTypes ); 
  eleChooserRollout = glui->add_rollout("    Hide Element Type    ");
 
  // go through loop creating checkboxes, and adding them to panel
  // we need a starting range, say 100. then handler says, for anything 
  // over 100, call the domaniviewer handler.
  
  for( i = 0; i < numTypes; i++ ) {
	nameBools[i] = 0;
   	glui->add_checkbox_to_panel( eleChooserRollout,eleStrings[names[i]], 
   		      &nameBools[i], ( SHOW_HIDE_ELEMENT_TYPE + i),
						  	 gluiCallback );
	//anything >=  SHOW_HIDE_ELEMENT_TYPE means call 
    // toggleShow function
  }
   
  selectElement =
   glui->add_edittext( "Select Element:", GLUI_EDITTEXT_INT, &selectedElement,
		       SELECT_ELEMENT, gluiCallback);

  glui->set_main_gfx_window( subWins[0].getWinNum() ); //tell the glui where 
                                           //the openGL window is

  // We register the idle callback with GLUI, *not* with GLUT 
  GLUI_Master.set_glutIdleFunc( myGlutIdle ); 

  // ********* COMMIT STEP LOADER GUI ***************** 
  // widgets loaded vertically unless add_column
  glui2 = GLUI_Master.create_glui( "Load Step",0,-5,-5 );

  // glui2->add_edittext("End",GLUI_EDITTEXT_INT, 
  //		    	  &commitStep, SET_COMMIT_STEP, gluiCallback);
  glui2->add_button("<<", BACKWARDS_MOVIE, gluiCallback);

  glui2->add_edittext("File",GLUI_EDITTEXT_TEXT, 
			  &filename, SET_FILE_NAME, gluiCallback);
  glui2->add_column(false);


  glui2->add_button("<", DECREMENT_STEP, gluiCallback);
  glui2->add_checkbox("Save Image", &saveImage);
  glui2->add_checkbox("Save as VTK", &saveVTK);
  glui2->add_checkbox("Save Rotate", &saveRotation);
  glui2->add_checkbox("Save Stress", &saveStress);
  glui2->add_checkbox("Save Eigentensor", &saveEigenTensor);
  glui2->add_checkbox("Save Eigen Lode Angle", &saveEigenLode);
  glui2->add_column(false);

  glui2->add_button("Load", LOAD_STEP, gluiCallback);
 
  glui2->add_edittext("",GLUI_EDITTEXT_INT, 
		    	  &commitStep, SET_COMMIT_STEP, gluiCallback);
  glui2->add_column(false);

  glui2->add_button(">", INCREMENT_STEP, gluiCallback);


  glui2->add_column(false);
  glui2->add_button(">>", MOVIE , gluiCallback);
  glui2->add_edittext("End Step",GLUI_EDITTEXT_INT, 
		    	  &endStep, SET_COMMIT_STEP, gluiCallback);





  
 
 
  //glui2->add_column(false);


  glui2->set_main_gfx_window( subWins[0].getWinNum() );
  
  // ********* ELEMENT VIEWER SUBWINDOW ************
   i = glutCreateSubWindow(mainWindow, (GAP + GAP + DOMAIN_WIN_SIZE) , GAP, 
				 ELEMENT_WIN_SIZE, ELEMENT_WIN_SIZE);
  // void GlutSubWindow::init(  int superWidth, int superHeight, int placeX, 
   //  int placeY, int width, int height, int num )   
   subWins[ELE_VIEWER].init(  mainWidth, mainHeight,
			    (GAP + GAP + DOMAIN_WIN_SIZE) , 
			   GAP, ELEMENT_WIN_SIZE, ELEMENT_WIN_SIZE, i);

   glutReshapeFunc(domainReshape);
   glutDisplayFunc(elementDisplay);
   glutKeyboardFunc(eleKeyboard);
   //glutMotionFunc( mainMouseMotion );
   glutMouseFunc( elementMouse );
   
   // ********* INTGRATION POINT SUBWINDOW *************
    i = glutCreateSubWindow(mainWindow, (GAP + GAP + DOMAIN_WIN_SIZE) , 
				    GAP+ GAP +ELEMENT_WIN_SIZE, 
				ELEMENT_WIN_SIZE, ELEMENT_WIN_SIZE);

    subWins[INTGR_PT_VIEWER].init(  mainWidth, mainHeight,
		  (GAP + GAP + DOMAIN_WIN_SIZE) , 
		   GAP+ GAP + ELEMENT_WIN_SIZE , 
		   ELEMENT_WIN_SIZE, ELEMENT_WIN_SIZE, i);
   glutReshapeFunc(domainReshape);
   glutDisplayFunc(glyphDisplay);	
   glutKeyboardFunc(intgrKeyboard);

 
   //************* GAUSS Detail GUI ************
   /*  glui4 = GLUI_Master.create_glui( "Integr. Pt. Visualization" );
   gaussDrawTypeList = glui4->add_listbox( "     Draw Function:", 
					   &gaussDrawType, GAUSS_DRAW_TYPE,
					   gluiCallback );
   for( i = 0; i < ( NUM_GAUSS_GLYPH_TYPES-1); i++ ) {
     // gaussDrawTypeList->add_item( i,  gaussColorSchemeNames[i] );
     }*/


     // ********** COLOR BAR SUBWINDOW *************
    i = glutCreateSubWindow(mainWindow, GAP , 
	        GAP+ GAP +DOMAIN_WIN_SIZE, 
	      	DOMAIN_WIN_SIZE, COLOR_BAR_SIZE);

    // void init(  int superWidth, int superHeight, int placeX, 
    //int placeY, int width, int height, int num );
    subWins[COLOR_BAR_VIEWER].init(  mainWidth, mainHeight,
				     GAP , GAP+ GAP +DOMAIN_WIN_SIZE,
				     DOMAIN_WIN_SIZE, COLOR_BAR_SIZE, i);
    glutReshapeFunc(domainReshape);
    glutDisplayFunc(colorDisplay);
    
    // ******** COLOR GLUI **********
    glui3 = GLUI_Master.create_glui_subwindow( mainWindow,
					       GLUI_SUBWINDOW_BOTTOM );
    
    glui3->set_main_gfx_window( subWins[0].getWinNum() ); //tell the glui where
    //the openGL window is 
    
    colorRangeNum = glui3->add_edittext("Color Map Range +/- ", 
		 GLUI_EDITTEXT_FLOAT,&colorMapRange, COLOR_MAP, gluiCallback ); 
    GLUI_RadioGroup * group = glui3->add_radiogroup( &logColorScale,
                                                   COLOR_MAP, gluiCallback );
    glui3->add_radiobutton_to_group(group, "Linear scale");
    glui3->add_radiobutton_to_group(group, "Log scale");
    
    glui3->add_column(false);

    // scalar coloring methods 
    scalarColorList =  glui3->add_listbox( "Function:", &colorType,
					COLOR_TYPE, gluiCallback );
    
    for( i = 0; i < (NUM_COLOR_SCHEMES-1); i++ ) {
      scalarColorList->add_item( i, colorSchemeNames[i] );
    }

  glutMainLoop();

}

 // write the image to file 
void GlutWin::writePPM() {
 
  FILE    *the_file;
  char filename[100];
  int x,y,width,height;
  unsigned char* rgb_data; 
  width = mainWidth*3; // r,g,b for each pixel
  height =  mainHeight;

  rgb_data =  new unsigned char[3*width];  
  sprintf(filename,"domain.%d.ppm",commitStep );
  if ( (the_file = fopen(filename, "w")) == NULL ) {
      fprintf( stderr,"Couldn't open %s for writing\n",filename );
      return;
  }
  else 
    fprintf( stderr,"\nOpened %s for writing.\n\n",filename );
  
  fprintf( the_file,"P6\n" );
  fprintf( the_file,"# CREATOR: XV Version 3.10a\n" );
  fprintf( the_file,"%d %d\n",width,height );
  fprintf( the_file,"255\n" );
  
  
  for( y = 0; y < height; y++ ) {// for each row
	glReadPixels(0, height - 1 - y, width, 
				 1, GL_RGB, GL_UNSIGNED_BYTE, rgb_data);
	
	// calculate local window position
	for(x =0;  x < width*3; x++) {
	  fputc( rgb_data[x], the_file );	  
	}
  }
  delete[] rgb_data;
  
  fclose(the_file);
  
  fprintf(stderr,"\nSaved %s \n",filename);
}

void GlutWin::writePPM2() {
 
  FILE    *the_file;
  char filename[100];
  int x,y, width;

  // we need a 2D array
  unsigned char **array;
  unsigned char  *pool;
  unsigned char *curPtr;

  width = currWidth*3;
  
  //(step 1) allocate memory for array of elements of column  
  array = new  unsigned char *[currHeight];
  
  //(step 2) allocate memory for array of elements of each row
  pool = new  unsigned char[width* currHeight];

  // Now point the pointers in the right place
  curPtr = pool;
  for( int i = 0; i < currHeight; i++) {
    *(array + i) = curPtr;
    curPtr += width;
  }

   for( y = 0; y < currHeight; y++ )
     for ( x = 0 ; x < width; x++ ) 
       array[y][x] = 200; // light grey
   
  
  sprintf(filename,"domain.%d.ppm",commitStep );
  if ( (the_file = fopen(filename, "w")) == NULL ) {
      fprintf( stderr,"Couldn't open %s for writing\n",filename );
      return;
  }
  else 
    fprintf( stderr,"\nOpened %s for writing.\n\n",filename );
  
  fprintf( the_file,"P6\n" );
  fprintf( the_file,"# CREATOR: XV Version 3.10a\n" );
  fprintf( the_file,"%d %d\n",currWidth,currHeight );
  fprintf( the_file,"255\n" );

  // for each window, get focus and write array

   glutSetWindow(subWins[DOMAIN_VIEWER].getWinNum() );//switch to elem
   drawDomain();
   subWins[DOMAIN_VIEWER].writeSubwindow(array);  

   glutSetWindow(subWins[ELE_VIEWER].getWinNum() );//switch to elem
   drawElement();
   subWins[ELE_VIEWER].writeSubwindow(array);  

   glutSetWindow(subWins[INTGR_PT_VIEWER].getWinNum() );//switch to elem
   drawGlyph();
   subWins[INTGR_PT_VIEWER].writeSubwindow(array);  

   for( y = 0; y < currHeight; y++ )
     for ( x = 0 ; x < width; x++ ) 
       fputc( array[y][x], the_file );

  
   // clean up
   delete [] *array;
   delete [] array;

   fclose(the_file);
  
   fprintf(stderr,"\nSaved %s \n",filename);
}


void GlutWin::writeSubwindow(int i, unsigned char ** array) {
  // get the width and height of the subwindow
  // get the offset of the subwindow (x,y)
  // note this is upper left

  int heightOffset, leftOffset, x, y;
  unsigned char* rgb_data; 

  rgb_data =  new unsigned char[3*subWins[i].getWidth()];  
  heightOffset =  subWins[i].getY();
  leftOffset = subWins[i].getX()*3;

  
   for( y = 0; y < subWins[i].getHeight(); y++ ) {// for each row
	glReadPixels(0, subWins[i].getHeight() - 1 - y, subWins[i].getWidth()*3, 
				 1, GL_RGB, GL_UNSIGNED_BYTE, rgb_data);
	
	// calculate local window position
	for(x =0;  x <subWins[i].getWidth()*3; x++) {
	  array[y + heightOffset][x + leftOffset] = rgb_data[x];
	}
  }
   
   delete[] rgb_data;
}

int GlutWin::snapshot() {
  return saveImage;
}

