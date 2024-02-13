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
** DomainViewer.cpp
** Wrapper class for OpenSees Domain to provide management for visualization
**/

#include "DomainViewer.h"

//using namespace std;

#include <Vector.h>
#include <map>
DomainViewer::DomainViewer( ) {
 
  theDomain = NULL;
  colorMax = 1.0;
  colorMin = 0.0;
  invert = false;
  maxEleIndex = maxNodeIndex = numGaussPts = 0;
  theFileChannel = NULL;

}

/**
 * instatiate file channel or check for matching filename.
 *  set commit time. update Domain state.
 * update color list(s). should we save the old domain state first?
 * yes. allow undo.
 */ 

int DomainViewer::updateState( const char * filename, int commitStep)
{
  int success;
  //save state
  printf("DomainViewer::updateState starts\n");
  printf("filename %s\n",filename); 
  fflush(stdout);
 
  if (theFileChannel != NULL ) {
    printf("theFileChannel != NULL\n");
    fflush(stdout);
    
    delete theFileChannel;

    printf("deleted\n");
    fflush(stdout);
    //theFileChannel = new File_Channel("VEESDomain.bak", FILE_WRITE );
    //theFileChannel->sendDomain( 0, theDomain );
    //delete theFileChannel;
  }

  // now set the new Domain state
  theFileChannel = new File_Channel( filename, FILE_READ );
  printf(" DomainViewer::updateState, trying theFileChannel->recvDomain\n");
  fflush(stdout);
  success =  theFileChannel->recvDomain( commitStep, theDomain );
    printf(" DomainViewer::updateState, theFileChannel->recvDomain returned %d, should be 0\n" , success);
  fflush(stdout);
  if(success == 0) {
  //reset local data values

 // we also need to recalculate stiffness stuff
 // updateStiffness();
  updateRanges();

 
  setDrawMethodForAll( drawMethod );

  }
  return success;
}



/**
 * Volume dimensions needed by VisWin for centering and scaling
 * in window
 **/
double DomainViewer::getCenterAndMaxDim( double center[3] ) {
  double max = 0;
  Vector v = theDomain->getPhysicalBounds();
  // order: minX, minY,minZ,maxX,maxY,maxZ

  center[X] = (v[X] + v[X+3])/2;
  center[Y] =  (v[Y] + v[Y+3])/2;
  center[Z] =  (v[Z] + v[Z+3])/2;
  if ( fabs(v[X+3] - v[X]) > max )
	max =  fabs(v[X+3] -v[X]);
  if ( fabs(v[Y+3] - v[Y]) > max )
	max =  fabs(v[Y+3] - v[Y]);
  if ( fabs(v[Z+3] - v[Z]) > max )
	max =  fabs(v[Z+3] - v[Z]);
 
  return max; 

}

/**
 * @returns the number of finite elements in the domain
 */
int DomainViewer::getNumElements() {
  return  maxEleIndex;
} 

/**
 * Set drawing parameters for axis legends, centering
 * pre!: Domain is set
 */
void DomainViewer::setBasicBounds() {
  double max,min;
  max = 0;
  min = 1000000000;
  Vector v = theDomain->getPhysicalBounds();
  // order: minX, minY,minZ,maxX,maxY,maxZ

  center[X] = (v[X] + v[X+3])/2;
  center[Y] =  (v[Y] + v[Y+3])/2;
  center[Z] =  (v[Z] + v[Z+3])/2;
  if ( fabs(v[X+3] - v[X]) > max )
	max =  fabs(v[X+3] -v[X]);
  if ( fabs(v[Y+3] - v[Y]) > max )
	max =  fabs(v[Y+3] - v[Y]);
  if ( fabs(v[Z+3] - v[Z]) > max )
	max =  fabs(v[Z+3] - v[Z]);

  if ( fabs(v[X+3] - v[X]) < min )
	min =  fabs(v[X+3] -v[X]);
  if ( fabs(v[Y+3] - v[Y]) < min )
	min =  fabs(v[Y+3] - v[Y]);
  if ( fabs(v[Z+3] - v[Z]) < min )
	min =  fabs(v[Z+3] - v[Z]);

  compassSize = min/3;
  
  compassLoc[X] = v[X] - min/10; //shift compass slightly off origin
  compassLoc[Y] = v[Y] - min/10;
  compassLoc[Z] = v[Z] - min/10;

}

/**
 * which values within 0..1 range would we like to view? */
void DomainViewer::setColorRange( double min, double max, int invertRange ) {
  bool invert;
  if ( invertRange )
	invert = true;
  else invert = false;
  ranges[currRange].setRange( min, max, invert );

}

/**
 * set domain and allocate scalar arrays for 
 * visualizing scalars, specifically color. Note: setting vis
 * type is element-oriented but actual data is node-oriented
 **/
void DomainViewer::setDomain( Domain  *d ) {

  theDomain = d;

  drawStyle = BY_NODE;
  rg = new ReynoldsGlyph();


  initRanges(); 

  currRange = BY_NODE_NUM;

  drawElementTypes[DRAW_WIRE] = new DrawWireFrame();
  drawElementTypes[DRAW_WIRE_DISPLACED] = new DrawDisplacedWireFrame();
  drawElementTypes[DRAW_SOLID] = new DrawSolid();

  drawElementTypes[DRAW_SOLID_DISPLACED] = new DrawSolidDisplaced();
  drawElementTypes[HIDE_ELEMENT] = new HideElement();
  drawElementTypes[DRAW_GAUSS_PTS] = new DrawGaussPoints();

  drawElementTypes[DRAW_STIFFNESS] = new DrawStiffness();
  //assign all elements a default drawing scheme
  drawMethods = new DrawElement*[maxEleIndex];

 
  setDrawMethodForAll( DRAW_WIRE);


  setBasicBounds();
  //  saveRegularVolumeStiffness( 1499 );

}


/**
 * Allocate ColorRanges and data (color maps)
 * Initialize data: BY_NODE_NUM, BY_ELEMENT_NUM, BY_ELEMENT_TYPE, 
 * BY_NODE_DISPLACEMENT, BY_MEAN_STRESS
 */
void DomainViewer::initRanges() {
  int i,num, type, prev;
  double prev_val;
  Element * theElement;
  Node * theNode;
  double * data0;
  double * data1;
  double * data2;
  double * data3;

  NodeIter &theNodes = theDomain->getNodes();  
  num = 0;
  while (( theNode = theNodes() ) != 0) {
	i = theNode->getTag();
	if (i > num)
	  num = i;
  }
  maxNodeIndex = num+1;

  //------ Element & Material Initialization ----
  ElementIter &theElements = theDomain->getElements();
  num = 0;
  while (( theElement = theElements() ) != 0) {
	i = theElement->getTag();
	if (i > num)
	  num = i;
  }
  maxEleIndex = num+1;


  // -------- BY_NODE_NUM (0 to 1) --------------
  ranges[BY_NODE_NUM].setDataName( BY_NODE_NUM );


  try {
   data0 = new double[maxNodeIndex];
  }
  catch( std::bad_alloc) {
    printf("DomainViewer::initRanges, Unable to allocate node color array 0\n");
    exit(0);
  }

  for (i = 0; i < (maxNodeIndex); i++ )
    data0[i] = ((double) i)/((double)maxNodeIndex);
  ranges[BY_NODE_NUM].setData( data0, maxNodeIndex );

  //-------- BY_ELEMENT_NUM (0 to 1)  --------------
  try {
    data1 = new double[maxEleIndex];
  }
 catch( std::bad_alloc) {
    printf("DomainViewer::initRanges, Unable to allocate element color array 1\n");
    exit(0);
  }
  for (i = 0; i < (maxEleIndex); i++ )
    data1[i] = ((double) i)/((double)maxEleIndex);
  ranges[BY_ELEMENT_NUM].setData( data1, maxEleIndex );
  ranges[BY_ELEMENT_NUM].setDataName( BY_ELEMENT_NUM );


  //-------- BY_ELEMENT_TYPE  (0 to 1) -----------
  try {
    data2 = new double[maxEleIndex];
  }
  catch( std::bad_alloc) {
    printf("DomainViewer::initRanges, Unable to allocate color array 2\n");
    exit(0);
  }
  prev_val = prev = -1;
  
  theElements = theDomain->getElements();
  //each color is a random number, seeded by element's class tag
  while ( ( theElement = theElements() ) != 0 ) {
    i = theElement->getTag();
    type = theElement->getClassTag(); 
    
    if ( type == prev ) {
      data2[i] = prev_val;
    }
    else {   //cache the generated color for efficiency
      srand( type );
      prev = type;
      data2[i] = ( (double)rand() )/( (double)RAND_MAX );
      prev_val = data2[i];
    } 
  } //for each element
  ranges[BY_ELEMENT_TYPE].setData( data2, maxEleIndex );
  ranges[BY_ELEMENT_TYPE].setDataName( BY_ELEMENT_TYPE );
 
  //------ BY_NODE_DISPLACEMENT -----------
  try {
    data3 = new double[maxNodeIndex];
  }
  catch( std::bad_alloc ) {
    printf("DomainViewer::initRanges, Unable to allocate color array 3\n");
    exit(0);
  }
  ranges[BY_NODE_DISPLACEMENT].setData(data3,maxNodeIndex);
  ranges[BY_NODE_DISPLACEMENT].setDataName( BY_NODE_DISPLACEMENT );
  if ( setDisplacement( data3 ) == true ) { //displacement has occurred 
    ranges[BY_NODE_DISPLACEMENT].scaleAndSetRange( data3, maxNodeIndex,
						   true, false);
  }

  //---- Gauss Pt Location ----------------
  try {
    gsdraw = new DrawGaussPoints[maxEleIndex];
  }
  catch (std::bad_alloc) {
    printf("DomainViewer::initRanges, Unable to allocate gauss color array \n");
    exit(0);
  }

  //--------- gauss stiffness ---------------
  try {
    stfdraw = new DrawStiffness[maxEleIndex];
  }
  catch (std::bad_alloc) {
     printf("DomainViewer::initRanges, Unable to allocate gauss color array \n");
    exit(0);
  }

  initGauss(); // for both point and stiffness

  //---- BY_MEAN_STRESS ---------------------------
  ranges[BY_MEAN_STRESS].setData( gaussColorPtr, numGaussPts );
  ranges[BY_MEAN_STRESS].setDataName( BY_MEAN_STRESS ); //el by el

  
  if ( setMeanStress() == true ) { //stress has occured

    setMeanStress();
    saveMeanStress();
    ranges[BY_MEAN_STRESS].scaleAndSetRange( gaussColorPtr, 
                             numGaussPts, false, false );
  }

  //---- BY_DEVIATORIC --------------------
  ranges[BY_DEVIATORIC].setData( devGaussColorPtr, numGaussPts );
  ranges[BY_DEVIATORIC].setDataName( BY_DEVIATORIC ); //el by el
  if ( setMeanStress() == true ) { //stress has occured		    
     
    setDeviatoricStress();
    //  saveDeviatoricStress();
    //ranges[BY_DEVIATORIC].scaleAndSetRange( devGaussColorPtr, 
    //                         numGaussPts, false, false );
    printf("setDeviatoricStress();\n");
    fflush(stdout);
    ranges[BY_DEVIATORIC].scaleAndSetRange( devGaussColorPtr, 
                             numGaussPts, false, false );
  }

  //----- BY_STIFFNESS ------------------
  updateStiffness();
  setMinStiffnessEigenvalue();
  //  ranges[BY_STIFFNESS].setData( stiffEigenColorPtr, numGaussPts );
  ranges[BY_STIFFNESS].scaleAndSetRange( stiffEigenColorPtr, 
                             numGaussPts, false, false );
  ranges[BY_STIFFNESS].setDataName( BY_STIFFNESS ); 

  ranges[BY_STIFF_EIGENMODE].setDataName(BY_STIFF_EIGENMODE);
  
}

/**
 * We have loaded a new state for all the finite elements.
 * Assumption: the number and type of elements has not changed.
 * We only need to update changed things such as displacement.
 * We keep the color map the same. There is a seperate method 
 * to update the color map
 */
void DomainViewer::updateRanges() {
  printf("DomainViewer::updateRanges\n");
  if ( setDisplacement( ranges[BY_NODE_DISPLACEMENT].getData() ) == true ) { //displacement has occurred 
  
    // check whether there is an old scale value. If so use it.
    // otherwise create a new scale value
    if ( ranges[BY_NODE_DISPLACEMENT].getBaseValue() == 1.0 ) //default
      ranges[BY_NODE_DISPLACEMENT].scaleAndSetRange(
	   ranges[BY_NODE_DISPLACEMENT].getData(), maxNodeIndex,
				      true, false);
    else {
      printf("\nScaling Node Displacement:");
      ranges[BY_NODE_DISPLACEMENT].scaleDataSet();
    }
  }

  
  bool set;
  set = setMeanStress(); //element by element
  printf("Mean stress set, scaling:");
  ranges[BY_MEAN_STRESS].scaleDataSet();

  printf("set deviatoric stress & scale\n");
  setDeviatoricStress();
  ranges[BY_DEVIATORIC].scaleDataSet();

  printf("setMinStiffness & scale\n");
  setMinStiffnessEigenvalue();
  ranges[BY_STIFFNESS].scaleDataSet();
}

/**
 * Displacement is based on magnitude of the first 2
 * or 3 DOFs, depending on how many coordinates the node has.
 * The assumption is that the first DOFs are displacement.
 */ 
bool DomainViewer::setDisplacement( double * colors) {
  int i,j;
  Node * theNode;
  Vector d,c;
  bool isDisplaced;

  NodeIter & theNodes = theDomain->getNodes();
  isDisplaced = false; //flag; there is some displacement. scale colors.

  for (i = 0; i < maxNodeIndex; i++ ) 
	colors[i] = 0.0;

   while ( ( theNode = theNodes() ) != 0 ) {
  
	c = theNode->getCrds();
	d = theNode->getDisp();
	for( j = 0; ( j < c.Size() && j < d.Size() ); j++ ) {
	  if ( d[j] != 0.0 ) {
		isDisplaced = true;
		colors[theNode->getTag()] += (d[j]*d[j]);
	  }
	}
	if ( colors[theNode->getTag()] > 0.0 )
	  colors[theNode->getTag()] = 
		sqrt( colors[theNode->getTag()] );
  }
  return isDisplaced;
}


void DomainViewer::setDrawMethodForAll( int theDrawMethod ) {
  int i;
 
  drawMethod = theDrawMethod;
  if (drawMethod != DRAW_GAUSS_PTS && drawMethod != DRAW_STIFFNESS )
    for (i = 0; i < maxEleIndex; i++ ) {
      drawMethods[i] =  drawElementTypes[drawMethod];
    }
  else { 
    for (i = 0; i < maxEleIndex; i++ ) {

      if( drawMethod == DRAW_GAUSS_PTS) {
	//	drawMethods[i] =  (DrawElement * )&gsdraw[i];
	drawMethods[i] =  (DrawElement * )&stfdraw[i];
	stfdraw[i].setToStress();
      }
      else if ( drawMethod == DRAW_STIFFNESS ){
	drawMethods[i] =  (DrawElement * )&stfdraw[i];
	stfdraw[i].setToStiffness();
      }
    } // for each element
  }
  printf("DomainViewer::setDrawMethodForAll %d\n",  drawMethod);

}

void DomainViewer::setDrawMethodForOne( int eleNum, int drawMethod ) {
  drawMethods[eleNum] =  drawElementTypes[drawMethod];
}

/** Set the coloring technique for all elements
 * @return whether the color method matches the current draw style,
 * e.g. gauss color scheme for gauss points, node or element for 
 * anything else 
*/
bool DomainViewer::setScalarColorMethodForAll( int colorMethod ) {

  // check whether the color method matches the draw style
  /* if ( colorMethod == BY_MEAN_STRESS && drawStyle != BY_GAUSS_PT ) {
    printf("Set the draw method to gauss points to see mean stress\n");
    return false;
    }*/

  currRange = colorMethod;
 
  printf("DomainViewer::setScalarColorMethodForAll %d\n", colorMethod );
  if (currRange == BY_MEAN_STRESS)
    printf("DomainViewer::setScalarColorMethodForAll to BY_MEAN_STRESS %f \n", ranges[currRange].getBaseValue());

  if (currRange == BY_DEVIATORIC)
    printf("DomainViewer::setScalarColorMethodForAll to BY_DEVIATORIC %f\n",ranges[currRange].getBaseValue());
  switch ( colorMethod ) {
	
  case ( BY_NODE_NUM ): 
      
	drawStyle = BY_NODE;
	break;
  
  case ( BY_ELEMENT_NUM ): 

	drawStyle = BY_ELEMENT;
	break;	
  
  case ( BY_NODE_DISPLACEMENT ): 

	drawStyle = BY_NODE;
	break;

  case ( BY_ELEMENT_TYPE ): // seed random number w/element type
	drawStyle = BY_ELEMENT;   // & take first rand in series for color
	break;
  
  case ( BY_MEAN_STRESS ): {

	drawStyle = BY_GAUSS_PT;
	break;
  }
  case ( BY_DEVIATORIC ): {

	drawStyle = BY_GAUSS_PT;
	break;
  }
  case( BY_STIFF_EIGENMODE) : // same as DRAW_STIFFNESS, just drop thru
  case (BY_STIFFNESS) : {
	drawStyle = BY_GAUSS_PT;
	break;
  }
  default:
	printf("Unknown color scheme\n");
  } // switch
  
  
  //for some reason we need to 
  //ensure we don't use a cached value
  fflush(stdout);
  glFlush();
  return true;
}



/** test function: draw nodes as LINE LOOP */ 
void DomainViewer::drawNodes() {
  if( theDomain != NULL ) {
	NodeIter &theNodesS = theDomain->getNodes();
	Node * aNode;
	
	glBegin(GL_POINTS); 
	while (( aNode=theNodesS() ) != 0) {
	  
	  Vector v = aNode->getCrds(); 
	  if (v.Size() == 3 )
		glVertex3f(v[0], v[1],v[2]);
	  else if (v.Size() == 2 )
		glVertex2f(v[0], v[1] );
	} 
	glEnd();
  }
  else printf("Error in DomainViewer::drawNodes: Domain is NULL\n");
}




/* 
 * Iterate through elements, drawing each using selected draw method.
 * pre!: draw methods assigned to all elements
 */
void DomainViewer::draw() {
 
  Element *theElement;
  int pick;

 
  double * black = new double[1];
  black[0] = 0;

  ElementIter &theElements = theDomain->getElements();
 
  drawCompass();
  
  // save pick state
  glGetIntegerv(GL_RENDER_MODE,&pick); 
  glLineWidth(2); // let's try it!


  while (( theElement = theElements() ) != 0) {

    switch( drawMethods[theElement->getTag()]->getClassTag() ) {
      //case  DRAW_STIFFNESS:
    case  DRAW_GAUSS_PTS: {
    
      if( pick == GL_SELECT ){ // pick element, not gauss point
	drawElementTypes[DRAW_SOLID]->drawElement( 
	       theElement,ranges[currRange].getData(), 
	       ranges[currRange],drawStyle );
	 
	} 
      else { 
       
	if ( currRange == BY_MEAN_STRESS) {
	drawMethods[theElement->getTag()]->drawElement( 
	     theElement,  gaussColor[theElement->getTag()], 
	       ranges[currRange],drawStyle );

	}
	else if ( currRange == BY_DEVIATORIC ) 
	drawMethods[theElement->getTag()]->drawElement( 
	     theElement,  devGaussColor[theElement->getTag()], 
	       ranges[currRange],drawStyle );

      } 

    
    }
    break;
   
    
    case DRAW_STIFFNESS: {
      // only draw if stiffness changed
   
     
      if (drawMethods[theElement->getTag()]->canDraw(theElement)) {
	if( pick == GL_SELECT ){ // pick element, not gauss point
	  drawElementTypes[DRAW_SOLID]->drawElement( 
	       theElement,ranges[currRange].getData(), 
	       ranges[currRange],drawStyle );
	 
	} 
	else { 

	  glLineWidth(1); 

	  drawElementTypes[DRAW_WIRE]->drawElement( 
	       theElement,black, 
	        ranges[currRange],SINGLE_COLOR );
	 
	  if (currRange == BY_DEVIATORIC)
	    drawMethods[theElement->getTag()]->drawElement( 
	       theElement,  
                devGaussColor[theElement->getTag()], 
	       ranges[currRange],drawStyle );
	  else if ( currRange == BY_MEAN_STRESS )
	    drawMethods[theElement->getTag()]->drawElement( 
	       theElement,  
                gaussColor[theElement->getTag()], 
	       ranges[currRange],drawStyle );

	  // for stiffness or eigenmode, use stiffness filter
	  else if ( currRange == BY_STIFFNESS )
	    drawMethods[theElement->getTag()]->drawElement( 
	       theElement,  
               stiffEigenColor[theElement->getTag()], 
	       ranges[currRange],drawStyle );
	  else if (  currRange == BY_STIFF_EIGENMODE ) 
	    drawMethods[theElement->getTag()]->drawElement( 
	       theElement,  
               stiffEigenColor[theElement->getTag()], 
	       ranges[BY_STIFFNESS],  BY_STIFF_EIGENMODE );
	 

	} 
      } // selectable element

    } //DRAW_STIFFNESS
      break;
    

    default:
      drawMethods[theElement->getTag()]->drawElement( 
	     theElement,ranges[currRange].getData(), 
	       ranges[currRange],drawStyle );

    } // end switch

  }// while more Elements
  delete []black;
}

 
/**
 *  show X,Y,Z axis legend with icon in upper left corner */
void  DomainViewer::drawCompass() {
  double base,height, len;
 
  GLUquadric * cylinder;
  cylinder = gluNewQuadric();
  height = compassSize;
  base = compassSize/8;

  glLineWidth(3);

  len = height/3;
 
  Vector v = theDomain->getPhysicalBounds();

  //------- Z axis ---------
 set_color_by_param_shiny( 0.0, 1.0 ); //blue

  
  glPushMatrix();
  
  glTranslatef( compassLoc[X], compassLoc[Y], compassLoc[Z]);
  gluCylinder( cylinder , base/2 , base/2, 
	       height*4 , 5 , 1 );

  /*  //letter Z
  glBegin(GL_LINE_STRIP);
    glVertex3f( base, len, len+height*3.5);
    glNormal3f(.4,.5,.4);
  
    glVertex3f(base ,len,base+height*3.5);
    glNormal3f(.4,.5,.4);
  
    glVertex3f(base,base,len+height*3.5);
    glNormal3f(.4,.5,.4);
  
  
    glVertex3f(base, base, base+height*3.5);
    glNormal3f(.4,.5,.4);
    glEnd();*/
  
  glPopMatrix();
  
  //---------- X axis ---------
 set_color_by_param_shiny( 1.0, 1.0 ); //red

  glPushMatrix();
  glTranslatef( compassLoc[X], compassLoc[Y], compassLoc[Z]);
  
  glRotatef(90, 0.0, 1.0, 0.0); //rotate around Y axis, swap x & z
  gluCylinder( cylinder , base/2 , base/2, 
	       height*4 , 5 , 1 );

  /*  glBegin(GL_LINES); //-- X  char --
    glVertex3f( base,base, len*2+height*3.5);
    glNormal3f(.4,.5,.4);

    glVertex3f(base,len,base*2+height*3.5);
    glNormal3f(.4,.5,.4);

    glVertex3f(base,len,len*2+height*3.5);
    glNormal3f(.4,.5,.4);


    glVertex3f(base,base,base*2+height*3.5);
    glNormal3f(.4,.5,.4);
   

    glEnd();*/

  glPopMatrix();

  //------- Y axis --------
  
  
  set_color_by_param_shiny( 0.5, 1.0 ); //green
  glPushMatrix();
  glTranslatef( compassLoc[X], compassLoc[Y], compassLoc[Z]);
  
  glRotatef(270, 1.0, 0.0, 0.0); //rotate around X axis
  gluCylinder( cylinder , base/2 , base/2, 
	       height*4 , 5 , 1);

  /* glBegin(GL_LINES);
    glVertex3f( base,len,base*4 +height*3.5 );
    glNormal3f(.4,.5,.4);

    glVertex3f( len, len, base*4 + len  +height*3.5);
    glNormal3f(.4,.5,.4);   // the long slash /

	//the short slash `
	            // x  z     y   
    glVertex3f( base,len, base*4 + len +height*3.5 );    
    glNormal3f(.4,.5,.4); 

     
    // mid point of long slash 
    glVertex3f( (base +len)/2, len, base*4 + len/2+height*3.5 );    
    glNormal3f(.4,.5,.4);
    glEnd();*/
 
  glPopMatrix();
  gluDeleteQuadric(cylinder); // clean up

  
}



/**
 * Make draw method be "Hide" but save old drawing method so you can "undo"
 */
void DomainViewer::toggleShowElement( int classTag, int hideElement ) {
  Element * theElement;
  ElementIter &theElements = theDomain->getElements();
  HideElement * prev; 

  while (( theElement = theElements() ) != 0) {
	
	if( hideElement ) {

	  if( theElement->getClassTag() == eleTags[classTag] && 
	   drawMethods[theElement->getTag()]->getClassTag() != HIDE_ELEMENT ) {

		HideElement * hide = new HideElement();
		hide->setLastDrawnAs( drawMethods[theElement->getTag()] );
		drawMethods[theElement->getTag()] = hide;

	  }
	} // if hide
	else { // undo hide 

	  if( theElement->getClassTag() == eleTags[classTag] && 
	   drawMethods[theElement->getTag()]->getClassTag() == HIDE_ELEMENT ) {

		prev = (HideElement *)drawMethods[theElement->getTag()];
   
		drawMethods[theElement->getTag()] = prev->revert();
		delete prev; // no leaks
	  }
	}

  }//while (elements)
 

}



/**
 * GUI would like a list of element types in Domain
 * so usr can show/hide by type. Meanwhile, the Domainviewer
 * needs a hashmap to translate the string to a classtag
 * !pre: names have been allocated by caller
 */
void DomainViewer::populateElementShowList( int names[], int & numTypes )
{
  int tag,i,j, cachedVal;
  int * classes;
  Element *theElement;  
  std::map<const char *, int > eleChoiceTags; //hashmap helps us count strings
  
  ElementIter &theElements = theDomain->getElements();
  cachedVal = -1;
  
  //---- count number of types and build hash table ---
  while (( theElement = theElements() ) != 0) {
	tag = theElement->getClassTag();

	if ( tag != cachedVal ) {
	  cachedVal = tag;

	  for( i = 0; i < numEleTypes; i++ ) {
		if ( eleTags[i] == tag  ) {
		  eleChoiceTags[eleStrings[i]] = tag; // add to hash table
		  i = numEleTypes;
		}
	  }// iterate through list to find right element
	} //if tag != cached
  }

// then allocate

  numTypes = eleChoiceTags.size();


  classes = new int[numTypes]; // classTags

  for( i = 0; i < numTypes; i++ ) {
	names[i] = -1;
	classes[i] = -1;
  }

  cachedVal = -1;

  //-- build the index list which will be used 
  //-- to build the string list (map iterator not working)
  theElements = theDomain->getElements(); // rewind iterator
  while (( theElement = theElements() ) != 0) {
	tag = theElement->getClassTag();
	if ( tag != cachedVal ) {
	  cachedVal = tag;

	  for( i = 0; i < numEleTypes; i++ ) { //find class tag
		if ( eleTags[i] == tag  ) {

		  // does array contain tag? 
		  for( j = 0; j < numTypes; j++ ) {
			if( classes[j] == tag ) // already in list
			  j = numTypes;  

			else if( classes[j] == -1) { // add to list
			  classes[j] = tag;
			  names[j] = i;
			  j = numTypes;
			  
			} 
		  } // does array contain tag? 
		  i = numEleTypes; //finished with 1 index
		}
	  }// find class tag
	} //if tag != cached
  }//for elemants 

  free( classes );
}


const ColorRange& DomainViewer::getColorRange(int colorScheme) {

  if ( colorScheme >= 0 && colorScheme <  NUM_COLOR_SCHEMES )
    return ranges[colorScheme];

  else {
    printf("DomainViewer::getColorRange, unknown color scheme %d requested\n",
	   colorScheme);
    printf("returning current color range\n");
    return ranges[currRange];
  }
  
}
  
/**
 * Update the current color map, if it is mutable.
 * Process: 1. recalculate data (displacement. stress, etc.)
 *          2. assign new range for color values
 *          3. scale the data
 */
void DomainViewer::updateColorMap( bool isLogScale, double dataRange){ 

  if( currRange == BY_NODE_DISPLACEMENT) {
    setDisplacement( ranges[BY_NODE_DISPLACEMENT].getData() );
     ranges[BY_NODE_DISPLACEMENT].setBasis( dataRange, true,isLogScale );
     ranges[BY_NODE_DISPLACEMENT].scaleDataSet();
  }
  else if (currRange == BY_MEAN_STRESS ) {
    printf("DomainViewer::updateColorMap BY_MEAN_STRESS\n");
    setMeanStress();
    //stress is not positve definite
    ranges[BY_MEAN_STRESS].setBasis( dataRange, false, isLogScale );
    ranges[BY_MEAN_STRESS].scaleDataSet();
  }
  else if ( currRange == BY_DEVIATORIC ) {
    setDeviatoricStress();
    ranges[BY_DEVIATORIC].setBasis( dataRange, false, isLogScale );
    ranges[BY_DEVIATORIC].scaleDataSet();
  }
  else if ( currRange == BY_STIFFNESS ) {
    setMinStiffnessEigenvalue();
    ranges[BY_STIFFNESS].setBasis( dataRange, false, isLogScale );
    ranges[BY_STIFFNESS].scaleDataSet();
  }
}

/**
 * Calculate gauss locations and size to draw spheres once.
 * Alllocate gauss pt color arrays (2D mixed length array)
 */
void DomainViewer::initGauss() {
  int i;
  ElementIter &theElements = theDomain->getElements();
  Element * theElement;

  double * currPtr;

  numGaussPts = 0;
  theElements = theDomain->getElements();
  while (( theElement = theElements() ) != 0) {
     if( gsdraw[theElement->getTag()].setElement( theElement ) ) {
       numGaussPts += gsdraw[theElement->getTag()].getNumPts();
       
     }
     stfdraw[theElement->getTag()].setElement( theElement ); // same for stiff
 }

  // allocate gauss color array of pointers
  gaussColor = new double *[maxEleIndex];

  if (!gaussColor ) {
    printf("Unable to allocate color pointers for gauss points\n");
    exit(0);
  }


  gaussColorPtr = new double[numGaussPts];
  if (!gaussColorPtr ) {
    printf("Unable to allocate color array for gauss points\n");
    exit(0);
  }

  for(i = 0; i < numGaussPts; i++ )
    gaussColorPtr[i] = 0.0;

  currPtr = gaussColorPtr;

  theElements = theDomain->getElements();
  while (( theElement = theElements() ) != 0) {
     *(gaussColor + theElement->getTag()) = currPtr;
     currPtr +=  gsdraw[theElement->getTag()].getNumPts(); 
  }

  //---- deviatoric stress ------------
 devGaussColor = new double *[maxEleIndex];

  if (!devGaussColor ) {
    printf("Unable to allocate color pointers for gauss points (mean stress)\n");
    exit(0);
  }


  devGaussColorPtr = new double[numGaussPts];
  if (!devGaussColorPtr ) {
    printf("Unable to allocate color array for gauss points (deviatoric)\n");
    exit(0);
  }

  for(i = 0; i < numGaussPts; i++ )
    devGaussColorPtr[i] = 0.0;

  currPtr = devGaussColorPtr;

  theElements = theDomain->getElements();
  while (( theElement = theElements() ) != 0) {
     *(devGaussColor + theElement->getTag()) = currPtr;
     currPtr +=  gsdraw[theElement->getTag()].getNumPts(); 
  }

  //---------- MIN STIFFNESS -------------

 stiffEigenColor = new double *[maxEleIndex];

  if (!stiffEigenColor ) {
    printf("Unable to allocate color pointers for gauss points (mean stress)\n");
    exit(0);
  }


  stiffEigenColorPtr = new double[numGaussPts];
  if (!stiffEigenColorPtr ) {
    printf("Unable to allocate color array for gauss points (deviatoric)\n");
    exit(0);
  }

  for(i = 0; i < numGaussPts; i++ )
    stiffEigenColorPtr[i] = 0.0;

  currPtr = stiffEigenColorPtr;

  theElements = theDomain->getElements();
  while (( theElement = theElements() ) != 0) {
     *(stiffEigenColor + theElement->getTag()) = currPtr;
     currPtr +=  gsdraw[theElement->getTag()].getNumPts(); 
  }
}

/**
 * Reset elements using new data
 */
void DomainViewer::updateStiffness() {
  int i;
  ElementIter &theElements = theDomain->getElements();
  Element * theElement;

  theElements = theDomain->getElements();
  while (( theElement = theElements() ) != 0) {
     stfdraw[theElement->getTag()].setElement( theElement ); // same for stiff
 }

}




/**
 * Calculate mean stress for every element's gauss pt set
 * @returns true if stress has occured.
 */  
bool DomainViewer::setMeanStress() {
  ElementIter &theElements = theDomain->getElements();
  Element * theElement;

  bool stressed = false;
  while (( theElement = theElements() ) != 0) {
    gsdraw[theElement->getTag()].setElement( theElement );

    if( gsdraw[theElement->getTag()].setMeanStress( 
		   gaussColor[theElement->getTag()] ))

      stressed = true;
  }
  double max = -50000;
  double min = 200000;
  for (int i =0; i < numGaussPts; i++ ) {
    if (gaussColorPtr[i] > max)
      max = gaussColorPtr[i];
    if(gaussColorPtr[i] < min )
      min = gaussColorPtr[i];
  }
  printf("max stress %f min stress %f\n", max , min);
  fflush(stdout);
  return stressed;

}

/**
 * Calculate deviatoric stress for every element's gauss pt set
 * @returns true if stress has occured.
 */  
void DomainViewer::setDeviatoricStress() {
  ElementIter &theElements = theDomain->getElements();
  Element * theElement;

 
  while (( theElement = theElements() ) != 0) {
    gsdraw[theElement->getTag()].setElement( theElement );

    gsdraw[theElement->getTag()].setDeviatoricStress( 
    		 devGaussColor[theElement->getTag()] );

 
  }
  double max = -50000;
  double min = 200000;
  for (int i =0; i < numGaussPts; i++ ) {
    if (devGaussColorPtr[i] > max)
      max = devGaussColorPtr[i];
    if(devGaussColorPtr[i] < min )
      min = devGaussColorPtr[i];
  }
  printf("max deviatoric stress %f min deviatoric stress %f\n", max , min);
  fflush(stdout);

}

void DomainViewer::setMinStiffnessEigenvalue() {
  ElementIter &theElements = theDomain->getElements();
  Element * theElement;

 
  while (( theElement = theElements() ) != 0) {
    stfdraw[theElement->getTag()].setElement( theElement );
    stfdraw[theElement->getTag()].setMinStiffness(
				  stiffEigenColor[theElement->getTag()] );

  }
}

/**
 * Saves gauss point mean stress to a file
 * precondition: mean stress has been calculated */
void DomainViewer::saveMeanStress() {
  ElementIter &theElements = theDomain->getElements();
  Element * theElement;
  Matrix points;
  FILE * theFile;
  char fileName[210];
  int i, j, curr;


  sprintf(fileName, "MeanStress");
  theFile = fopen(fileName, "w");
  curr = 0;
  while (( theElement = theElements() ) != 0) {
     points = gsdraw[theElement->getTag()].getMatrix();
   
     for (i = 0; i < points.noRows(); i++) {
       for (j = 0; j < 3; j++) {
	 // print XYZ
	 fprintf(theFile, "%f ", points(i,j));
       }
       //print gaussColorPtr[curr]
       fprintf(theFile, "%f\n", gaussColorPtr[curr] );
       curr++;
     } // for each gauss point  
  }// for each element

  fclose(theFile);
}
/**
 * Saves gauss point mean stress to a file
 * precondition: deviatoric stress has been calculated */
void DomainViewer::saveDeviatoricStress() {
  ElementIter &theElements = theDomain->getElements();
  Element * theElement;
  Matrix points;
  FILE * theFile;
  char fileName[210];
  int i, j, curr;


  sprintf(fileName, "DeviatoricStress");
  theFile = fopen(fileName, "w");
  curr = 0;
  while (( theElement = theElements() ) != 0) {
     points = gsdraw[theElement->getTag()].getMatrix();
   
     for (i = 0; i < points.noRows(); i++) {
       for (j = 0; j < 3; j++) {
	 // print XYZ
	 fprintf(theFile, "%f ", points(i,j));
       }
       //print gaussColorPtr[curr]
       fprintf(theFile, "%f\n", gaussColorPtr[curr] );
       curr++;
     } // for each gauss point  
  }// for each element

  fclose(theFile);
}

/**
 * File format: 
 * current elastic stiffness
 * current stiffness 
 * yield surface normal 
 * plastic flow normal
 * deformation
 */
void DomainViewer::saveStiffness( int commitStep ) {
  ElementIter &theElements = theDomain->getElements();
  Element * theElement;
  FILE * theFile;
  char fileName[210];
  

  sprintf(fileName, "StiffB%05d.txt", commitStep );
  theFile = fopen(fileName, "w");
  if (!theFile) {
    printf( "DomainViewer::saveStiffness, unable to open file for writing\n");
    return;
  }

  while (( theElement = theElements() ) != 0) {
    stfdraw[theElement->getTag()].printStiffness( theFile );
    // go 0 thru 7 on the brick, get the things & print them out
    

  }
  fprintf(theFile, "EOF%d\n", commitStep);
  fclose(theFile);
}

void DomainViewer::saveVTKdata( int commitStep ) {
  ElementIter &theElements = theDomain->getElements();
  Element * theElement;
  FILE * theFile;
  char fileName[210];
  int numPts = 0; // number of gauss points 
  int i,j;

  sprintf(fileName, "StiffB%05d.vtk", commitStep );
  theFile = fopen(fileName, "w");
  if (!theFile) {
    printf( "DomainViewer::saveVTKdata, unable to open file for writing\n");
    return;
  }

  // print the vtk header
  fprintf(theFile,"# vtk DataFile Version 2.0\n");
  fprintf(theFile, "Gauss point locations min stiffness and eigenmode\n");
  fprintf(theFile,"ASCII\n");
  fprintf(theFile, "DATASET UNSTRUCTURED_GRID\n");
  
 
  // POINTS 27 float
  // print gauss pt locations

  // POINT_DATA 27
  //SCALARS minStiffness float 1 ==> SCALARS dataName dataType numComp
  //LOOKUP_TABLE default

  // count number of gauss points
  while (( theElement = theElements() ) != 0) {
    if (theElement->getClassTag() ==  ELE_TAG_EightNodeBrick )
      numPts += 8;
  }
  fprintf( theFile, "POINTS %d float\n", numPts );

  theElements = theDomain->getElements();
  while (( theElement = theElements() ) != 0) {
    if  (theElement->getClassTag() ==  ELE_TAG_EightNodeBrick ) {
      // note the points must be in isoparametric order
      // not element order
      stfdraw[theElement->getTag()].printPointsVTK( theFile );
    }
  }

 
  fprintf( theFile, "CELLS %d %d\n", numPts/8, ((numPts/8)+numPts)  );
  j = 0;
  for (i = 0; i < (numPts/8); i++) {
    fprintf( theFile, "8 %d %d %d %d %d %d %d %d\n", j,
	     j+1, j+2, j+3, 
       j+4, j+5, j+6, j+7);
    j+=8;
  }

  fprintf( theFile, "CELL_TYPES %d\n", numPts/8);
  for (i = 0; i < numPts/8; i++)
    fprintf( theFile, "12\n");

  fprintf( theFile, "POINT_DATA %d\n", numPts );
  fprintf( theFile, "SCALARS eigenMode float 1\n" );
  fprintf( theFile, "LOOKUP_TABLE default\n");
  theElements = theDomain->getElements();
  while (( theElement = theElements() ) != 0) {
    if  (theElement->getClassTag() ==  ELE_TAG_EightNodeBrick ) {
      stfdraw[theElement->getTag()].printEigenmodeVTK( theFile );
    }
  }
  fprintf( theFile,"\n");

  fclose(theFile);

}

void DomainViewer::saveRotationAccuracy( int commitStep ) {
  ElementIter &theElements = theDomain->getElements();
  Element * theElement;
   FILE * theFile;
  char fileName[210];
  int numElements = 0;
   sprintf(fileName, "Rotation%05d.txt", commitStep );
  theFile = fopen(fileName, "w");
  if (!theFile) {
    printf( "DomainViewer::saveRotationAccuracy, unable to open file for writing\n");
    return;
  }

  
  theElements = theDomain->getElements();

  while (( theElement = theElements() ) != 0) {
    if (theElement->getClassTag() ==  ELE_TAG_EightNodeBrick )
      numElements++;
  }

  fprintf(theFile, "Elements %d\n",  numElements );

  theElements = theDomain->getElements();
      while (( theElement = theElements() ) != 0) {
	if  (theElement->getClassTag() ==  ELE_TAG_EightNodeBrick ) {
   
      stfdraw[theElement->getTag()].printRotationTest( theFile );

    }
  }


   fclose(theFile);
}

/**
 * VTK file format */
void DomainViewer::saveRegularVolumeStiffness( int commitStep) {
  double deltaX, deltaY, deltaZ;
  double loc[3];
  int numX, numY, numZ; // number of spatial increments that keep 
                        // the aspect ratio
  double a;
  int i,j, k; // regular grid locations
  int m,n,p,q; // tensor or matrix entries
  int currElement = 0;

  double incr; // spatial increment, regular grid
  bool zero = false;
  int desiredPoints = 3000; // desired number of points


  FILE * theFile;
  char fileName[210];
  int numPts = 0; // number of gauss points 

  sprintf(fileName, "Stiff_T%05d.vtk", commitStep );
  theFile = fopen(fileName, "w");
  if (!theFile) {
    printf( "DomainViewer::saveRegularVolumeStiffness, unable to open file for writing\n");
    return;
  }

  int dim[4] = {3,3,3,3}; // will get copied

  Tensor t_stiff( 4, dim, 0.0 );

  Matrix m_stiff(6,6);
  Matrix m_rotate(6,6);
  Matrix m_stretch(6,6);


  currElement = 1;

  Vector bounds = theDomain->getPhysicalBounds();
  
  deltaX = fabs(bounds(3) - bounds(0));
  deltaY = fabs(bounds(4) - bounds(1));
  deltaZ = fabs(bounds(5) - bounds(2));

  incr =  minDistBetweenPts()/2.0;
  numX = deltaX/incr;
  numY = deltaY/incr;
  numZ = deltaZ/incr;

  if ( ((numX+1)*(numY+1)*(numZ+1))  < desiredPoints ) {
 
    a = desiredPoints/( deltaX*deltaY*deltaZ ); // a^3 = desiredpts/x*y*z
 
    a = cubeRoot(a); 

    numX = a*deltaX;
    numY = a*deltaY;
    numZ = a*deltaZ;

    incr = deltaX/numX;
   printf("Interpolation scheme 1: number of points %d\n", numX*numY*numZ);
  }
  else
  printf("Interpolation scheme 2: number of points %d\n", numX*numY*numZ);



  fprintf(theFile, "# vtk DataFile Version 2.0\n");
  fprintf( theFile,
	  "3x3x3x3 Stiffness Tensors on regular grid, timestep %d, aneeman@cse.ucsc.edu\n",
	   commitStep);
  fprintf( theFile,"ASCII\n\n");
  fprintf( theFile,"DATASET STRUCTURED_POINTS\n");
  fprintf( theFile, "DIMENSIONS %d %d %d\n",numX+1, numY+1,numZ+1);
  fprintf( theFile,"ORIGIN %f %f %f\n", bounds(0), bounds(1), bounds(2) );
  fprintf( theFile, "SPACING %f %f %f\n", incr, incr, incr);
  fprintf( theFile, "POINT_DATA %d\n", (numX+1)* (numY+1)*(numZ+1));

  //  fprintf("TENSORS fourthOrder float\n");
  fprintf( theFile, "SCALARS fourthOrder float 81\n");
  fprintf( theFile, "LOOKUP_TABLE default\n");

 
  for ( k = 0; k <= numZ; k++ ) { // VTK volume order is Z Y X
    for ( j = 0; j <= numY; j++ ) {
      for ( i = 0; i <= numX; i++) {

  
	zero = false;
	loc[0] = bounds(0) + (i*incr);
	loc[1] = bounds(1) + (j*incr);
	loc[2] = bounds(2) + (k*incr);

	
	// check for inside current element & interpolate
	if (  currElement != -1  &&
	      stfdraw[currElement].insideElement( loc ) ) {
	    stfdraw[currElement].interpolateStiffness( m_stiff, loc );
	}
	
        else {  // find element & interpolate
	  currElement = findElement( loc );
	  if ( currElement != -1  ) {
	    stfdraw[currElement].interpolateStiffness( m_stiff, loc );
	  }

      
	  else { // not found, zero matrix
	    zero = true;

	    for ( m = 0; m < 6; m++ ) {
	      for (n = 0; n < 6; n++ ) {
		m_stiff(k,m) = 0;
	      }
	    }

	  } // no such location 
	} // find element & interpolate
	double T[4][4][4][4];
       	Matrix2TensorSysR4(m_stiff, T);

	if ( i ==3 && j == 3 && k == 3) {
	  t_stiff.print();
	  opserr << m_stiff;
	  
	}
	
	for ( m = 1; m < 4; m++ ) {
	  for (n = 1; n < 4; n++ ) {
	    for ( p = 1; p < 4; p++ ) {
	      for ( q = 1; q < 4; q++ ) {
	      
		 
		 fprintf(theFile, "%.16g ", T[m][n][p][q] );
	      } // q
	    }//p
	    // fprintf(theFile, "\n");
	  }//n
	}//m
	fprintf(theFile,"\n\n\n");
	
      } //k
      fprintf(theFile, "\n");
    }//j
    fprintf(theFile, "\n");
  } // i
  

  fclose(theFile);

  
}

/**
 * VTK file format */
void DomainViewer::saveStressTensor( int commitStep) {
  double deltaX, deltaY, deltaZ;
  double loc[3];
  int numX, numY, numZ; // number of spatial increments that keep 
                        // the aspect ratio
  double a;
  int i,j, k; // regular grid locations
  int m,n; // matrix entries
  int currElement = 0;

  double incr; // spatial increment, regular grid
  bool zero = false;
  int desiredPoints = 200; // desired number of points

  Matrix m_stress(3,3);

  FILE * theFile;
  char fileName[210];
  int numPts = 0; // number of gauss points 



  sprintf(fileName, "Stress%05d.vtk", commitStep );
  theFile = fopen(fileName, "w");
  if (!theFile) {
    printf( "DomainViewer::saveStress, unable to open file for writing\n");
    return;
  }


  currElement = 1;


  Vector bounds = theDomain->getPhysicalBounds();
  
  deltaX = fabs(bounds(3) - bounds(0));
  deltaY = fabs(bounds(4) - bounds(1));
  deltaZ = fabs(bounds(5) - bounds(2));

  incr =  minDistBetweenPts()/2.0;
  numX = deltaX/incr;
  numY = deltaY/incr;
  numZ = deltaZ/incr;

  if ( ((numX+1)*(numY+1)*(numZ+1))  < desiredPoints ) {
 
    a = desiredPoints/( deltaX*deltaY*deltaZ ); // a^3 = desiredpts/x*y*z
 
    a = cubeRoot(a); 

    numX = a*deltaX;
    numY = a*deltaY;
    numZ = a*deltaZ;

    incr = deltaX/numX;
   printf("Interpolation scheme 1: number of points %d\n", numX*numY*numZ);
  }
  else
  printf("Interpolation scheme 2: number of points %d\n", numX*numY*numZ);



 

  fprintf(theFile, "# vtk DataFile Version 3.0\n");
  fprintf( theFile,
	  "3x3 Stress Tensors on regular grid, timestep %d\n",
	   commitStep);
  fprintf( theFile,"ASCII\n");
  fprintf( theFile,"DATASET STRUCTURED_POINTS\n");
  fprintf( theFile, "DIMENSIONS %d %d %d\n",numX+1, numY+1,numZ+1);
  fprintf( theFile,"ORIGIN %f %f %f\n", bounds(0), bounds(1), bounds(2) );
  fprintf( theFile, "SPACING %f %f %f\n", incr, incr, incr);
  fprintf( theFile, "POINT_DATA %d\n", (numX+1)* (numY+1)*(numZ+1));

  //  fprintf("TENSORS fourthOrder float\n");
  // fprintf( theFile, "SCALARS stress_lode_angle float 1\n");
  //fprintf( theFile, "LOOKUP_TABLE default\n"); 

  fprintf( theFile, "TENSORS stress float\n");
 
  for ( k = 0; k <= numZ; k++ ) {
    for ( j = 0; j <= numY; j++ ) {
      for ( i = 0; i <= numX; i++) {

	zero = false;
	loc[0] = bounds(0) + (i*incr);
	loc[1] = bounds(1) + (j*incr);
	loc[2] = bounds(2) + (k*incr);

	
	// check for inside current element & interpolate
	if (  currElement != -1  &&
	      stfdraw[currElement].insideElement( loc ) ) {
	    stfdraw[currElement].interpolateStress( m_stress, loc );
	}
	
        else {  // find element & interpolate
	   currElement = findElement( loc );
	  if ( currElement != -1  ) {
	    stfdraw[currElement].interpolateStress( m_stress, loc );
	  }

      
	  else { // not found, zero matrix
	    zero = true;

	  } // no such location 
	} // find element & interpolate
	
       	// if zero, print zero eigentensor
	if ( zero ) {
	  for (m = 0; m < 3; m++) {
	    for (n = 0; n < 3; n++ ) {
	      fprintf(theFile, "%.16g ", 0.00 );
	    }
	    fprintf(theFile,"\n");
	  }
	  fprintf(theFile,"\n");
	}
	else {
	  for (m = 0; m < 3; m++) {
	    for (n = 0; n < 3; n++ ) { 
	      fprintf(theFile, "%.16g ", m_stress(m,n ));
	    }
	    fprintf(theFile,"\n");
	  }
	  fprintf(theFile,"\n");
	}

      } //k
    
    }//j
  
  } // i
  
  fclose(theFile);

  
}





/* Intensity = amplitude^2 sin^2(delta/2)sin^2 (2 theta)
 * delta = 2 PI/lamda (deviatoric) thickness */

double stressContourColor( double psi) {
  double delta_n;
  double two_PI_over_lamda;
  double polycarbonate;
  double polar_light; // incoming angle in radians
  double angle_diff; // angle between e_vector1 and incoming light

  polar_light = 3.14159/4; //hardcoded to 45 degrees

  /* Stress optic coefficient for polycarbonate is -78/Pa ,
     Pa is m^3/s (so how fast light is retarded) */
  polycarbonate = -78;
  delta_n = psi *polycarbonate; //use constant as (n1-n2)
  two_PI_over_lamda = 12566360;//2PI over wavelength of visible light
  
  delta_n = delta_n*two_PI_over_lamda;
  delta_n = sin( delta_n/16.0 ); // sin (delta/2)
  delta_n = delta_n*delta_n; //sin^2(expr)
  
  
  /* angle_diff =  ((theta[x][z]*3.14159)/180); //radians instead of degrees
  angle_diff = angle_diff -polar_light;
  
  // angle_diff = -3.14159/4;// abs(PI/4) is white, abs(PI/2) black
  angle_diff = sin(2*angle_diff);
  angle_diff = angle_diff*angle_diff;//sin^2(expr)
  
  //  angle_diff = 0.7; //fix to constant
  // return (delta_n*angle_diff);*/
  return delta_n;

}


/**
 * Contour based on unique eigenvector and lode angle. 
 *
*/
void DomainViewer::saveStressLodeAngle( int commitStep) {
  double deltaX, deltaY, deltaZ;
  double loc[3];
  int numX, numY, numZ; // number of spatial increments that keep 
                        // the aspect ratio
  double a;
  int i,j, k; // regular grid locations
  int m,n; // matrix entries
  int currElement = 0;

  double incr; // spatial increment, regular grid
  bool zero = false;
  int desiredPoints = 200; // desired number of points

  Matrix m_stress(3,3);

  FILE * theFile;
  char fileName[210];
  int numPts = 0; // number of gauss points 

  double radius, theta;
  double * e_vectors;  // eigenvectors and eigenvalues of stress
  double * e_values;

  try {
   e_vectors = new double[9];
   e_values = new double[3];
  }
  catch( std::bad_alloc) {
    printf("DomainViewer::saveEigenLodeAngle, unable to allocate eigenvector data\n");
    return;
  } 

 radius = theta = 0.0;



  sprintf(fileName, "StressLode%05d.vtk", commitStep );
  theFile = fopen(fileName, "w");
  if (!theFile) {
    printf( "DomainViewer::saveStress, unable to open file for writing\n");
    return;
  }




  currElement = 1;


  Vector bounds = theDomain->getPhysicalBounds();
  
  deltaX = fabs(bounds(3) - bounds(0));
  deltaY = fabs(bounds(4) - bounds(1));
  deltaZ = fabs(bounds(5) - bounds(2));

  incr =  minDistBetweenPts()/2.0;
  numX = deltaX/incr;
  numY = deltaY/incr;
  numZ = deltaZ/incr;

  if ( ((numX+1)*(numY+1)*(numZ+1))  < desiredPoints ) {
 
    a = desiredPoints/( deltaX*deltaY*deltaZ ); // a^3 = desiredpts/x*y*z
 
    a = cubeRoot(a); 

    numX = a*deltaX;
    numY = a*deltaY;
    numZ = a*deltaZ;

    incr = deltaX/numX;
   printf("Interpolation scheme 1: number of points %d\n", numX*numY*numZ);
  }
  else
  printf("Interpolation scheme 2: number of points %d\n", numX*numY*numZ);



 

  fprintf(theFile, "# vtk DataFile Version 3.0\n");
  fprintf( theFile,
	  "3x3 Stress Tensors on regular grid, timestep %d\n",
	   commitStep);
  fprintf( theFile,"ASCII\n");
  fprintf( theFile,"DATASET STRUCTURED_POINTS\n");
  fprintf( theFile, "DIMENSIONS %d %d %d\n",numX+1, numY+1,numZ+1);
  fprintf( theFile,"ORIGIN %f %f %f\n", bounds(0), bounds(1), bounds(2) );
  fprintf( theFile, "SPACING %f %f %f\n", incr, incr, incr);
  fprintf( theFile, "POINT_DATA %d\n", (numX+1)* (numY+1)*(numZ+1));

  //  fprintf("TENSORS fourthOrder float\n");
  fprintf( theFile, "SCALARS stress_lode_angle float 1\n");
  fprintf( theFile, "LOOKUP_TABLE default\n"); 

  // fprintf( theFile, "TENSORS stress float\n");
 
  for ( k = 0; k <= numZ; k++ ) {
    for ( j = 0; j <= numY; j++ ) {
      for ( i = 0; i <= numX; i++) {

	zero = false;
	loc[0] = bounds(0) + (i*incr);
	loc[1] = bounds(1) + (j*incr);
	loc[2] = bounds(2) + (k*incr);

	
	// check for inside current element & interpolate
	if (  currElement != -1  &&
	      stfdraw[currElement].insideElement( loc ) ) {
	    stfdraw[currElement].interpolateStress( m_stress, loc );
	}
	
        else {  // find element & interpolate
	   currElement = findElement( loc );
	  if ( currElement != -1  ) {
	    stfdraw[currElement].interpolateStress( m_stress, loc );
	  }

      
	  else { // not found, zero matrix
	    zero = true;

	  } // no such location 
	} // find element & interpolate
	
       	// if zero, print zero eigentensor
	if ( zero ) {
	  /*  for (m = 0; m < 3; m++) {
	    for (n = 0; n < 3; n++ ) {
	      fprintf(theFile, "%.16g ", 0.00 );
	    }
	    fprintf(theFile,"\n");
	  }
	  */
	  fprintf(theFile, "%.16g ", -1.00 );
	  fprintf(theFile,"\n");
	}
	else {
	  // for (m = 0; m < 3; m++) {
	  // for (n = 0; n < 3; n++ ) { 
	  //  fprintf(theFile, "%.16g ", m_stress(m,n ));
	  // }
	  // fprintf(theFile,"\n");
	  // }
	
	  rg->lodeAngle(m_stress, radius, theta);
	  if (radius <= 0.0 )
	    fprintf(theFile, "%.16g ", -1.0); // isotropic
	  else {
	    theta = (theta + 30)/60;
	    fprintf(theFile, "%.16g ", theta); // values between zero and one
	  }
	  fprintf(theFile,"\n");
	}

      } //k
    
    }//j
  
  } // i
  
 


// Now save the associated eigenvector.!!!!!!!!!!!!!

fprintf( theFile, "\nVECTORS stress_eigenmode_eigenvector float\n");
  for ( k = 0; k <= numZ; k++ ){
    for ( j = 0; j <= numY; j++ ) {
      for ( i = 0; i <= numX; i++) {

	zero = false;
	loc[0] = bounds(0) + (i*incr);
	loc[1] = bounds(1) + (j*incr);
	loc[2] = bounds(2) + (k*incr);

	
	// check for inside current element & interpolate
	if (  currElement != -1  &&
	      stfdraw[currElement].insideElement( loc ) ) {
	    stfdraw[currElement].interpolateStress( m_stress, loc );

	}
	
        else {  // find element & interpolate
	   currElement = findElement( loc );
	  if ( currElement != -1  ) {
	    stfdraw[currElement].interpolateStress( m_stress, loc );

	  }

      
	  else { // not found, zero matrix
	    zero = true;

	  } // no such location 
	} // find element & interpolate
	
       	// if zero, print zero eigentensor
	if ( zero ) {
	  /*  for (m = 0; m < 3; m++) {
	    for (n = 0; n < 3; n++ ) {
	      fprintf(theFile, "%.16g ", 0.00 );
	    }
	    fprintf(theFile,"\n");
	  }
	  */
	  fprintf(theFile, "%.16g %.16g %.16g", 0.0,0.0,0.0 );
	  fprintf(theFile,"\n");
	}
	else {

	  
	  rg->lodeAngle(m_stress, radius, theta);
	  if (radius <= 0.0 ) {	  
	    fprintf(theFile, "%.16g %.16g %.16g", 0.0,0.0,0.0 );
	    //fprintf(theFile,"\n");
	  }
	   
	  else {
	    theta = (theta + 30)/60; // it is now zero to 1.
                                     // 0 to 1/3: TX Compression
	                             // 1/3 to 2/3 SHEAR
                                     // 2/3 to 1.0 TX tension
	  
	    if ( rg->calcSymmEigens( m_stress, 3, e_vectors, e_values ) != 0 ) 
	      {
		printf(" DomainViewer::saveStress, calculating eigenvectors failed\n");
		return;
	      }                 
	    // eigenvalues are sorted ascending!!

	    if ( theta <= 0.3333333333333 ) { 
	      // Triaxial compression use smallest eigenvalue (-2,+1,+1)
	      fprintf(theFile, "%.16g %.16g %.16g", e_vectors[0],
		      e_vectors[1], e_vectors[2]);
	    }
         
           
	    else if ( theta > 0.3333333333333 &&
		      theta <= 0.6666666666666) {
	      // pure shear, use middle eigenvalue (-1, 0, +1)
	      fprintf(theFile, "%.16g %.16g %.16g", e_vectors[3],
		      e_vectors[4], e_vectors[5]);

	    }
	    
            else {
	      // Triaxial tension, eigenvalues -1,-1,2 use biggest
	      fprintf(theFile, "%.16g %.16g %.16g", e_vectors[6],
		      e_vectors[7], e_vectors[8]);

	   
	    }

	  }
	  fprintf(theFile,"\n");
	}

      } //k
    
    }//j
  
  } //i


// clean up 
 fclose(theFile);
  delete []e_vectors;
  delete []e_values;
  
}
/**
 * Lode angle of stiffness, eigenvector that goes with it.
 *
*/
void DomainViewer::saveEigenLodeAngle( int commitStep) {
  double deltaX, deltaY, deltaZ;
  double loc[3];
  int numX, numY, numZ; // number of spatial increments that keep 
                        // the aspect ratio
  double a;
  int i,j, k; // regular grid locations
  int m,n; // matrix entries
  int currElement = 0;

  double incr; // spatial increment, regular grid
  bool zero = false;
  int desiredPoints = 200; // desired number of points

  Matrix m_stress(3,3);

  FILE * theFile;
  char fileName[210];
  int numPts = 0; // number of gauss points 

  double * e_vectors;
  double * e_values;

  try {
   e_vectors = new double[9];
   e_values = new double[3];
  }
  catch( std::bad_alloc) {
    printf("DomainViewer::saveEigenLodeAngle, unable to allocate eigenvector data\n");
    return;
  }


  double radius, theta;
  radius = theta = 0.0;

  sprintf(fileName, "EigenLode%05d.vtk", commitStep );
  theFile = fopen(fileName, "w");
  if (!theFile) {
    printf( "DomainViewer::saveEigenLodeAngle, unable to open file for writing\n");
    return;
  }

  
  currElement = 1;


  Vector bounds = theDomain->getPhysicalBounds();
  
  deltaX = fabs(bounds(3) - bounds(0));
  deltaY = fabs(bounds(4) - bounds(1));
  deltaZ = fabs(bounds(5) - bounds(2));

  incr =  minDistBetweenPts()/2.0;
  numX = deltaX/incr;
  numY = deltaY/incr;
  numZ = deltaZ/incr;

  if ( ((numX+1)*(numY+1)*(numZ+1))  < desiredPoints ) {
 
    a = desiredPoints/( deltaX*deltaY*deltaZ ); // a^3 = desiredpts/x*y*z
 
    a = cubeRoot(a); 

    numX = a*deltaX;
    numY = a*deltaY;
    numZ = a*deltaZ;

    incr = deltaX/numX;
   printf("Interpolation scheme 1: number of points %d\n", numX*numY*numZ);
  }
  else
  printf("Interpolation scheme 2: number of points %d\n", numX*numY*numZ);



 

  fprintf(theFile, "# vtk DataFile Version 3.0\n");
  fprintf( theFile,
	  "3x3 Stress Tensors on regular grid, timestep %d\n",
	   commitStep);
  fprintf( theFile,"ASCII\n");
  fprintf( theFile,"DATASET STRUCTURED_POINTS\n");
  fprintf( theFile, "DIMENSIONS %d %d %d\n",numX+1, numY+1,numZ+1);
  fprintf( theFile,"ORIGIN %f %f %f\n", bounds(0), bounds(1), bounds(2) );
  fprintf( theFile, "SPACING %f %f %f\n", incr, incr, incr);
  fprintf( theFile, "POINT_DATA %d\n", (numX+1)* (numY+1)*(numZ+1));

  //  fprintf("TENSORS fourthOrder float\n");
  fprintf( theFile, "SCALARS eigenmode_lode_angle float 1\n");
  fprintf( theFile, "LOOKUP_TABLE default\n"); 

  // 
  for ( k = 0; k <= numZ; k++ ) {
    for ( j = 0; j <= numY; j++ ) {
      for ( i = 0; i <= numX; i++) {

	zero = false;
	loc[0] = bounds(0) + (i*incr);
	loc[1] = bounds(1) + (j*incr);
	loc[2] = bounds(2) + (k*incr);

	
	// check for inside current element & interpolate
	if (  currElement != -1  &&
	      stfdraw[currElement].insideElement( loc ) ) {
	    stfdraw[currElement].interpolateEigenTensor( m_stress, loc );

	}
	
        else {  // find element & interpolate
	   currElement = findElement( loc );
	  if ( currElement != -1  ) {
	    stfdraw[currElement].interpolateEigenTensor( m_stress, loc );

	  }

      
	  else { // not found, zero matrix
	    zero = true;

	  } // no such location 
	} // find element & interpolate
	
       	// if zero, print zero eigentensor
	if ( zero ) {
	  /*  for (m = 0; m < 3; m++) {
	    for (n = 0; n < 3; n++ ) {
	      fprintf(theFile, "%.16g ", 0.00 );
	    }
	    fprintf(theFile,"\n");
	  }
	  */
	  fprintf(theFile, "%.16g ", -1.00 );
	  fprintf(theFile,"\n");
	}
	else {
	  // for (m = 0; m < 3; m++) {
	  // for (n = 0; n < 3; n++ ) { 
	  //  fprintf(theFile, "%.16g ", m_stress(m,n ));
	  // }
	  // fprintf(theFile,"\n");
	  // }
	  
	  rg->lodeAngle(m_stress, radius, theta);
	  if (radius <= 0.0 )
	    fprintf(theFile, "%.16g ", -1.0);
	  else {
	    theta = (theta + 30)/60;
	    fprintf(theFile, "%.16g ", theta);
	  }
	  fprintf(theFile,"\n");
	}

      } //k
    
    }//j
  
  } // i
  
  // second iteration. Calculate lode angle and eigenvectors/eigenvalues
  // depending on lode angle, write an eigenvector to file.
  //calcSymmEigens(  Matrix  theMatrix, int leadingDim, 
  //				   double * a (1x9), 
  //                               a[k] =  e_vectors(i,j);
  //double * w (1x3, e_values sorted ascending))

  // if lode angle = -1, print 0,0,0
  // otherwise
  fprintf( theFile, "\nVECTORS stiff_eigenmode_eigenvector float\n");
  for ( k = 0; k <= numZ; k++ ){
    for ( j = 0; j <= numY; j++ ) {
      for ( i = 0; i <= numX; i++) {

	zero = false;
	loc[0] = bounds(0) + (i*incr);
	loc[1] = bounds(1) + (j*incr);
	loc[2] = bounds(2) + (k*incr);

	
	// check for inside current element & interpolate
	if (  currElement != -1  &&
	      stfdraw[currElement].insideElement( loc ) ) {
	    stfdraw[currElement].interpolateEigenTensor( m_stress, loc );

	}
	
        else {  // find element & interpolate
	   currElement = findElement( loc );
	  if ( currElement != -1  ) {
	    stfdraw[currElement].interpolateEigenTensor( m_stress, loc );

	  }

      
	  else { // not found, zero matrix
	    zero = true;

	  } // no such location 
	} // find element & interpolate
	
       	// if zero, print zero eigentensor
	if ( zero ) {
	  /*  for (m = 0; m < 3; m++) {
	    for (n = 0; n < 3; n++ ) {
	      fprintf(theFile, "%.16g ", 0.00 );
	    }
	    fprintf(theFile,"\n");
	  }
	  */
	  fprintf(theFile, "%.16g %.16g %.16g", 0.0,0.0,0.0 );
	  fprintf(theFile,"\n");
	}
	else {

	  
	  rg->lodeAngle(m_stress, radius, theta);
	  if (radius <= 0.0 ) {	  
	    fprintf(theFile, "%.16g %.16g %.16g", 0.0,0.0,0.0 );
	    //fprintf(theFile,"\n");
	  }
	   
	  else {
	    theta = (theta + 30)/60; // it is now zero to 1.
                                     // 0 to 1/3: TX Compression
	                             // 1/3 to 2/3 SHEAR
                                     // 2/3 to 1.0 TX tension
	  
	    if ( rg->calcSymmEigens( m_stress, 3, e_vectors, e_values ) != 0 ) 
	      {
		printf(" DomainViewer::saveEigenLodeAngle, calculating eigenvectors failed\n");
		return;
	      }                 
	    // eigenvalues are sorted ascending!!

	    if ( theta <= 0.3333333333333 ) { 
	      // Triaxial compression use smallest eigenvalue (-2,+1,+1)
	      fprintf(theFile, "%.16g %.16g %.16g", e_vectors[0],
		      e_vectors[1], e_vectors[2]);
	    }
           
           
	    else if ( theta > 0.3333333333333 &&
		      theta <= 0.6666666666666) {
	      // pure shear, use middle eigenvalue (-1, 0, +1)
	      fprintf(theFile, "%.16g %.16g %.16g", e_vectors[3],
		      e_vectors[4], e_vectors[5]);

	    }
	    
            else {
	      // Triaxial tension, efprintf( theFile, "\nVECTORS stiff_eigenmode_eigenvector float\n");
  for ( k = 0; k <= numZ; k++ ){
    for ( j = 0; j <= numY; j++ ) {
      for ( i = 0; i <= numX; i++) {

	zero = false;
	loc[0] = bounds(0) + (i*incr);
	loc[1] = bounds(1) + (j*incr);
	loc[2] = bounds(2) + (k*incr);

	
	// check for inside current element & interpolate
	if (  currElement != -1  &&
	      stfdraw[currElement].insideElement( loc ) ) {
	    stfdraw[currElement].interpolateEigenTensor( m_stress, loc );

	}
	
        else {  // find element & interpolate
	   currElement = findElement( loc );
	  if ( currElement != -1  ) {
	    stfdraw[currElement].interpolateEigenTensor( m_stress, loc );

	  }

      
	  else { // not found, zero matrix
	    zero = true;

	  } // no such location 
	} // find element & interpolate
	
       	// if zero, print zero eigentensor
	if ( zero ) {
	  /*  for (m = 0; m < 3; m++) {
	    for (n = 0; n < 3; n++ ) {
	      fprintf(theFile, "%.16g ", 0.00 );
	    }
	    fprintf(theFile,"\n");
	  }
	  */
	  fprintf(theFile, "%.16g %.16g %.16g", 0.0,0.0,0.0 );
	  fprintf(theFile,"\n");
	}
	else {

	  
	  rg->lodeAngle(m_stress, radius, theta);
	  if (radius <= 0.0 ) {	  
	    fprintf(theFile, "%.16g %.16g %.16g", 0.0,0.0,0.0 );
	    //fprintf(theFile,"\n");
	  }
	   
	  else {
	    theta = (theta + 30)/60; // it is now zero to 1.
                                     // 0 to 1/3: TX Compression
	                             // 1/3 to 2/3 SHEAR
                                     // 2/3 to 1.0 TX tension
	  
	    if ( rg->calcSymmEigens( m_stress, 3, e_vectors, e_values ) != 0 ) 
	      {
		printf(" DomainViewer::saveEigenLodeAngle, calculating eigenvectors failed\n");
		return;
	      }                 
	    // eigenvalues are sorted ascending!!

	    if ( theta <= 0.3333333333333 ) { 
	      // Triaxial compression use smallest eigenvalue (-2,+1,+1)
	      fprintf(theFile, "%.16g %.16g %.16g", e_vectors[0],
		      e_vectors[1], e_vectors[2]);
	    }
           
           
	    else if ( theta > 0.3333333333333 &&
		      theta <= 0.6666666666666) {
	      // pure shear, use middle eigenvalue (-1, 0, +1)
	      fprintf(theFile, "%.16g %.16g %.16g", e_vectors[3],
		      e_vectors[4], e_vectors[5]);

	    }
	    
            else {
	      // Triaxial tension, eigenvalues -1,-1,2 use biggest
	      fprintf(theFile, "%.16g %.16g %.16g", e_vectors[6],
		      e_vectors[7], e_vectors[8]);

	   
	    }

	  }
	  fprintf(theFile,"\n");
	}

      } //k
    
    }//j
  
  } //iigenvalues -1,-1,2 use biggest
	      fprintf(theFile, "%.16g %.16g %.16g", e_vectors[6],
		      e_vectors[7], e_vectors[8]);

	   
	    }

	  }
	  fprintf(theFile,"\n");
	}

      } //k
    
    }//j
  
  } //i

  fclose(theFile);
  delete []e_vectors;
  delete []e_values;
  
}



/**
 * VTK file format */
void DomainViewer::saveEigenTensor( int commitStep) {
  double deltaX, deltaY, deltaZ;
  double loc[3];
  int numX, numY, numZ; // number of spatial increments that keep 
                        // the aspect ratio
  double a;
  int i,j, k; // regular grid locations
  int m,n; // matrix entries
  int currElement = 0;

  double incr; // spatial increment, regular grid
  bool zero = false;
  int desiredPoints = 200; // desired number of points

  Matrix m_eigen(3,3);

  FILE * theFile;
  char fileName[210];
  int numPts = 0; // number of gauss points 



  sprintf(fileName, "Eigentensor%05d.vtk", commitStep );
  theFile = fopen(fileName, "w");
  if (!theFile) {
    printf( "DomainViewer::saveEigenTensor, unable to open file for writing\n");
    return;
  }


  currElement = 1;


  Vector bounds = theDomain->getPhysicalBounds();
  
  deltaX = fabs(bounds(3) - bounds(0));
  deltaY = fabs(bounds(4) - bounds(1));
  deltaZ = fabs(bounds(5) - bounds(2));

  incr =  minDistBetweenPts()/2.0;
  numX = deltaX/incr;
  numY = deltaY/incr;
  numZ = deltaZ/incr;

  if ( ((numX+1)*(numY+1)*(numZ+1))  < desiredPoints ) {
 
    a = desiredPoints/( deltaX*deltaY*deltaZ ); // a^3 = desiredpts/x*y*z
 
    a = cubeRoot(a); 

    numX = a*deltaX;
    numY = a*deltaY;
    numZ = a*deltaZ;

    incr = deltaX/numX;
   printf("Interpolation scheme 1: number of points %d\n", numX*numY*numZ);
  }
  else
  printf("Interpolation scheme 2: number of points %d\n", numX*numY*numZ);



 

  fprintf(theFile, "# vtk DataFile Version 3.0\n");
  fprintf( theFile,
	  "3x3 Stress Tensors on regular grid, timestep %d\n",
	   commitStep);
  fprintf( theFile,"ASCII\n");
  fprintf( theFile,"DATASET STRUCTURED_POINTS\n");
  fprintf( theFile, "DIMENSIONS %d %d %d\n",numX+1, numY+1,numZ+1);
  fprintf( theFile,"ORIGIN %f %f %f\n", bounds(0), bounds(1), bounds(2) );
  fprintf( theFile, "SPACING %f %f %f\n", incr, incr, incr);
  fprintf( theFile, "POINT_DATA %d\n", (numX+1)* (numY+1)*(numZ+1));

  //  fprintf("TENSORS fourthOrder float\n");
  // fprintf( theFile, "SCALARS stress_lode_angle float 1\n");
  //  fprintf( theFile, "LOOKUP_TABLE default\n"); 

  fprintf( theFile, "TENSORS eigenmode float\n");
 
  for ( k = 0; k <= numZ; k++ ) {
    for ( j = 0; j <= numY; j++ ) {
      for ( i = 0; i <= numX; i++) {

	zero = false;
	loc[0] = bounds(0) + (i*incr);
	loc[1] = bounds(1) + (j*incr);
	loc[2] = bounds(2) + (k*incr);

	
	// check for inside current element & interpolate
	if (  currElement != -1  &&
	      stfdraw[currElement].insideElement( loc ) ) {
	    stfdraw[currElement].interpolateEigenTensor( m_eigen, loc );
	}
	
        else {  // find element & interpolate
	   currElement = findElement( loc );
	  if ( currElement != -1  ) {
	    stfdraw[currElement].interpolateEigenTensor( m_eigen, loc );
	  }

      
	  else { // not found, zero matrix
	    zero = true;

	  } // no such location 
	} // find element & interpolate
	
       	// if zero, print zero eigentensor
	if ( zero ) {
	  for (m = 0; m < 3; m++) {
	    for (n = 0; n < 3; n++ ) {
	      fprintf(theFile, "%.16g ", 0.00 );
	    }
	    fprintf(theFile,"\n");
	  }
	  fprintf(theFile,"\n");
	}
	else {
	  for (m = 0; m < 3; m++) {
	    for (n = 0; n < 3; n++ ) { 
	      fprintf(theFile, "%.16g ", m_eigen(m,n ));
	    }
	    fprintf(theFile,"\n");
	  }
	  fprintf(theFile,"\n");
	}

      } //k
    
    }//j
  
  } // i
  
  fclose(theFile);

  
}



/**
 * Brute force.
 * @returns the element number or -1 if not contained
 * in any brick
*/
int DomainViewer::findElement(double loc[3] ) {
  int tag = -1;
    // finally brute force
  ElementIter &theElements = theDomain->getElements();
  Element * theElement;

  while (( theElement = theElements() ) != 0) {
     if (theElement->getClassTag() == ELE_TAG_EightNodeBrick &&  
	 stfdraw[theElement->getTag()].insideElement(loc) )
       tag = theElement->getTag();
  }
   
  // otherwise, bail out
  return tag;
}


/**
 * Ask each brick the min distance between gauss pts
 */
double DomainViewer::minDistBetweenPts() {
  double best = 1.0e16;
  double curr;
  ElementIter &theElements = theDomain->getElements();
  Element * theElement;

  while (( theElement = theElements() ) != 0) {
    if (theElement->getClassTag() == ELE_TAG_EightNodeBrick) {
      curr = stfdraw[theElement->getTag()].getMinDistanceBetweenPts();
      if ( curr > 0.0 && curr < best )
	best = curr;
    }
  }
  return best;
}


/**
 * Halley's Rational Formula for cube roots 
 * citation: G. Alefeld, "On the convergence of Halley's method", 
 * Amer. Math. Monthly, 88 (1981) 530--536. 
 * see http://www.mathpath.org/Algor/cuberoot/algor.cube.root.halley.htm
 */
double DomainViewer::cubeRoot(double val) {
  int i,j;
  int notFound = 1;
  double guess, next, cube;
  
  if ( val < 1 ) {
    printf("Error, DomainViewer::cubeRoot, val must be >= 1, val = %f\n", val);
    return 0;
  }
  
  i = 1;
  

  // find guess
  while( notFound ) {

    guess = i*i*i;
    if ( guess >= val ) {

      if (guess == val) {
	return i;
      }

      else {  // we have the range, refine with Halley

	guess = (i-1) + 0.5; // choose midpoint between curr and previous
 
	for (j = 0; j < 3; j++) {

	  cube = guess*guess*guess;
 
	  next = guess*( cube + (2*val) );

	  next = next/( 2*cube + val );

	  guess = next;
	}
	return guess;
      }
    } // i^3 >= val

    else 
      i++;

  }
  
  
  
}

int DomainViewer::Matrix2TensorSysR4(const Matrix& M, double T[4][4][4][4])

{
 

  int nr = M.noRows();
  int nc = M.noCols();
  if (nr < 6 || nc < 6) {
    opserr << "NewTemplate3Dep::Matrix2TensorSysR4 - matrix must be no less than (6, 6)" << endln;
    return 1;
  }

  double sqrthalf = sqrt(0.5);
  double half = 0.5;
  
  // Adopt method from Helnwein (2001):  
  
  T[1][1][1][1] = M(0,0);
  T[1][1][2][2] = M(0,1);
  T[1][1][3][3] = M(0,2);      
  T[1][1][1][2] = T[1][1][2][1] = M(0,3) *sqrthalf;
  T[1][1][2][3] = T[1][1][3][2] = M(0,4) *sqrthalf;
  T[1][1][1][3] = T[1][1][3][1] = M(0,5) *sqrthalf;      

  T[2][2][1][1] = M(1,0);
  T[2][2][2][2] = M(1,1);
  T[2][2][3][3] = M(1,2);      
  T[2][2][1][2] = T[2][2][2][1] = M(1,3) *sqrthalf;
  T[2][2][2][3] = T[2][2][3][2] = M(1,4) *sqrthalf;
  T[2][2][1][3] = T[2][2][3][1] = M(1,5) *sqrthalf; 

  T[3][3][1][1] = M(2,0);
  T[3][3][2][2] = M(2,1);
  T[3][3][3][3] = M(2,2);      
  T[3][3][1][2] = T[3][3][2][1] = M(2,3) *sqrthalf;
  T[3][3][2][3] = T[3][3][3][2] = M(2,4) *sqrthalf;
  T[3][3][1][3] = T[3][3][3][1] = M(2,5) *sqrthalf;

  T[1][2][1][1] = T[2][1][1][1] = M(3,0) *sqrthalf;
  T[1][2][2][2] = T[2][1][2][2] = M(3,1) *sqrthalf;
  T[1][2][3][3] = T[2][1][3][3] = M(3,2) *sqrthalf;      
  T[1][2][1][2] = T[2][1][1][2] = T[1][2][2][1] = T[2][1][2][1] = M(3,3) *half;
  T[1][2][2][3] = T[2][1][2][3] = T[1][2][3][2] = T[2][1][3][2] = M(3,4) *half;
  T[1][2][1][3] = T[2][1][1][3] = T[1][2][3][1] = T[2][1][3][1] = M(3,5) *half;

  T[2][3][1][1] = T[3][2][1][1] = M(4,0) *sqrthalf;
  T[2][3][2][2] = T[3][2][2][2] = M(4,1) *sqrthalf;
  T[2][3][3][3] = T[3][2][3][3] = M(4,2) *sqrthalf;      
  T[2][3][1][2] = T[3][2][1][2] = T[2][3][2][1] = T[3][2][2][1] = M(4,3) *half;
  T[2][3][2][3] = T[3][2][2][3] = T[2][3][3][2] = T[3][2][3][2] = M(4,4) *half;
  T[2][3][1][3] = T[3][2][1][3] = T[2][3][3][1] = T[3][2][3][1] = M(4,5) *half;

  T[1][3][1][1] = T[3][1][1][1] = M(5,0) *sqrthalf;
  T[1][3][2][2] = T[3][1][2][2] = M(5,1) *sqrthalf;
  T[1][3][3][3] = T[3][1][3][3] = M(5,2) *sqrthalf;      
  T[1][3][1][2] = T[3][1][1][2] = T[1][3][2][1] = T[3][1][2][1] = M(5,3) *half;
  T[1][3][2][3] = T[3][1][2][3] = T[1][3][3][2] = T[3][1][3][2] = M(5,4) *half;
  T[1][3][1][3] = T[3][1][1][3] = T[1][3][3][1] = T[3][1][3][1] = M(5,5) *half;

  return 0;
}
