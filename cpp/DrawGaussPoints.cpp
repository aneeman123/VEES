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
** DrawGaussPoints.cpp is a class to draw a finite element's
** integration points. The data ( points and sphere size ) can be kept
** either inside or outside the object.
**/
#include "DrawGaussPoints.h"
#include "GaussCoord3D.h"
#include <OPS_Stream.h>
#include <StandardStream.h>
#include <Recorder.h>
#include <ElementRecorder.h>
#include <Information.h>
#include <Response.h>

#include <DispBeamColumn3d.h>

//#include <EightNodeBrick.h>
#include <Matrix.h>

DrawGaussPoints::DrawGaussPoints() {
  sphere = gluNewQuadric();
  el = NULL;
  drawByMatrixAlone = false;
  sphereSize = 0;
  numPts = 0;
  selected = -1;
}
 

DrawGaussPoints::~DrawGaussPoints() {
  gluDeleteQuadric( sphere );
}


/**
 * highlight a selected point */
void DrawGaussPoints::setSelected(int gaussNum) {
  selected = gaussNum;
}

/**
 * For time efficiency, save gauss point locations and sphere size.
 * draw many times after setting.
 * @return 0 for error, non-zero for success
 */
int DrawGaussPoints::setElement( Element * ele ) {
  int i, success;
  double curr;
  GaussCoord3D gcd; 
  Vector v;//good
  Information info(v);//good

  DispBeamColumn3d * beam;

  if ( ele != el && canDraw( ele ) ) {
	el = ele;

	if ( el->getNumExternalNodes() == 8 ) {

	  success = gcd.computeGaussPoint8Node( ele->getNodePtrs(), info );
  
	  if ( success != 0 ) {
		  printf( "DrawGaussPoints::computeGaussPoint8Node didn't work\n");
		  el = NULL;
		  return 0; //error
	  }
	  const Vector & gaussData  = info.getData();
 
	  if ( &gaussData == NULL ) {
	    printf( "DrawGaussPoints::No gauss data found in info Vector\n");
	    el = NULL;
	    return 0; //error
	  }
	  
	  numPts = (int)gaussData(0);
	  points.resize( numPts, 3 );
	  // copy points to matrix
	  for (i = 0; i < numPts; i++ ) {
	    points(i,0) = gaussData[i*3+1];
	    points(i,1) = gaussData[i*3+2];
	    points(i,2) = gaussData[i*3+3];
	  }
	  sphereSize = 10000000; 
	  for ( i = 1; i < numPts; i++ ) {
	    
	    curr = Mdistance( points, 0, i );
	    if ( curr < sphereSize )
	      sphereSize = curr;
	  }
	  sphereSize = sphereSize/10;
	 
	  
	} // 8 nodes

	else if ( el->getNumExternalNodes() == 2 ) {
	  printf("DispBeamColumn3d!!!!!!!!!!!!!!!!\n");
	  fflush(stdout);
	  beam = (DispBeamColumn3d *) ( el );
	  numPts = beam->getNumGaussPoints();
	  points.resize( numPts, 3 );
	  printf("fetching beam gauss points\n");
	  beam->getGaussPoints( points );
	  
	  sphereSize = 10000000; 
	  for ( i = 1; i < numPts; i++ ) {
	    
	    curr = Mdistance( points, 0, i );
	    if ( curr < sphereSize )
	      sphereSize = curr;
	  }
	  sphereSize = sphereSize/2;
	  
	} // 2 nodes
	drawByMatrixAlone = false;
	return 1;
  }//canDraw
  return 0; // can't draw

} 


/**
 * pre! data has the correct number of gauss points 
 * @returns true if stress has occured
*/
bool DrawGaussPoints::setMeanStress( double * data ) {
  int i;
  Response * theResponse = NULL;
  Information eleInfo(1.0);
  bool stressed = false;

  char * argStress[] = { "stress" };
  
 
  if (el != NULL ) {
    //   printf("good mean stress!!!!!!!!!!\n");
    if ( el->getNumExternalNodes() == 8 ) {
    //try get stress
      theResponse = el->setResponse((const char **)argStress, 1, eleInfo); 
      if (  theResponse ) {
	theResponse->getResponse();

	Information &theInfo = theResponse->getInformation();
	const Vector &eleData = theInfo.getData();

	if ( eleData.Size() == 48 ) {
	  for( i = 0; i < 8; i++) {
            // 11 22 33 12 13 23
	    if( eleData(i*6) == 0.0 && eleData(i*6 +1) == 0.0 && 
		eleData(i*6 +2) == 0.0 )
	      
	      data[i] = 0.0;
	    
	    else {
	      data[i] = (eleData(i*6) + eleData(i*6 +1) + eleData(i*6 +2))/3.0;
	      stressed = true;
	    }
	      
	  }
	}

	else if ( eleData.Size() == 49 ) {
	  
		 
	  for( i = 0; i < 8; i++){
	    

	    if( eleData(i*6+1) == 0.0 && eleData(i*6 +2) == 0.0 && 
		eleData(i*6 +3) == 0.0 )
	      
	      data[i] = 0.0;
	    else {
	      data[i] = 
		(eleData(i*6+1) + eleData(i*6 + 2) + eleData(i*6 +3))/3.0;
	    
	      stressed = true;
	    }
	  }
	}
      } //non-null response
    }//8 nodes
    //    printf("data : %f %f %f \n", data[0],data[1],data[2]);
  }// non-null element

  return stressed;
}



/**
 * pre! data has the correct number of gauss points 
 * @returns true if stress has occured
*/
void DrawGaussPoints::setDeviatoricStress( double * data ) {
  int i;
  Response * theResponse = NULL;
  Information eleInfo(1.0);
  bool stressed = false;

  char * argStress[] = { "pqall" };
  
 
  if (el != NULL ) {
    //   printf("good mean stress!!!!!!!!!!\n");
    if ( el->getNumExternalNodes() == 8 ) {
    //try get stress
      theResponse = el->setResponse((const char **)argStress, 1, eleInfo); 
      if (  theResponse ) {
	theResponse->getResponse();

	Information &theInfo = theResponse->getInformation();
	const Vector &eleData = theInfo.getData();
       
	for( i = 0; i < 8; i++) {    
	  data[i] = eleData(i*2+1); // q for each gauss point
	}

       
      } //non-null response
    }//8 nodes
    //    printf("data : %f %f %f \n", data[0],data[1],data[2]);
  }// non-null element
}




int DrawGaussPoints::getNumPts() {
  return numPts;
}

const Matrix &  DrawGaussPoints::getMatrix() {
  return points;
}

double  DrawGaussPoints::getSphereSize() {
  return sphereSize;
}

/**
 * DomainViewer may wish to keep Matrix (calculate once
 * and save data). So here we simply copy the points rather 
 * than the element
 */
void DrawGaussPoints::setPoints( const Matrix & other, double radius ) {

  points = other;

  numPts = points.noRows();

  sphereSize = radius;
 
  drawByMatrixAlone = true;
  
}

/**
 * @pre element has been set or sphere size and point set have
 * been set
 */
void DrawGaussPoints::drawElement(Element * ele, double * colors, 
			    ColorRange range, int colorStyleFlag )
{
  int i, pick;
  /* if (ele != el || ( !drawByMatrixAlone ) ) { 
	if ( ! setElement(ele) ) {
	  return;
	}
    }
  */
  if (el != NULL ) {
    gluQuadricDrawStyle(sphere, GLU_FILL );
    gluQuadricNormals( sphere, GLU_SMOOTH );
    
    //	set_color_by_param(0.0, 1.0);
    
    for ( i = 0; i < numPts; i++ ) {
     
      // check each color against colorrange
      // before draw
    
      if (colorStyleFlag ==  SINGLE_COLOR ) 
      { // SINGLE COLOR

	glGetIntegerv(GL_RENDER_MODE,&pick); 
	if( pick == GL_SELECT )
	  glPushName( i ); //make gauss point pickable
	
	set_color_by_param( 0.0, 1.0);	  
	
	glPushMatrix();
	glTranslatef( points(i,0), points(i,1), points(i,2) );
	gluSphere( sphere, sphereSize, 10,10);
	glPopMatrix();
	
	if ( i == selected ) {
	  glPushMatrix();
	  glTranslatef( points(i,0), points(i,1), points(i,2) );
	  set_color_by_param( 1.0, 1.0);
	  glutWireSphere(sphereSize*1.2, 4, 4);
	  glPopMatrix();
	}

	glGetIntegerv(GL_RENDER_MODE,&pick); 
	if( pick == GL_SELECT )
	  glPopName();

     }
      
      else if (  range.inRange( colors[i] ) ) {

	  glGetIntegerv(GL_RENDER_MODE,&pick); 
	  if( pick == GL_SELECT )
	    glPushName( i ); //make gauss point pickable
	  
	  set_color_by_param(colors[i], 1.0);	  
	  
	  glPushMatrix();
	  glTranslatef( points(i,0), points(i,1), points(i,2) );
	  gluSphere( sphere, sphereSize, 10,10);
	  glPopMatrix();
	  
	  if ( i == selected ) {
	    glPushMatrix();
	    glTranslatef( points(i,0), points(i,1), points(i,2) );
	    
	    set_saturation_by_param( colors[i], 1.0); 
	    glutWireSphere(sphereSize*1.2, 4, 4);
	    glPopMatrix();
	  }
	  
	  glGetIntegerv(GL_RENDER_MODE,&pick); 
	  if( pick == GL_SELECT )
	    glPopName();
	  
	} // if in range
     
    }// for each gauss pt
  }//non-null element

}

  
/** return if this method is capable of drawing 
 * this element type */
bool  DrawGaussPoints::canDraw( Element * el )
{ 
 
  if ( drawByMatrixAlone )
    return true;

  return (el->getNumExternalNodes() == 8  ||
          el->getClassTag() ==  ELE_TAG_DispBeamColumn3d ); 
}

/**
 * Find out what DrawMethod this is */
int DrawGaussPoints::getClassTag(){
  return DRAW_GAUSS_PTS;
}


/**
 * Distance between 2 gauss points( used for deciding sphere size)
 */
double  DrawGaussPoints::Mdistance( Matrix pts, int j, int k ) {
  double delta[3];
  int i;
  double sum = 0;
  int dim = pts.noCols();

  for(i = 0; i < dim; i++ ) {

	delta[i] = pts(j,i) - pts(k,i);
	delta[i] = delta[i]*delta[i];
	sum += delta[i];
  }


  return (sqrt( sum )); 
}

