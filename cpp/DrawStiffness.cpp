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
** DrawStiffness.cpp is a class to draw a finite element's
** integration points. The data ( points and sphere size ) can be kept
** either inside or outside the object.
**/

#include "DrawStiffness.h"
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

#ifndef WIN32
#include <sys/time.h>
#include <sys/resource.h>
#endif


DrawStiffness::DrawStiffness() {

  sphere = gluNewQuadric();
  el = NULL;
  drawByMatrixAlone = false;
  sphereSize = 0;
  numPts = 0;
  selectable = true; 
  minDistBetweenPts = -1;

  printf("*");
  try {
    rg = new ReynoldsGlyph[8];

  }
  catch( std::bad_alloc & oops) {
    printf("DrawStiffness, Unable to allocate reynolds glyphs %s\n",
   oops.what());


    exit(0);
  }

}
 


DrawStiffness::~DrawStiffness() {
  delete []rg;
  gluDeleteQuadric( sphere );
}


/**
 *@pre: element has been set */
void DrawStiffness::setToStress() {
  int i;
  if ( el != NULL &&  el->getNumExternalNodes() == 8  ) {
    for ( i = 0; i < 8; i++)
      rg[i].setDrawMethod(STRESS);
  }
}


void DrawStiffness::setToStiffness() {

  int i;
  if ( el != NULL &&  el->getNumExternalNodes() == 8  ) {
     for ( i = 0; i < 8; i++)
      rg[i].setDrawMethod(PETAL0R);
  }
  
}

double DrawStiffness::getMinDistanceBetweenPts() {
  return  minDistBetweenPts;
}

/**
 * For time efficiency, save gauss point locations and sphere size.
 * draw many times after setting.
 * @return 0 for error, non-zero for success
 */
int DrawStiffness::setElement( Element * ele ) {
  int i, success;
  double curr = 1000000;
  GaussCoord3D gcd; 
  Vector v;//good
  Information info(v);//good
  int changedStiffness = 0;
  // DispBeamColumn3d * beam;

  if ( true ) {
	el = ele;
	numPts = 0; 

	if ( el->getClassTag() == ELE_TAG_EightNodeBrick ) {
	  numPts = el->getNumExternalNodes();
	  for( i = 0; i < el->getNumExternalNodes(); i++) {
	   
	    rg[i].setElement(ele);
	    rg[i].setIntegrPoint(i);
	    changedStiffness += rg[i].getChangedStiffness();
	  }
	  if (changedStiffness > 0)
	    selectable = true;
	  else selectable = false;
	 
	  setFacetSign(); 
	  

	  success = gcd.computeGaussPoint8Node( ele->getNodePtrs(), info );
  
	  if ( success != 0 ) {
		  printf( "DrawStiffness::computeGaussPoint8Node didn't work\n");
		  el = NULL;
		  return 0; //error
	  }
	  const Vector & gaussData  = info.getData();
 
	  if ( &gaussData == NULL ) {
	    printf( "DrawStiffness::No gauss data found in info Vector\n");
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

          // test
	  /*
	  double pt[3];
	  for (i = 0; i < numPts; i++ ) {
	    
	    pt[0] =  points(i,0);
	    pt[1] = points(i,1);
	    pt[2] = points(i,2);
   
	   
	    insideElement( pt );

	  }
	  // another test
	
	  pt[0] = 12;
	  pt[1] = 12;
	  pt[2] = 12;

	  insideElement( pt );
	  fflush(stdout);
	  */
	  sphereSize = 10000000; 
	  for ( i = 1; i < numPts; i++ ) {
	    
	    curr = Mdistance( points, 0, i );
	    if ( curr < sphereSize )
	      sphereSize = curr;
	  }
	  sphereSize = sphereSize;
	  minDistBetweenPts = curr;
	 

       

	} // 8 nodes

	else if ( el->getNumExternalNodes() == 2 )
	  ;/*{
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
	  minDistBetweenPts = curr;
	  
	  }*/ // 2 nodes
	drawByMatrixAlone = false;
	return 1;
  }//canDraw
  return 0; // can't draw

} 


/**
 * pre! data has the correct number of gauss points 
 * @returns true if stress has occured
*/
bool DrawStiffness::setMeanStress( double * data ) {
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


int DrawStiffness::getNumPts() {
  return numPts;
}

const Matrix &  DrawStiffness::getMatrix() {
  return points;
}

double  DrawStiffness::getSphereSize() {
  return sphereSize;
}

/**
 * DomainViewer may wish to keep Matrix (calculate once
 * and save data). So here we simply copy the points rather 
 * than the element
 */
void DrawStiffness::setPoints( const Matrix & other, double radius ) {

  points = other;

  numPts = points.noRows();

  sphereSize = radius;
 
  drawByMatrixAlone = true;
  
}

/**
 * @pre element has been set or sphere size and point set have
 * been set
 */
void DrawStiffness::drawElement(Element * ele, double * colors, 
			    ColorRange range, int colorStyleFlag )
{
  int i;
   
  if (el != NULL ) {
     
    //	set_color_by_param(0.0, 1.0);
    
    for ( i = 0; i < numPts; i++ ) {  
    
       if ( colorStyleFlag == BY_STIFF_EIGENMODE ) {

	glPushMatrix();
	glTranslatef( points(i,0), points(i,1), points(i,2) );	

	//glScalef(0.1,0.1,0.1); // original size is diameter of 2
	glScalef(sphereSize/2, sphereSize/2, sphereSize/2);
		rg[i].drawEigenMode();
	glPopMatrix();
      }
      // for stiffness or eigenmode, use stiffness filter
      
      else if ( range.inRange( colors[i] ) ) {
      
	glPushMatrix();
	glTranslatef( points(i,0), points(i,1), points(i,2) );	

	//glScalef(0.1,0.1,0.1); // original size is diameter of 2
	glScalef(sphereSize/2, sphereSize/2, sphereSize/2);


	  
	 if (range.getDataName() == BY_STIFFNESS || 
	    range.getDataName() ==  BY_MEAN_STRESS ||
	    range.getDataName() ==  BY_DEVIATORIC ) {
	  rg[i].draw(colors[i]);
	 
	}

	glPopMatrix();
      }
	        
    }// for each gauss pt
    fflush(stdout);
  }//non-null element

}

  


/** return if this method is capable of drawing 
 * this element type */
bool  DrawStiffness::canDraw( Element * el )
{ 
  /*
  if ( drawByMatrixAlone )
    return true;

  return (el->getNumExternalNodes() == 8  ||
  el->getClassTag() ==  ELE_TAG_DispBeamColumn3d ); */
  return selectable;
}

/**
 * Find out what DrawMethod this is */
int DrawStiffness::getClassTag(){
  return DRAW_STIFFNESS;
}


/**
 * Distance between 2 gauss points( used for deciding sphere size)
 */
double  DrawStiffness::Mdistance( Matrix pts, int j, int k ) {
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

/**
 * Calculate which normal direction is inside element
 * for a fixed order of traversing a facet's nodes
 */
void DrawStiffness::setFacetSign() {
  int i,j;//, k;
  Vector v; // OpenSees vector
  double pts[3][3];   // vector to array transform 
  double crossPt[3];

  Node ** theNodes; // nodes for a single element

  theNodes = el->getNodePtrs();

  //----  4 sides of parallelpiped ---
  for ( i = 0; i < 4; i++ ) {
    //  product = 1.0;
 
    facetSign[i] = 0.0;
    for ( j = 0; j < 3; j++ ) { // single facet
      
      v = theNodes[brick8Loops[i][j]]->getCrds();
      pts[j][X] = v[X];
      pts[j][Y] = v[Y];
      pts[j][Z] = v[Z];
    }
    
    // now, pick a point from the opposite face to test
    v =  theNodes[oppose8Loops[i][0]]->getCrds();
    crossPt[X] = v[X];
    crossPt[Y] = v[Y];
    crossPt[Z] = v[Z];
    
    // now, test against the opposing face
    facetSign[i] =  orient3D( pts[0], pts[1], pts[2], crossPt );
 
  }

  //--- top quad. ----------
    for ( j = 0; j < 3; j++) {
    v = theNodes[j]->getCrds();
    pts[j][X] = v[X];
    pts[j][Y] = v[Y];
    pts[j][Z] = v[Z];
  }

  v =  theNodes[6]->getCrds();
  crossPt[X] = v[X];
  crossPt[Y] = v[Y];
  crossPt[Z] = v[Z];
  
  facetSign[4] = orient3D( pts[0], 
			  pts[1], pts[2], crossPt );


   //---- bottom quad ---
  for ( j = 4; j < 7; j++) {
    v = theNodes[j]->getCrds();
    pts[j-4][X] = v[X];
    pts[j-4][Y] = v[Y];
    pts[j-4][Z] = v[Z];
  }
  
  v =  theNodes[0]->getCrds();
  
  crossPt[X] = v[X];
  crossPt[Y] = v[Y];
  crossPt[Z] = v[Z];
  
  facetSign[5] = orient3D( pts[0], 
			   pts[1], pts[2], crossPt );

  /*  printf("orient3D returns %f, positive for below\n", orient3D( pts[0], 
			  pts[1], pts[2], crossPt) );


  printf("============facet signs (should be negative) ==========\n");
  for (i = 0; i < 6; i++ )
    printf("%f ", facetSign[i]);
  printf("\n========================\n");
  fflush(stdout);*/
}


/**
 * Check whether the point is inside the 8 node element
 * inside includes on a facet. 
 * We depend on the point order as given in OpenSees. 
*/
bool DrawStiffness::insideElement( double thePoint[3] ) {
  int i,j;
  Vector v; // OpenSees vector
  double pts[3][3];   /// vector to array transform 
  double sign;
  Node ** theNodes; // nodes for a single element

  if (el ->getClassTag() != ELE_TAG_EightNodeBrick )
    return false;

  theNodes = el->getNodePtrs();
  //  printf("=========thePoint!!!: %f %f %f --------------\n",
  // thePoint[0], thePoint[1], thePoint[2]);
 

  //----  4 sides of parallelpiped ---
  for ( i = 0; i < 4; i++ ) {
 
    for ( j = 0; j < 3; j++ ) { // single facet
      v = theNodes[brick8Loops[i][j]]->getCrds();
      pts[j][X] = v[X];
      pts[j][Y] = v[Y];
      pts[j][Z] = v[Z];
    }

    //check facet sign, test against oriented points
    if ( facetSign[i] >= 0 )
      sign = orient3D(pts[0], pts[1], pts[2], thePoint );
    else 
       sign = orient3D(pts[2], pts[1], pts[0], thePoint );
    
    if (sign >= 0)
      ;
      //	printf("insideElement::sign %d is positive\n", i);
    else {
      //   printf("insideElement: sign %d negative, outside %f %f %f\n", i,
      //      thePoint[0], thePoint[1], thePoint[2]);
      return false;

    }
  }// 4 sides of parallelpiped

  //--- top quad. ----

  for ( j = 0; j < 3; j++) {
    v = theNodes[j]->getCrds();
    pts[j][X] = v[X];
    pts[j][Y] = v[Y];
    pts[j][Z] = v[Z];
  }

 //check facet sign, test against oriented points
    if ( facetSign[4] >= 0 )
      sign = orient3D(pts[0], pts[1], pts[2], thePoint );
    else 
       sign = orient3D(pts[2], pts[1], pts[0], thePoint );
    
    if (sign >= 0)
      ;
    //	printf("insideElement::sign %d is positive\n", 4);
    else {
      //  printf("insideElement: sign %d negative, outside  %f %f %f\n", 4,
      //     thePoint[0], thePoint[1], thePoint[2]);
      return false;

    }

 
   //---- bottom quad ---
   for ( j = 4; j < 7; j++) {
     v = theNodes[j]->getCrds();
     pts[j-4][X] = v[X];
     pts[j-4][Y] = v[Y];
     pts[j-4][Z] = v[Z];
   }

 //check facet sign, test against oriented points
    if ( facetSign[5] >= 0 )
      sign = orient3D(pts[0], pts[1], pts[2], thePoint );
    else 
       sign = orient3D(pts[2], pts[1], pts[0], thePoint );
    
    if (sign >= 0)
      ; //printf("insideElement::sign %d is positive\n", 5);
    else {
      //  printf("insideElement: sign %d negative, outside %f %f %f\n", 5,
      //    thePoint[0], thePoint[1], thePoint[2]);
      return false;

    } 

    return true;
}



/**
 * @returns positive value if query point, d, lies "below"
 * oriented plane passing through a,b,c, where a,b,c
 * appear in counterclockwise order when viewed from "above"
 *
 * Source: Jonathan Richard Shewchuk, Adaptive Precision
 *  Floating-Point Arithmetic and Fast Robust Geometric Predicates, 
 * Discrete & Computational Geometry 18:305-363, 1997
 */ 
double  DrawStiffness::orient3D( double a[3], double b[3], double c[3], double d[3] )

{
  double result = 0.0;
  double mx[3][3];
  int i,j;

  // build the matrix
  for( i = 0; i < 3; i++) {
    mx[0][i] = a[i] - d[i]; 
    mx[1][i] = b[i] - d[i];
    mx[2][i] = c[i] - d[i];
  }

  // calculate the determinant
 
  result = result + mx[0][0]*( (mx[1][1]*mx[2][2]) - (mx[1][2]*mx[2][1]) );

  result = result - mx[0][1]*( (mx[1][0]*mx[2][2]) - (mx[1][2]*mx[2][0]) );

  result = result + mx[0][2]*( (mx[1][0]*mx[2][1]) - (mx[1][1]*mx[2][0]) );

  return result;
  
} 


/**
 * go through gauss points, print
 * File format: 
 * current elastic stiffness
 * current stiffness 
 * yield surface normal 
 * plastic flow normal
 * deformation
 */
void DrawStiffness::printStiffness( FILE * theFile ) {
  int i,j, k, node;
  Matrix m_elastic(6,6);
  Matrix m_curr(6,6);
  Matrix m_yieldNorm(3,3);
  Matrix m_flowNorm(3,3);
  Node ** theNodes = el->getNodePtrs();
  int indeterminate,sameStiff;  

  for (k = 0; k < 8; k++ ) {   // for each gauss point
  
    rg[k].getStiffData(m_elastic,  m_curr, m_yieldNorm, m_flowNorm,
		       indeterminate,sameStiff);

    //original elastic 6x6
    for ( i = 0; i < 6; i++ ) {
      for (j = 0; j < 6; j++ ) {
	fprintf(theFile, "%.16g ", m_elastic(i,j));
      }
      fprintf(theFile, "\n");
    }
    fprintf(theFile, "\n\n");
    
    // current elastic-plastic
    for ( i = 0; i < 6; i++ ) {
      for (j = 0; j < 6; j++ ) {
	fprintf(theFile, "%.16g ", m_curr(i,j));
      }
      fprintf(theFile, "\n");
    }
    fprintf(theFile, "\n\n"); 

   
    // yield surface normal
    for ( i = 0; i < 3; i++ ) {
      for (j = 0; j < 3; j++ ) {
	fprintf(theFile, "%.16g ", m_yieldNorm(i,j));
      }
      fprintf(theFile, "\n");
    }
    fprintf(theFile, "\n\n"); 

    // plastic flow normal
    for ( i = 0; i < 3; i++ ) {
      for (j = 0; j < 3; j++ ) {
	fprintf(theFile, "%.16g ", m_flowNorm(i,j));
      }
      fprintf(theFile, "\n");
    }
    fprintf(theFile, "\n\n");  
   
    
    node = node2Gauss[k]; // closest node index to gauss location
    Vector  v =  theNodes[node]->getDisp();
    for ( i = 0; i < v.Size(); i++ )
      fprintf(theFile, "%.16g ",v(i));
    fprintf(theFile, "\n\n");  

    

  }

}


void DrawStiffness::printRotationTest(FILE * theFile) {
  int i,j,k,m,x,indeterminate,sameStiff;
  Matrix m_rotate(6,6);
  int dim[4] = {3,3,3,3};

  Tensor t_elastic( 4,dim, 0.0 );
  Tensor t_rotate( 4,dim, 0.0 );

  Matrix m_elastic(6,6);
  Matrix m_curr(6,6);
  Matrix m_EeMinusEep(6,6);
  Matrix m_yieldNorm(3,3);
  Matrix m_flowNorm(3,3);

  if (el->getClassTag()  ==  ELE_TAG_EightNodeBrick ) {
    fprintf(theFile, "Element %d\n", el->getTag());

    for (x = 0; x < 8; x++ ) {

      fprintf(theFile,"gausspt %d\n", x);  
      fprintf(theFile, "%.16g %.16g %.16g\n", 
	      points(gauss2Node[x],0), points(gauss2Node[x],1), 
	      points(gauss2Node[x],2));

      rg[x].getStiffData(m_elastic, m_curr, m_yieldNorm, m_flowNorm, 
			 indeterminate,sameStiff);

      fprintf(theFile,"ys_normal\n");
      // yield surface normal
      for ( i = 0; i < 3; i++ ) {
	for (j = 0; j < 3; j++ ) {
	  fprintf(theFile, "%.16g ", m_yieldNorm(i,j));
	}
	fprintf(theFile, "\n");
      }
      fprintf(theFile, "\n\n"); 
      
      fprintf(theFile,"pf_normal\n");
      // plastic flow normal
      for ( i = 0; i < 3; i++ ) {
	for (j = 0; j < 3; j++ ) {
	  fprintf(theFile, "%.16g ", m_flowNorm(i,j));
	}
	fprintf(theFile, "\n");
      }
      fprintf(theFile, "\n\n"); 

      fprintf(theFile,"m_elastic_curr\n");
      // current elastic
      for ( i = 0; i < 6; i++ ) {
	for (j = 0; j < 6; j++ ) {
	  fprintf(theFile, "%.16g ", m_curr(i,j));
	}
	fprintf(theFile, "\n");
      }
      fprintf(theFile, "\n\n"); 

      rg[x].getRotationMatrix(m_rotate);
      
      fprintf(theFile,"m_rotate\n");
      for ( i = 0; i < 6; i++ ) {
	for (j = 0; j < 6; j++ ) {
	  fprintf(theFile, "%.16g ", m_rotate(i,j));
	}
	fprintf(theFile, "\n");
      }
	 
      fprintf(theFile, "\n\n"); 

      
      
      fprintf(theFile, "indeterminate %d\n\n", indeterminate); 
      fprintf(theFile, "sameStiff %d\n\n", sameStiff);

      // C-T rotation, run vangeldarpolar( C-T,R,U,toler)
      /*   m_EeMinusEep = m_elastic - m_curr;
      fprintf(theFile,"m__EeMinusEep\n");
      if( rg[x].vanGelderPolar(  m_EeMinusEep, m_rotate, m_elastic, 1.0e12) != 0 )
	{
	  for ( i = 0; i < 6; i++ ) {
	    for (j = 0; j < 6; j++ ) {
	      fprintf(theFile, "%.16g ", 0.0);
	    }
	    fprintf(theFile, "\n");
	  }
	  
	  fprintf(theFile, "\n\n"); 
	}
      else {
	for ( i = 0; i < 6; i++ ) {
	  for (j = 0; j < 6; j++ ) {
	    fprintf(theFile, "%.16g ", m_rotate(i,j));
	  }
	  fprintf(theFile, "\n");
	}
	
	fprintf(theFile, "\n\n"); 
      }
      */
    }// for each gauss point
  }// if 8 node brick
}




/** c = A:B/ |A| |B| 
 *@pre: A, B same dimensions
*/
double DrawStiffness::matrixDifference(Matrix m_A, Matrix m_B) {
  double magA,magB,AcolonB;
  int i,j;
  magA = magB = AcolonB = 0;
    
  if (m_A.noRows() != m_B.noRows() ||
      m_A.noCols() != m_B.noCols() )
    return -1;

    for (i = 0; i < m_A.noRows(); i++ ) {
      for (j = 0; j < m_A.noCols(); j++ ) {
	printf("m");
	magA = magA + (m_A(i,j)*m_A(i,j));
	magB = magB + (m_B(i,j)*m_B(i,j));	
	AcolonB = AcolonB + (m_A(i,j)*m_B(i,j));
      }
    }

    magA = sqrt(magA);
    magB = sqrt(magB);

    AcolonB =  AcolonB/(magA*magB);
    if (AcolonB < .0000000000001)
      printf("too tiny!!\n");

    return AcolonB;
}

void DrawStiffness::normalizeMatrix( Matrix & theMatrix ) {
  int i, j;
  double norm;
   norm = 0;
    
    for (i = 0; i < theMatrix.noRows(); i++ ) {
      for (j = 0; j < theMatrix.noCols(); j++ ) {
	norm = norm + ( theMatrix(i,j)* theMatrix(i,j) );
      }
    }
    norm = sqrt(norm);

    for (i = 0; i < 6; i++ ) {
      for (j = 0; j < 6; j++ ) {
	theMatrix(i,j) =  theMatrix(i,j)/norm ;
      }
    }

}


// gauss points
void DrawStiffness::printPointsVTK( FILE * theFile ) {
  int i;
  if (numPts == 8) {
    for (i = 0; i < numPts; i++ ) {
      fprintf(theFile, "%.16g %.16g %.16g\n", 
	      points(gauss2Node[i],0), points(gauss2Node[i],1), 
	      points(gauss2Node[i],2));
      
    }
  }
}

void DrawStiffness::printMinStiffVTK( FILE * theFile ) {
  int i;
  // ask each reynolds glyph for its min stiffness 
  if (numPts == 8) {
    for (i = 0; i < numPts; i++ ) {
      fprintf(theFile, "%.16g ",  rg[i].getLowestEigen());
    }
  }
}

/*
 * red: isotropic part ( color[0] )
 * green: shear part ( color[1] )
 * blue: axisymmetric part ( color[2] ) 
 * convert into scalar between 0 and 1
 */
void DrawStiffness::printEigenmodeVTK( FILE * theFile ) {
  int i;
 float color[4];
 float value;
 float cubed = 16777216; // 256*256*256
  // ask each reynolds glyph for its eigenmode
  if (numPts == 8) {
    for (i = 0; i < numPts; i++ ) {
      
      rg[gauss2Node[i]].eigenMode(color);
      
      // red = red*255; green = green*255; blue = blue*255 
      color[0] = color[0]*255;
      color[1] = color[1]*255;
      color[2] = color[2]*255;
      
      //int.value=256*256*red+256*green+blue     
      value = (256*256*color[0]) + (256*color[1]) + color[2];

      //float.value = int.value/256*256*256    
      value = value/cubed; // a number between 0 and .996

      fprintf(theFile, "%.16g ", value );
    }
  }
}

/**
 * precondition: state has been updated */
void DrawStiffness::setMinStiffness( double * data ) {
  int i;
 
  if (numPts > 0 ) { 
    for( i = 0; i < numPts; i++ ) {
      data[i] = rg[i].getLowestEigen();
    }
  }
}



/**
 * If element is not an eight node brick, zero matrix will be returned.
 * Otherwise, the result is Shepard's interpolation 
 * of gauss point stiffness within element
 * m_stiff is 6x6 matrix 
 */
void DrawStiffness::interpolateStiffness(Matrix & m_stiff, double loc[3]) {
  double sumOneODistSquared = 0; // SUM (distance^2)
  double distSquared[8];
  int i, j, k;
  double dist;
  int numGauss = 8;

  for ( i = 0; i < 6; i++ ) {
    for ( j = 0; j < 6; j++ ) {
      m_stiff(i,j) = 0;
    }
  }

  if (el->getClassTag()  !=  ELE_TAG_EightNodeBrick ) {
     return;
   }

  Matrix m_zero(6,6);
  Matrix m_one(6,6);
  Matrix m_two(6,6);
  Matrix m_three(6,6);
  Matrix m_four(6,6);
  Matrix m_five(6,6);
  Matrix m_six(6,6);
  Matrix m_seven(6,6);

  Matrix m_min(6,6);
  Matrix m_max(6,6);

  Matrix * matrices[8];
  matrices[0] = &m_zero;
  matrices[1] = &m_one;
  matrices[2] = &m_two;
  matrices[3] = &m_three;
  matrices[4] = &m_four;
  matrices[5] = &m_five;
  matrices[6] = &m_six;
  matrices[7] = &m_seven;
  
   //-------- denominator --------
  for ( i = 0; i < numGauss; i++ )
    distSquared[i] = 0;

  // --- distance to each gauss, squared -----
  for ( i = 0; i < numGauss; i++ ) {
    for (j = 0; j < 3; j++) { // x,y,z
  
      dist = loc[j] - points(i,j);
    
      distSquared[i] =  distSquared[i] + (dist*dist);
   }
  
  }


  // sum of one over squared distances
  sumOneODistSquared = 0;
  for ( i = 0; i < numGauss; i++ ) {
    sumOneODistSquared = sumOneODistSquared + 1.0/distSquared[i];
  }

  // get the data
  for (i = 0; i < numGauss; i++)
    rg[i].getStiffness( *matrices[i] );
  
  //---- sanity test --------------
   for (i = 0; i < 6; i++) {
     for ( j = 0; j < 6; j++ ) {
       m_min(i,j) = m_max(i,j) = (*matrices[0])(i,j);
      
     }
   }

   for( k = 1; k < numGauss; k++ ) {
     for (i = 0; i < 6; i++) {
       for ( j = 0; j < 6; j++ ) {
	 if ( (*matrices[k])(i,j) < m_min(i,j) )
	   m_min(i,j) = (*matrices[k])(i,j);
	 if ( (*matrices[k])(i,j) > m_max(i,j) )
	   m_max(i,j) = (*matrices[k])(i,j);
       }
     }

   }
   //--------------------------

  // ok, now for the Shepard's calculation, point by point
   for (i = 0; i < 6; i++) {
     for ( j = 0; j < 6; j++ ) {

      for( k = 0; k < numGauss; k++ ) {

	m_stiff(i,j) = m_stiff(i,j) + ( (*matrices[k])(i,j)/distSquared[k] );

      }

      // now handle the denominator
      m_stiff(i,j) = m_stiff(i,j) /sumOneODistSquared;

     }//j
   } // i
}

//----------------------------------------------------
void DrawStiffness::interpolateStress(Matrix & m_stress, double loc[3]) {
  double sumOneODistSquared = 0; // SUM (distance^2)
  double distSquared[8];
  int i, j, k;
  double dist;
  int numGauss = 8;

  for ( i = 0; i < 3; i++ ) {
    for ( j = 0; j < 3; j++ ) {
      m_stress(i,j) = 0;
    }
  }

  if (el->getClassTag()  !=  ELE_TAG_EightNodeBrick ) {
     return;
   }

  Matrix m_zero(3,3);
  Matrix m_one(3,3);
  Matrix m_two(3,3);
  Matrix m_three(3,3);
  Matrix m_four(3,3);
  Matrix m_five(3,3);
  Matrix m_six(3,3);
  Matrix m_seven(3,3);

  Matrix m_min(3,3);
  Matrix m_max(3,3);

  Matrix * matrices[8];
  matrices[0] = &m_zero;
  matrices[1] = &m_one;
  matrices[2] = &m_two;
  matrices[3] = &m_three;
  matrices[4] = &m_four;
  matrices[5] = &m_five;
  matrices[6] = &m_six;
  matrices[7] = &m_seven;
  
   //-------- denominator --------
  for ( i = 0; i < numGauss; i++ )
    distSquared[i] = 0;

  // --- distance to each gauss, squared -----
  for ( i = 0; i < numGauss; i++ ) {
    for (j = 0; j < 3; j++) { // x,y,z
  
      dist = loc[j] - points(i,j);
    
      distSquared[i] =  distSquared[i] + (dist*dist);
   }
  
  }


  // sum of one over squared distances
  sumOneODistSquared = 0;
  for ( i = 0; i < numGauss; i++ ) {
    sumOneODistSquared = sumOneODistSquared + 1.0/distSquared[i];
  }

  // get the data
  for (i = 0; i < numGauss; i++)
    rg[i].getStress( *matrices[i] );
  
  //---- sanity test --------------
   for (i = 0; i < 3; i++) {
     for ( j = 0; j < 3; j++ ) {
       m_min(i,j) = m_max(i,j) = (*matrices[0])(i,j);
      
     }
   }

   for( k = 1; k < numGauss; k++ ) {
     for (i = 0; i < 3; i++) {
       for ( j = 0; j < 3; j++ ) {
	 if ( (*matrices[k])(i,j) < m_min(i,j) )
	   m_min(i,j) = (*matrices[k])(i,j);
	 if ( (*matrices[k])(i,j) > m_max(i,j) )
	   m_max(i,j) = (*matrices[k])(i,j);
       }
     }

   }
   //--------------------------

  // ok, now for the Shepard's calculation, point by point
   for (i = 0; i < 3; i++) {
     for ( j = 0; j < 3; j++ ) {

      for( k = 0; k < numGauss; k++ ) {

	m_stress(i,j) = m_stress(i,j) + ( (*matrices[k])(i,j)/distSquared[k] );

      }

      // now handle the denominator
      m_stress(i,j) = m_stress(i,j) /sumOneODistSquared;

     }//j
   } // i
}


void DrawStiffness::interpolateEigenTensor(Matrix & m_eigen, double loc[3]) {
  double sumOneODistSquared = 0; // SUM (distance^2)
  double distSquared[8];
  int i, j, k;
  double dist;
  int numGauss = 8;

  for ( i = 0; i < 3; i++ ) {
    for ( j = 0; j < 3; j++ ) {
      m_eigen(i,j) = 0;
    }
  }

  if (el->getClassTag()  !=  ELE_TAG_EightNodeBrick ) {
     return;
   }

  Matrix m_zero(3,3);
  Matrix m_one(3,3);
  Matrix m_two(3,3);
  Matrix m_three(3,3);
  Matrix m_four(3,3);
  Matrix m_five(3,3);
  Matrix m_six(3,3);
  Matrix m_seven(3,3);

  Matrix m_min(3,3);
  Matrix m_max(3,3);

  Matrix * matrices[8];
  matrices[0] = &m_zero;
  matrices[1] = &m_one;
  matrices[2] = &m_two;
  matrices[3] = &m_three;
  matrices[4] = &m_four;
  matrices[5] = &m_five;
  matrices[6] = &m_six;
  matrices[7] = &m_seven;
  
   //-------- denominator --------
  for ( i = 0; i < numGauss; i++ )
    distSquared[i] = 0;

  // --- distance to each gauss, squared -----
  for ( i = 0; i < numGauss; i++ ) {
    for (j = 0; j < 3; j++) { // x,y,z
  
      dist = loc[j] - points(i,j);
    
      distSquared[i] =  distSquared[i] + (dist*dist);
   }
  
  }


  // sum of one over squared distances
  sumOneODistSquared = 0;
  for ( i = 0; i < numGauss; i++ ) {
    sumOneODistSquared = sumOneODistSquared + 1.0/distSquared[i];
  }

  // get the data
  for (i = 0; i < numGauss; i++)
    rg[i].getEigenTensor( *matrices[i] );
  
  //---- sanity test --------------
   for (i = 0; i < 3; i++) {
     for ( j = 0; j < 3; j++ ) {
       m_min(i,j) = m_max(i,j) = (*matrices[0])(i,j);
      
     }
   }

   for( k = 1; k < numGauss; k++ ) {
     for (i = 0; i < 3; i++) {
       for ( j = 0; j < 3; j++ ) {
	 if ( (*matrices[k])(i,j) < m_min(i,j) )
	   m_min(i,j) = (*matrices[k])(i,j);
	 if ( (*matrices[k])(i,j) > m_max(i,j) )
	   m_max(i,j) = (*matrices[k])(i,j);
       }
     }

   }
   //--------------------------

  // ok, now for the Shepard's calculation, point by point
   for (i = 0; i < 3; i++) {
     for ( j = 0; j < 3; j++ ) {

      for( k = 0; k < numGauss; k++ ) {

	m_eigen(i,j) = m_eigen(i,j) + ( (*matrices[k])(i,j)/distSquared[k] );

      }

      // now handle the denominator
      m_eigen(i,j) = m_eigen(i,j) /sumOneODistSquared;

     }//j
   } // i
}

