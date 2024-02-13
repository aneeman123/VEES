#include "GlyphViewer.h"
#include "auxiliary.h"
#include "materials.h"
#include "Camera.h"
#include <Brick.h>



GlyphViewer::GlyphViewer() {
  int i,j;

  saveImage = 0;
  theElement = NULL;
  intgrPoint = -1;
  eleNum = -1;

  drawMethod = REYNOLDS_GLYPH4;
  sizeRange.setData(ext, NUM_PTS);


  // draw the sphere
  for ( i = 0; i < NUM_PTS; i++ ) {
	ext[i] = 1.0;
	for ( j = 0; j < 3; j++ ) {
	  extPts[i][j] = ptLoc[i][j];
	}
  }

  // initialize eigentensors eigenvalues
  for ( i = 0; i < 6; i++ ) {
    f_real_eigens[i] = f_im_eigens[i] = 0.0;
  }

    nonPositiveDefinite = 0;
    indeterminate = 0;

  /* Matrix iMatrix(identityMatrix, 6,6);
  setTangentM(iMatrix);
  printf("Identity!!!!\n");
  printTangent();
  calcTarantalo();
  */
 
  

  drawStretch = 1;
}

int GlyphViewer::getIntgrPoint() {
  return intgrPoint;
}

void GlyphViewer::setDrawMethod(unsigned char meth) {
  int method;
  printf("printing method!!\n");
  printf("%c !!!\n", meth);
  fflush(stdout);
  if (meth == 'S') {
    drawMethod = SYMM_PART;
     printf("DRAW SYMM PART OF STIFFNESSS EIGENTENSORS\n");
      fflush(stdout);
    setIntegrPoint( intgrPoint );
  drawStretch = 0;

  }
  else if (meth == 'G' ) {
    drawMethod = STRESS_SINGLE_GLYPH;
    setStress(intgrPoint);
    printf("Draw Stress!\n");
    return;
  }
  
  method = atoi((const char *) &meth);
  if (method >= PETAL0 && method <= REYNOLDS_GLYPH4)
    drawMethod = method;
  else { // calculate eigentensors for 
    if ( method ==  EIGEN_STIFFNESS ) { //8
      drawStretch = 0;
      printf("DRAW STIFFNESSS EIGENTENSORS\n");
      fflush(stdout);
    }
    if ( method == EIGEN_STRETCH ) { //9
      drawStretch = 1;
      printf("DRAW POLAR EIGENTENSORS\n");
      fflush(stdout); 
    }
    if ( method == 83 ) { // capital S
      printf("SYMMETRIC DECOMPOSITION");
      
    }
  
    setIntegrPoint( intgrPoint );
  } 
}



/**
 * precondition: color noormalized based on 4th order tensor
 */

void GlyphViewer::setFlower( double eigSet[6][9] ) {
  int i, j, k;
  double max = -100000000;
  double min =  100000000;
  
  bool positiveDef = true;
	

  for (i = 0; i < 6; i++) {
	// create matrix

	Matrix theMatrix( eigSet[i], 3,3);

	//calculate morph constants
	for( j = 0; j < NUM_PTS; j++ ) { //for each exterior point
	  f_ext[i][j] = einsteinSum2( theMatrix, ptLoc[j] );
	  //	  f_ext[i][j] =   (f_ext[i][j] + 1.0)/2.0; // normalize
	  // track range 
	  if (i > 0) {
	  if ( f_ext[i][j] > max )
		max =  f_ext[i][j];
	  if( f_ext[i][j] < min )
		min =  f_ext[i][j];
	  //if ( f_ext[i][j] < 0.0 )
	  // positiveDef = false;
	  }
	 
	}	
	 
  } //for each petal, calculate point location
  //  printf("min %f max %f\n", min, max);
  // enforce -1 to +1 range as thats how matlab "normalizes"?
  //  normalize f_ext[i][j]
 
  //morph along ray
  // Since the tensors are symmetric, 
  // we use the absolute value and avoid 
  // turning the polygons inside out
  for (i = 0; i < 6; i++) { // for each petal
	  for ( j = 0; j < NUM_PTS; j++ ) {
		for ( k = 0; k < 3; k++ ) {
		  f_extPts[i][j][k] = fabs(f_ext[i][j])*ptLoc[j][k]; 
		  if ( f_extPts[i][j][k] < 0.0 )
			positiveDef = false;
		} //k
	  } //j
	} //i


  //calculate triangle normals & point normals
  for (i = 0; i < 6; i++) { // for each petal
	setNormals(f_extPts[i], f_ptNorm[i], 
			   f_norm[i]); 
  }
}




void GlyphViewer::setDomain( Domain * d ) {
  theDomain = d;
}

/*
void GlyphViewer::setCubeNormals() {
  int y,x, i;
  double normA[3]; //lower left triangle
  double normB[3]; // upper right triangle

  double cube_norm[CUBE_Y][CUBE_X][3];   //patch normal (avg of 2 triangles)

  // calculate normals for each square

  for (y = 0; y < CUBE_Y -1; y++) {
    for ( x = 0; x < CUBE_X -1; x++ ) {
      calcNormal( cube_extPts[y][x],  cube_extPts[y+1][x], 
		  cube_extPts[y+1][x +1], normA);

      calcNormal( cube_extPts[y+1][x+1],  cube_extPts[y][x+1], 
		  cube_extPts[y][x], normB);

      for (i = 0; i < 3; i++ )
	cube_norm[y][x][i] = (normA[i] + normB[i])/2.0;

    }
  }

  // now calculate the point normals.

}

 */

/**
 * calculate normals for each point on the glyph */
void GlyphViewer::setNormals(double extPts[NUM_PTS][3],
		       double ptNorm[NUM_PTS][3], 
			  double norm[NUM_TRIANGLES][3]) {
  int i, j, count;
  //double center[3] = {0.01,0.01,0.01}; 
  
  for(i = 0; i < NUM_TRIANGLES; i++ ) {
		calcNormal( extPts[triPts[i][0]],extPts[triPts[i][1]], 
			extPts[triPts[i][2]], norm[i]);
	  }

  /* calculate normal for each point based on average
   * of surrounding triangles */
  for( i = 0; i < NUM_PTS; i++ ) {
   	for( j = 0; j < 6; j++ ) {
	  ptNorm[i][j] = 0.0;
	}
  }

  for( i = 0; i < NUM_PTS; i++ ) {
	count = 0;
	for( j = 0; j < 6; j++ ) {
	  if ( neighbors[i][j] != -1 ) {
		ptNorm[i][0] = 	ptNorm[i][0] + norm[neighbors[i][j]][0];
		ptNorm[i][1] = 	ptNorm[i][1] + norm[neighbors[i][j]][1];
		ptNorm[i][2] =  ptNorm[i][2] + norm[neighbors[i][j]][2];
		count++;
	  }
	} // 6 neighbors
 
	ptNorm[i][0] /= (double)count;
	ptNorm[i][1] /= (double)count;
	ptNorm[i][2] /= (double)count;
    
	if (ext[i] < 0.0 ) { //reverse the normal
	  ptNorm[i][0] = -ptNorm[i][0];
	  ptNorm[i][1] = -ptNorm[i][1];
	  ptNorm[i][2] = -ptNorm[i][2];

	}

  } //each point
}


/** 
 * draw the triangles 
 * @pre: color and mesh vs. solid already set */
void  GlyphViewer::draw(){
 
  switch( drawMethod ) {
  
  case REYNOLDS_GLYPH4:
        drawCompass();
 	drawE_ijkl();

	break;
  case 	FLOWER_GLYPH:
	drawFlower();
	break;
  case STRESS_SINGLE_GLYPH:{
    drawCompass();
    setStress(intgrPoint);
    drawStressPetal();
   
  }
    break;
  case  PETAL0:
  case PETAL1:
  case PETAL2:
  case PETAL3:
  case PETAL4:
  case PETAL5: 
        drawCompass();  
	drawPetal(drawMethod );
   
  }

} 


/**
 * Draw a single eigentensor as a Reynolds Glyph */
void GlyphViewer::drawPetal(int j) {
  int i;
   float color[4];
  double theta, radius;


  glPolygonMode(GL_FRONT_AND_BACK, GL_FILL);
  glEnable(GL_CULL_FACE);
  glCullFace(GL_BACK);

  eigenColor( j, color);
  glMaterialfv( GL_FRONT_AND_BACK, GL_AMBIENT, color );
  glMaterialfv( GL_FRONT_AND_BACK, GL_DIFFUSE, color );

  glBegin( GL_TRIANGLES );

	for(i = 0; i < NUM_TRIANGLES; i++){

	  glNormal3d(f_ptNorm[j][triPts[i][0]][0], 
		     f_ptNorm[j][triPts[i][0]][1],
			   f_ptNorm[j][triPts[i][0]][2]);
	  
	  
	  // set_color_minus_plus( f_ext[j][triPts[i][0]], 1.0 );
	  glVertex3d(f_extPts[j][triPts[i][0]][0],f_extPts[j][triPts[i][0]][1],
		     f_extPts[j][triPts[i][0]][2] );
	  
	  glNormal3d(f_ptNorm[j][triPts[i][1]][0], 
		     f_ptNorm[j][triPts[i][1]][1],
		     f_ptNorm[j][triPts[i][1]][2]);
	  
	  // set_color_minus_plus( f_ext[j][triPts[i][1]], 1.0 );
	  glVertex3d(f_extPts[j][triPts[i][1]][0],
		     f_extPts[j][triPts[i][1]][1],
		     f_extPts[j][triPts[i][1]][2] );
	  
	  glNormal3d(f_ptNorm[j][triPts[i][2]][0], 
		     f_ptNorm[j][triPts[i][2]][1],
		     f_ptNorm[j][triPts[i][2]][2]);
	  
	  // set_color_minus_plus( f_ext[j][triPts[i][2]], 1.0 );
	  glVertex3d(f_extPts[j][triPts[i][2]][0],
		     f_extPts[j][triPts[i][2]][1],
		     f_extPts[j][triPts[i][2]][2] );	
	}
   glEnd();

   printf("Eigenvalue for petal %d: %f\n",j, f_real_eigens[j]); 

}


/**
 * Draw a single eigentensor as a Reynolds Glyph */
void GlyphViewer::drawStressPetal() {
  int i;

  printf("drawing stress petal!!\n");
  glPolygonMode(GL_FRONT_AND_BACK, GL_FILL);
  glEnable(GL_CULL_FACE);
  glCullFace(GL_BACK);

  glBegin( GL_TRIANGLES );

	for(i = 0; i < NUM_TRIANGLES; i++){

	  glNormal3d(s_ptNorm[triPts[i][0]][0], 
		     s_ptNorm[triPts[i][0]][1],
			   s_ptNorm[triPts[i][0]][2]);
	  
	  
	  set_color_minus_plus( s_ext[triPts[i][0]], 1.0 );
	  glVertex3d(s_extPts[triPts[i][0]][0],s_extPts[triPts[i][0]][1],
		     s_extPts[triPts[i][0]][2] );
	  
	  glNormal3d(s_ptNorm[triPts[i][1]][0], 
		     s_ptNorm[triPts[i][1]][1],
		     s_ptNorm[triPts[i][1]][2]);
	  
	  set_color_minus_plus( s_ext[triPts[i][1]], 1.0 );
	  glVertex3d(s_extPts[triPts[i][1]][0],
		     s_extPts[triPts[i][1]][1],
		     s_extPts[triPts[i][1]][2] );
	  
	  glNormal3d(s_ptNorm[triPts[i][2]][0], 
		     s_ptNorm[triPts[i][2]][1],
		     s_ptNorm[triPts[i][2]][2]);
	  
	  set_color_minus_plus( s_ext[triPts[i][2]], 1.0 );
	  glVertex3d(s_extPts[triPts[i][2]][0],
		     s_extPts[triPts[i][2]][1],
		     s_extPts[triPts[i][2]][2] );	
	}
   glEnd();

    

}

void GlyphViewer::drawFlower() {
  int i, j;
  
  glPolygonMode(GL_FRONT, GL_FILL); 
  glBegin( GL_TRIANGLES );
  for (j = 1; j < 6; j++) {
    set_color_by_param((double)j/6.0, 1.0);
 
	for(i = 0; i < NUM_TRIANGLES; i++) {

	glNormal3d(f_ptNorm[j][triPts[i][0]][0], f_ptNorm[j][triPts[i][0]][1],
			   f_ptNorm[j][triPts[i][0]][2]);
    

	set_color_minus_plus( f_ext[j][triPts[i][0]], 1.0 );
	glVertex3d(f_extPts[j][triPts[i][0]][0],f_extPts[j][triPts[i][0]][1],
			   f_extPts[j][triPts[i][0]][2] );

	glNormal3d(f_ptNorm[j][triPts[i][1]][0], f_ptNorm[j][triPts[i][1]][1],
		 f_ptNorm[j][triPts[i][1]][2]);

	set_color_minus_plus( f_ext[j][triPts[i][1]], 1.0 );
	glVertex3d(f_extPts[j][triPts[i][1]][0],f_extPts[j][triPts[i][1]][1],
			   f_extPts[j][triPts[i][1]][2] );

	glNormal3d(f_ptNorm[j][triPts[i][2]][0], f_ptNorm[j][triPts[i][2]][1],
			   f_ptNorm[j][triPts[i][2]][2]);

	set_color_minus_plus( f_ext[j][triPts[i][0]], 1.0 );
	glVertex3d(f_extPts[j][triPts[i][2]][0],f_extPts[j][triPts[i][2]][1],
			   f_extPts[j][triPts[i][2]][2] );	
	}
  }
  glEnd();

}

/**
 * @pre: orthographic projection to avoid squeezing the glyph */
void GlyphViewer::drawE_ijkl() {
  int i;
 
  glPolygonMode(GL_FRONT_AND_BACK, GL_LINE); 
  glBegin( GL_TRIANGLES );
  for(i = 0; i < (NUM_TRIANGLES); i++ ) {

	glNormal3d(ptNorm[triPts[i][0]][0], ptNorm[triPts[i][0]][1],
			   ptNorm[triPts[i][0]][2]);
   	set_color_by_param_shiny( ext[triPts[i][0]], 0.1 );
	
	glVertex3d(extPts[triPts[i][0]][0],extPts[triPts[i][0]][1],
			   extPts[triPts[i][0]][2] );

	//--------------------------------------

	glNormal3d(ptNorm[triPts[i][1]][0], ptNorm[triPts[i][1]][1],
		 ptNorm[triPts[i][1]][2]);
	set_color_by_param_shiny( ext[triPts[i][1]], 0.1 );
	glVertex3d(extPts[triPts[i][1]][0],extPts[triPts[i][1]][1],
			   extPts[triPts[i][1]][2] );


	//--------------------------------------

	glNormal3d(ptNorm[triPts[i][2]][0], ptNorm[triPts[i][2]][1],
			   ptNorm[triPts[i][2]][2]);
	set_color_by_param_shiny( ext[triPts[i][2]], 0.1);
	glVertex3d(extPts[triPts[i][2]][0],extPts[triPts[i][2]][1],
			   extPts[triPts[i][2]][2] );
	
  }
  glEnd();
 
}



double  GlyphViewer::getCenterAndMaxDim( double center[3] ){
  center[0] = center[1] = center[2] = 0;

  return 3.0; 
  
}

int  GlyphViewer::getNumElements() {
  return 1;
}

void  GlyphViewer::setElement( int ele ) {
  printf("GlyphViewer::setElement %d\n", ele);
  fflush(stdout);
  theElement = NULL;
  theElement = theDomain->getElement( ele );
  eleNum = ele;
}


void  GlyphViewer::setStress( int theNum ) {
  int i;
  Response * theResponse = NULL;
  Information eleInfo(1.0);
  char * argStress[] = { "stress" };
  double norm;

  if (theNum == -1 || theElement == NULL ) 
    return;

  // actual hit
    intgrPoint  = theNum;
    
    for ( i = 0; i < 9; i++ ) {
      stress[i]  = 0.0;
    }
    
 
    if ( theElement->getNumExternalNodes() == 8 ) {
	//try get stress
	theResponse = 
	  theElement->setResponse((const char **)argStress, 1, eleInfo); 
	if (  theResponse ) {
	  theResponse->getResponse();
	  
	  Information &theInfo = theResponse->getInformation();
	  const Vector &eleData = theInfo.getData();
	  
	  if ( eleData.Size() == 48 ) {
	    // xx yy zz xy xz yz
	    
	    stress[0] = eleData(intgrPoint*6);
	    stress[4] = eleData(intgrPoint*6+1);
	    stress[8] = eleData(intgrPoint*6+2);
	    stress[1] = stress[3] = eleData(intgrPoint*6+3);
	    stress[2] = stress[6] = eleData(intgrPoint*6+4);
	    stress[5] = stress[7] = eleData(intgrPoint*6+5);
	    
	  } // size 48
	  
	  else if ( eleData.Size() == 49 ) {
	    
	    stress[0] = eleData(intgrPoint*6+1);
	    stress[4] = eleData(intgrPoint*6+2);
	    stress[8] = eleData(intgrPoint*6+3);
	    stress[1] = stress[3] = eleData(intgrPoint*6+4);
	    stress[2] = stress[6] = eleData(intgrPoint*6+5);
	    stress[5] = stress[7] = eleData(intgrPoint*6+6);
	    
	    
	  } // size 49
	} // non-null response
	// NORMALIZE THE TENSOR, OTHERWISE THE GLYPH WILL 
	// BE HUGE
	norm = 0;
	for ( i = 0; i < 9; i++)
	  norm = norm + (stress[i]*stress[i]);
	norm = sqrt(norm);
	for ( i = 0; i < 9; i++)
	  stress[i] = stress[i]/norm;

	printf("\nSTRESS--------\n");
	for ( i = 0; i < 9; i++) {
	  if (i%3 == 0)
	    printf("\n");
	  printf("%f ", stress[i]);
	}
	printf("\n-------------\n");
	setPetal(stress);
    } // 8 nodes
}



/**
 * precondition: color noormalized based on 4th order tensor
 */

void GlyphViewer::setPetal( double eigSet[9] ) {
  int j, k;
  double max = -100000000;
  double min =  100000000;
  
  bool positiveDef = true;
	


  // create matrix
  Matrix theMatrix( eigSet, 3,3);
  
  //calculate morph constants
  for( j = 0; j < NUM_PTS; j++ ) { //for each exterior point
    s_ext[j] = einsteinSum2( theMatrix, ptLoc[j] );

   
    if ( s_ext[j] > max )
      max =  s_ext[j];
    if( s_ext[j] < min )
      min =  s_ext[j];
  }		

  //morph along ray
  // Since the tensors are symmetric, 
  // we use the absolute value and avoid 
  // turning the polygons inside out

  for ( j = 0; j < NUM_PTS; j++ ) {
    for ( k = 0; k < 3; k++ ) {
      s_extPts[j][k] = fabs(s_ext[j])*ptLoc[j][k]; 
      if (  s_extPts[j][k] < 0.0 )
	positiveDef = false;
    } //k
  } //j



  //calculate triangle normals & point normals
  setNormals(s_extPts, s_ptNorm, s_norm);//, s_norm); 
}




/**
 * Set the integration point, get its stiffness tensor, 
 * calculate eigentensors and their eigenvalues, update glyph
 */
void  GlyphViewer::setIntegrPoint( int theNum ) {
  int i,j;

  printf(" GlyphViewer::setIntegrPoint %d\n", theNum);
  if (theNum != -1 ) { // actual hit
	intgrPoint  = theNum;

	// I'm guessing the default rotation (zero degrees)
	// is the identity. (true for 3x3)
	// not sure if it si same for 6x6
	for (i = 0; i < 6; i++) {
	  for (j = 0; j < 6; j++) {
	    if ( i != j )
	      rotation[i][j] = 0;
	    else 
	      rotation[i][j] = 1;
	  }
	}
	  
	// do the big cast
       if ( theElement->getClassTag() == ELE_TAG_EightNodeBrick ) {
	 setStress(theNum);

	  EightNodeBrick * el = (EightNodeBrick *) ( theElement );
	 
	  tensor theTensor = (el->getTangentTensor( intgrPoint ));
	  //	  theTensor.print();
	 
	  setTangentT( theTensor );
	 
	
	  morphGlyphPoints4();
	 
	  setNormals(extPts, ptNorm, norm);

	  Matrix * theMatrix = new Matrix(6,6);
	  Matrix * theTranspose = new Matrix(6,6);
	  Tensor2MatrixSysR4(theTensor, *theMatrix); // tensor to a helnwein matrix
	  if (drawMethod == SYMM_PART ) {
	    transposeMatrix( *theMatrix, *theTranspose );
	    (*theMatrix) = (*theMatrix) + (*theTranspose); 
	    (*theMatrix) = (*theMatrix)*0.5;
	  }

	  //  printf("Stiffness Eigens\n\n\n");
	  double * vr = new double[36];
	  double * wr = new double[6];
	  double * wi = new double[6];
	  printf("\n\nEigentensors of symmetric part of stiffness\n");
	  if ( calcSymmEigens( * theMatrix, 6, vr, wr) != 0 ) //{
	    printf("GlyphViewer::setIntegrPoint: error calculating eigentensors\n");
	  
	  for (int ii = 0; ii <6; ii++) {
	    f_real_eigens[ii] = wr[ii];
	  }

	
	    if (! drawStretch ) {	 // stiffness eigentensors
	    
	    sixVectorToSymMatrix( vr );
	    setFlower(f_matrices);
	    delete[] vr;
	    delete[] wr;
	    delete[] wi;
	    //testSmallest(theTensor);
	  }
	  
	    else { // stretch eigentensors
	    
	    // check for singularity, check whether we are still at initial
	    // stiffness
	    
	    int dim[4] = {3,3,3,3};
	    Tensor init( 4,dim, 0.0 );
	    //	    Tensor inverse( 4,dim, 0.0 );
	    Tensor curr( 4,dim, 0.0 );;
	    
	    Matrix m_init(6,6);
	    Matrix m_elastic(6,6);
	    Matrix m_curr(6,6);
	    Matrix m_inverse(6,6);
	    Matrix m_symm(6,6);
	    /*
	    init =  el->getInitialTangent();
	    Tensor2MatrixSysR4(init, m_init);

	    printf("Initial stiffness!!\n");
	    opserr << m_init;*/
	    

	    // curr = el->ElasticStiffnessTensor( intgrPoint );
	    
	    //Tensor2MatrixSysR4(curr, m_elastic);

	    //  printf("Current Elastic Stiffness!!!!!!!!!!!!!!!\n");
	    // opserr << m_elastic;

	    curr = el->getTangentTensor( intgrPoint );
	    // curr.print();
	    
	    Tensor2MatrixSysR4( curr, m_curr );

	    printf("Current stiffness!\n");
	    opserr << m_curr;

	    symmetrizeMatrix(m_curr, m_symm ); // get symmetric part

	    //  printf("Symmetric part of curr. stiffness!\n");
	    //opserr << m_symm;

	    //  theTensor = el->dFods( intgrPoint );
	    //theTensor.printshort("dFods!!!!!!!!\n");

	   
	    //theTensor = el->dQods( intgrPoint );
	    //theTensor.printshort("dQods!!!!!!!!\n");

	    //m_diff = m_curr*m_inverse; // C.inverse(C_0); 
	    // we really need matrix multiply. Ho do we do it for tensors?
	    //printf("Difference from original!!\n");
	    //opserr << m_diff;
	    
	    testBifurcation(m_symm);

	    printf("\nnonPositiveDefinite in symmetric part %d\n",nonPositiveDefinite);
	    printf("indeterminate in symmetric part %d\n",indeterminate) ;

	    testBifurcation(m_curr);

	    printf("\nnonPositiveDefinite in raw stiffness %d\n",nonPositiveDefinite);
	    printf("indeterminate in raw stiffness %d\n\n",indeterminate) ;

	    

	    //  printf("polar decomposition\n");
	    Matrix stretch;
	    Matrix rotate;
	    
	    vanGelderPolar( m_curr, rotate, stretch, 1e4);

	    for (i = 0; i < 6; i++) {
	      for (j = 0; j < 6; j++ ) {
		rotation[i][j] = rotate(i,j); // Matrices column-wise ordered?
	      }
	    }


	    // now, get eigentensors of stretch
	    printf("\nEigenTensors of Raw Stiffness Polar Stretch\n");
	    if ( calcSymmEigens( stretch, 6, vr, f_real_eigens ) == 0 ) {// success
	      
	   
	      
	      sixVectorToSymMatrix( vr ); // actually turns them into 3x3 matrices  
	     

	      double iso, txc, txe, shear;
	      //  eigenMode();
	      // may reverse direction of eigentensor
	      eigenMode2( iso, txc, txe, shear); 
	      setFlower(f_matrices);
	      
	    }
	    else printf("GlyphViewer::setIntegrPoint: error calculating eigentensors\n");	  
	    
	    }
	  
       }
       else if ( theElement->getClassTag() == ELE_TAG_Brick ) {
	 Brick * el = (Brick *) ( theElement );
	 
	 Matrix theMatrix = el->getTangent( intgrPoint );
	 setTangentM( theMatrix );
	 morphGlyphPoints4();
	 setNormals(extPts, ptNorm, norm);
       }
       else if ( theElement->getClassTag() == 
		 ELE_TAG_EightNodeBrick_u_p_U ) {
	 EightNodeBrick_u_p_U * el = (EightNodeBrick_u_p_U *) ( theElement );
	 /*  
	  Matrix theMatrix = el->getTangent( intgrPoint );
	  setTangentM( theMatrix );
	  morphGlyphPoints4();
	  setNormals(extPts, ptNorm, norm);*/
	  tensor theTensor = (el->getTangentTensor( intgrPoint ));
	  // theTensor.print();
	  setTangentT( theTensor );
	  morphGlyphPoints4();
	  setNormals(extPts, ptNorm, norm);


	  // now, convert the tensor to a helnwein matrix
	  Matrix * theMatrix = new Matrix(6,6);
	  Tensor2MatrixSysR4(theTensor, *theMatrix);
	  double * vr = new double[36];
	  if ( calcAsymmEigens( * theMatrix, 6, vr, f_real_eigens, 
				f_im_eigens) == 0 ) {// success
	    sixVectorToSymMatrix( vr );
	    setFlower(f_matrices);
	  }
	  else printf("GlyphViewer::setIntegrPoint: error calculating eigentensors\n");

	  printf("Success setting integrpt\n");
	  fflush(stdout);
	  }
	//make elasticcrossanisotropic material

  } //hit
  else { 
    printf("miss\n:");
  }

}


/**
 * Load data from new time step */
void GlyphViewer::updateState( int step ) {
  printf(" GlyphViewer::updateState ele %d intgr %d\n",  eleNum,  intgrPoint);
  fflush(stdout);
  if ( eleNum != -1 &&  intgrPoint > -1 ) {
    printf("GlyphViewer::updateState, setElement?\n");
    fflush(stdout);
    setElement( eleNum );
    setIntegrPoint( intgrPoint );
  }
  else {
    printf("no update for ele %d integration point %d\n", eleNum, intgrPoint);
    fflush(stdout);
  }
  commitStep = step;
}

/**
 * C_ij n_i n_j for appropriate number of indeces.
 * Constract the tensor to a scalar with each normalized
 * vector on the sphere. 
 */
void GlyphViewer::morphGlyphPoints4(){
  int i, j;

  for( i = 0; i < NUM_PTS; i++ ) {
    ext[i] = einsteinSum4( ptLoc[i] ); 
  }


  sizeRange.setData( ext, NUM_PTS);

  //    sizeRange.setBasis(77774, true,false);
  // sizeRange.scaleDataSet(); // works!!
  

  // sizeRange.scaleAndSetRange( ext, NUM_PTS, true, false); 
  // check whether there is an old scale value. If so use it.
  // otherwise create a new scale value
  if ( sizeRange.getBaseValue() == 1.0 ){ //default
                               // scalar, length, posDef, logScaling
    sizeRange.scaleAndSetRange( ext, NUM_PTS, true, false ); 
    printf(
       "\n\n\n*******************************\nsizeRange.scaleAndSetRange\n");
    fflush(stdout);
    }
  else {
    printf("\nScaling By curr factor:");
    fflush(stdout);
    sizeRange.scaleDataSet();
    }
  
 
 
  // normalizeMinusPlus( ext, NUM_PTS ); // worry about   //singles laterl;l 
  for ( i = 0; i < NUM_PTS; i++ ) {
    for ( j = 0; j < 3; j++ ) {     
      extPts[i][j] = ext[i]*ptLoc[i][j];	       
    }

  } // i

}  

/**
 * Use for single Reynolds Stress/Strain glyph, not flower
 */
void GlyphViewer::morphGlyphPoints2(  Matrix theMatrix ){

  int i, j;

  for( i = 0; i < NUM_PTS; i++ ) {
	ext[i] = einsteinSum2( theMatrix, ptLoc[i] ); 
	// ext[i] = HWYglyph( theMatrix, ptLoc[i] ); 
  }
  normalizeColorScale( ext, NUM_PTS, true ); // worry about   //singles later


  
  for ( i = 0; i < NUM_PTS; i++ ) {
	for ( j = 0; j < 3; j++ )
	  extPts[i][j] = ext[i]*ptLoc[i][j];
  } 

}

/**
 * Tarantalo's unrolling of the 3x3x3x3 tensor into 9x9 
 * unrolls current tangent tensor
 */
void GlyphViewer::calcTarantalo() {
  int i,j,k,l;

  // clear tensor
 for( i = 0; i < 9; i++ ) {
	for ( j = 0; j < 9; j++ ) {  	
	  tarantalo[i][j] = 0;
	}
 }
   

  for( i = 0; i < 3; i++ ) {
	for ( j = 0; j < 3; j++ ) {
	  for ( k = 0; k < 3; k++ ) {
		for ( l = 0; l < 3; l++ ) {	

		  if( i == j && k == l ) {  //upper left sector 
			tarantalo[i][k] = tangent[i][j][k][l];
		  }

		  else if( i == j ) { // remaining upper sectors
            //i = row number  
			//if kl == 01, cols 3 & 6
            //"        12  cols cols 4 & 7
            //"        20  cols 5 & 8
			
			//calculate 1st col: 1/sqrt(2) tang[i][j][k][l]+tang[i][j][l][k]
			tarantalo[i][k+3] = ONE_OVER_SQRT_2*(2*tangent[i][j][k][l]);
           
            //2nd col:zero!!
			// 1/sqrt(2) tang[i][j][k][l]-tang[i][j][l][k] = 0 
			//because of minor symmetry
		  }

		  else if (k == l) { //left hand column
            //k = column number
            //if ij = 01 rows 3 & 6
            // "     =12  rows 4 & 7
            // "     = 20 rows 5 & 8 
            //calculate 1st row 1/sqrt(2)  tang[i][j][k][l]+ tang[j][i][k][l]
            tarantalo[i+3][k] = ONE_OVER_SQRT_2*(2*tangent[i][j][k][l]);
			//2nd row: zero. tang[i][j][k][l] == tang[j][i][k][l] 
			//so there differnce is 0.
		  }

		  else  { // lower right quad
	
            // 1/2 (tangent[i][j][k][l] + tangent[i][j][l][k] + tangent[j][i][k][l] + tangent[j][i][l][k])
            // equiv. to 2*tangent[i][j][k][l] because of minor symmetry
			//placing: 1212 1223 1231 2312 2323 2331 3112 3123 3131
            if (i== 0 && j == 1) {
			  if (( k == 0 && l == 1) || (  k == 1 && l == 2) 
				  ||( k == 2 && l == 0 ))
				tarantalo[i+3][k+3] = 2*tangent[i][j][k][l];
			} // row 1 of lower right 3 x 3
			else if ( i == 1 && j == 2 ) {
			  if (( k == 0 && l == 1) || (  k == 1 && l == 2) 
				  ||( k == 2 && l == 0 ))
				tarantalo[i+3][k+3] = 2*tangent[i][j][k][l];		
			}//middle row
			else if ( i == 2 && j == 0 ) {
			  if (( k == 0 && l == 1) || (  k == 1 && l == 2) 
				  ||( k == 2 && l == 0 ))
			  	tarantalo[i+3][k+3] = 2*tangent[i][j][k][l];
			}//bottom row
		  }
		  
		}//i
	  }//j
	}//k
  }//l
  /*
    printf("Tarantalo\n");
 for( i = 0; i < 6; i++ ) {
	for ( j = 0; j < 6; j++ ) {  	
	  printf("%f ", tarantalo[i][j]);
	}
	printf("\n");
	}*/

}

/**
 * Stress or strain tensor. Hashah's method to show shear
 * sigma_S = sqrt{T^2 - sigma^2_N}  = sqrt{T_iTi - sigma^2_N} 
 * T_i = sigma_{ij}n_j
 * @returns the result of the formula or 0 if not a 2nd order tensor
 */
double GlyphViewer::HWYglyph( Matrix theMatrix,double n[3] ) {
  int nRows;
  int nCols;
  int i, j;
  double sum, Tdot;
  double T[3];

  nRows = theMatrix.noRows();
  nCols = theMatrix.noCols();
  sum = T[0] = T[1] = T[2] = Tdot = 0.0;
  

  if ( nRows == 3 && nCols == 3 ) { // 3 x 3 matrix
	sum = einsteinSum2( theMatrix, n );

	//get vector T
	for( i = 0; i < nRows; i++ ) {
	  for ( j = 0; j < nCols; j++ ) {
		T[i] = T[i] + theMatrix(j,i)*n[j];//matrix class indexing reversed
	  }
	 
	}

    //dot T with T	
	for( i = 0; i < nRows; i++ ) {
	  Tdot = Tdot + (T[i]*T[i]);
	}
	// sqrt of difference of squares
	Tdot = Tdot - (sum*sum);

	Tdot = sqrt(Tdot);

	return Tdot;
  }
  return 0;
}



/** Einstein summation over second order tensor 
 * with a vector 
 * @returns the result of the formula or 0 if not a 2nd order tensor
*/
double GlyphViewer::einsteinSum2( Matrix theMatrix, double n[3] ) {
  int nRows, nCols, i, j;
  double sum;
  nRows = theMatrix.noRows();
  nCols = theMatrix.noCols();
  sum = 0.0;
  
  if ( nRows == 3 && nCols == 3 ) { // 3 x 3 matrix
	for( i = 0; i < nRows; i++ ) {
	  for ( j = 0; j < nCols; j++ ) {
		sum = sum + theMatrix(j,i)*n[i]*n[j];//Matrix class indexing reversed
	  }
	}
	return sum;
  } // 3 x 3 matrix
  return 0; // error
}  




/** Einstein summation over fourth order tensor 
 * with a vector */
double GlyphViewer::einsteinSum4(  double n[3] ) {
  int  i, j, k,l;
  double sum;

  sum = 0.0;
  

  for( i = 0; i < 3; i++ ) {
	for ( j = 0; j < 3; j++ ) {
	  for ( k = 0; k < 3; k++ ) {
		for ( l = 0; l < 3; l++ ) {	
		 
		   sum = sum + tangent[i][j][k][l]*n[i]*n[j]*n[k]*n[l];
		  // vanGelder Sum
		  // sum = sum + tangent[i][j][k][l]*n[i]*(j==i)*n[k]*(l==k);
		} // l
	  } // k 
	} // j
  } // i
 
 
  return sum;
}  




/**
 * copy from BJtensor into standard 4D array 
 * ( 3 x 3 x 3 x 3 ) 
 */
void GlyphViewer::setTangentT( tensor theTensor )
{
  int i,j,k,l;
  //theTensor.print();
  // cval is the true tensor index, and it runs from 1 to dim
  for( i = 0; i < 3; i++ ) {
	for ( j = 0; j < 3; j++ ) {
	  for ( k = 0; k < 3; k++ ) {
		for ( l = 0; l < 3; l++ ) {
		  tangent[i][j][k][l] = theTensor.cval(i+1,j+1,k+1,l+1);
		}
	  }
	}
  }

}
 

/** map 6 x 6 Matrix => 3 x 3 x 3 x 3 tensor with minor symmetry */
void GlyphViewer::setTangentM( Matrix theMatrix ) {
  int i, j, k, l, ii, jj;
  //  printf("\n\n:setTangentM:\n");
  //theMatrix.Output( opserr );
  for ( ii = 0; ii < 6; ii++ ) {
     for ( jj = 0; jj < 6; jj++ )  {
 
           indexMap( ii, i, j ) ;
           indexMap( jj, k, l ) ;
 
           tangent[i][j][k][l]  = theMatrix( ii, jj ); 
	   //	   printf("Matrix[%d][%d]= %f = tangent[%d][%d][%d][%d]\n", ii, jj,  theMatrix( ii, jj ), i,j,k,l);
           //minor symmetries 
	   
	   tangent [j][i][k][l] = theMatrix(ii,jj);
	   tangent [i][j][l][k] = theMatrix(ii,jj);
           tangent [j][i][l][k] = theMatrix(ii,jj); 
	   
	   //major symmetry
	   tangent[i][j][k][l] = theMatrix( ii,jj );
	   tangent[k][l][j][i] = theMatrix( ii, jj );
	   tangent[l][k][j][i] = theMatrix( ii, jj );

    } // end for jj
   } // end for ii
  
 


}

//Matlab 9 x 9 matrix
void GlyphViewer::printTangent() {
  int i,j,k,l;
  printf("T = [");
  for (i = 0; i < 3; i++) {
	for ( j = 0; j < 3; j++ ) {
	  for (k = 0; k < 3; k++ ) {
		for (l = 0; l < 3; l++ ) {
	
		   printf("%f ", tangent[i][j][k][l]);
		}//l
	  }//k
	  printf(" ; "); //fgsbh
    }
  }
  printf("]");
}

//Matlab 6 x 6 matrix
void GlyphViewer::printMatrix(Matrix theMatrix) {
  int i,j;

  printf("M = [");
  for (i = 0; i < 6; i++) {
	for ( j = 0; j < 6; j++ ) {
	  printf("%f ", theMatrix(j,i));//Matrix class indexing reversed
	}
	printf("; ");
  }
  printf("]");
}

 // matrix to tensor mapping
  //  Matrix      Tensor
  // -------     -------
  //   0           0 0
  //   1           1 1
  //   2           2 2   
  //   3           0 1  ( or 1 0 )
  //   4           1 2  ( or 2 1 )
  //   5           2 0  ( or 0 2 ) 
void GlyphViewer::indexMap( int mIndex, int & i, int & j ) {

  switch( mIndex ) {
  case 0:
	i = j = 0; // xx
	break;

  case 1:
	i = j = 1; // yy
	break;

  case 2:
	i = j = 2; // zz
	break;

  case 3:   //  yx
	i = 0; //col
	j = 1; //row
	break;

  case 4: // yz, zy
	i = 1; //col
	j = 2;  //row
	break;

  case 5: // zx, xz
	i = 2;
	j = 0;
	break;
  default:
	printf ("Error GlyphViewer::index map: called with %d, should be 0 to 5\n",
			mIndex);
			
  }
}

/**
 * take data that has both positive & negative values, 
 * normalize and center over 0-1. Range may be preset or calculated
 * on the fly
 */
void GlyphViewer::normalizeColorScale( double * scalar, int length, bool setRange ) {
  int i;
  double glyphMax, glyphMin, range;
  bool positiveDefinite;
  positiveDefinite = false;
  if (setRange) {
	glyphMax = -100000000;
	glyphMin =  100000000;
	

	
	/* establish range */
	for ( i = 0; i < length; i++ ) {
	  if ( scalar[i] > glyphMax )
		glyphMax = scalar[i];
	  if( scalar[i] < glyphMin )
		glyphMin = scalar[i];
	  if (scalar[i] < 0.0 )
		positiveDefinite = false;
	} 
	
	
	if(fabs(glyphMin) > fabs(glyphMax))
	  range = fabs( glyphMin );
	else range = fabs(glyphMax);
  }

  for( i = 0; i < length; i++) {
	if( positiveDefinite ) {
	  scalar[i] = scalar[i]/range; /* normalize */
	if (scalar[i] < 0.0 )
	  scalar[i] = 0.0;  // repair error

	printf("FOUND BROKEN!!!!!!!!!!!!!!!!!!!!!!!!\n");
	}
	else {
	  scalar[i] = scalar[i] + range; /* shift over 0-1 */
	  scalar[i] = scalar[i]/(range*2); /* normalize */
	}

  }
}

/**
 * normalize scale over +1 to -1. Assume abs(max) = 1.0 and
 * -(abs(max)) = -1.0
 */

void GlyphViewer::normalizeMinusPlus( double * scalar, int length) {
  int i;
  double glyphMax, glyphMin, range;
  glyphMax = -100000000;
  glyphMin =  100000000;
	
  bool positiveDefinite = true;

  /* establish range */
  for ( i = 0; i < length; i++ ) {
	if ( scalar[i] > glyphMax )
	  glyphMax = scalar[i];
	if( scalar[i] < glyphMin )
	  glyphMin = scalar[i];
   
	if (scalar[i] < 0.0 )
	  positiveDefinite = false;
	
  }
  printf("\n*****************\nglyph max = %f min = %f\n",glyphMax, glyphMin);
	
  if ( positiveDefinite ) {
	range = glyphMax;
  }
	
  else {
	if(fabs(glyphMin) > fabs(glyphMax)){
	  range = fabs( glyphMin );
	}
	else{
	  range = glyphMax;
	}
  }

  for( i = 0; i < length; i++) {
	scalar[i] = scalar[i]/range; // retain sign (+/-)
  }

}

/**
 *  show X,Y,Z axis legend with icon in upper left corner */
void  GlyphViewer::drawCompass() {
  double base,height, len;
  double compassSize = 1.3;

  GLUquadric * cylinder;
  cylinder = gluNewQuadric();
  glLineWidth(2);

  height = compassSize;
  base = compassSize/8;

  len = height/3;
  set_color_by_param_shiny( 1.0, 0.2 ); //red
  // set_saturation_by_param(0.0, 1.0); // black
  Vector v = theDomain->getPhysicalBounds();
  glPolygonMode(GL_FRONT, GL_FILL); 

  //------- Z axis ---------
   set_color_by_param_shiny( 0, 0.2); //blue
  glPushMatrix();
  
  //  glTranslatef( compassLoc[X], compassLoc[Y], compassLoc[Z]);
   gluCylinder( cylinder , base/7 , base/7, 
	       height*.85, 8, 1 );
   
 glPushMatrix();
  glTranslatef(0,0, height*.85);
  glutSolidCone(base/3, base/2,6,1);
  glPopMatrix();

 
  glPopMatrix();
  
  //---------- X axis ---------
  set_color_by_param_shiny( 1.0, 0.2 ); //red
  glPushMatrix();
  //  glTranslatef( compassLoc[X], compassLoc[Y], compassLoc[Z]);
  
  glRotatef(90, 0.0, 1.0, 0.0); //rotate around Y axis, swap x & z
  gluCylinder( cylinder , base/7 , base/7, 
	       height*.85, 8 , 1);
  glPushMatrix();
  glTranslatef(0,0, height*.85);
  glutSolidCone(base/3, base/2,6,1);
  glPopMatrix();
 

  glPopMatrix();

  //------- Y axis --------
  
   set_color_by_param_shiny( 0.5, 1.0 ); //green

  glPushMatrix();
  //  glTranslatef( compassLoc[X], compassLoc[Y], compassLoc[Z]);
  
  glRotatef(270, 1.0, 0.0, 0.0); //rotate around X axis
  gluCylinder( cylinder , base/7 , base/7, 
	       height*.85 , 8 , 1 );

  glPushMatrix();
  glTranslatef(0, 0,height*.85 );
  glutSolidCone(base/3, base/2,6,1);
  glPopMatrix();

 
  glPopMatrix();
  set_color_by_param_shiny(1.0,1.0); // change params re lighting
  
  delete cylinder; // clean up
  glLineWidth(1);
}


/**
 * precondition: tarantalo matrix has been calculated
 * from stiffness tensor. Note: eventually we will have a separate 
 * class to wrap the LAPACK/BLAS code. D for double precision
 * we want the "right" eigenvectors, i.e. A*v_j = lamda*v_j so JOBVL
 * is NO and JOBVR is YES

  resulting eigenvectors are normalized.

 */
#include <f2c.h>
extern "C" int dgeev_(  char * JOBVL, char * JOBVR, int * N, double * A, 
            int * LDA, double * WR, double * WI, double * VL, int * LDVL, 
           double * VR, int * LDVR, double * WORK, int * LWORK, int *INFO);
/**
 * calculate eigens for 6 x 6 or 3 x 3. method will allocate eigenvectors
 * real and imaginary parts of eigenvalues. Caller must free.
 * @param vr:  right eigenvectors
 * @param wr: real part of eigenvalues
 * @param wi: imaginary part of eigenvalues
 * @returns 0 for correct execution 
 */
int  GlyphViewer::calcAsymmEigens(  Matrix  theMatrix, int leadingDim, double * vr, double * wr, double * wi )
{
  /* now, we can allocate and call dgeev. so, let's try a little example that 
   * we know in adavance the answer
   *  |  5  7 |
   *  | -2 -4 | L1 = -2 EV1 = |1 -1|^T  L2 = 3 EV2 = |7 -2|^T
   */
  int i,j,k;
  int  n; // order of our (2D) matrix
  int lda; // leading dimension of A, in this case it's 2 (from 2 x 2)
  int ldvl; // leading dimension of vl, again 2
  int ldvr; // leading dimension of vl, again 2
  int lwork; // work array size, 4*n
  int info; // 0= success, nonzero failure
  
  // double * wr; // real part of eigenvalues
  //double * wi; // imaginary part of eigenvalues
  double * vl; // left eigenvectors
  // double * vr;  // right eigenvectors
  double * work; // scratch space
  double * a; // the initial matrix, that will be overwritten!!
  char jobvl = 'N';
  char jobvr = 'V';

  
  lda = ldvl = ldvr = n = leadingDim;
  lwork = 4*n;
  

  // allocate arrays
  //  wr = new double[lda];
  //  wi = new double[lda];
  vl = new double[lda*lda];
  // vr = new double[lda*lda];
  work = new double[lwork];


 /* LAPACK uses column first ordering, e.g.
   * |  5  7 |
   * | -2 -4 | == [ 5 -2 7 -4 ] */

  a = new double[lda*lda];
  k = 0;
  // copy the Matrix into a
  for( i = 0; i < lda; i++) {
    for(j = 0; j < lda; j++) {
      a[k] =  theMatrix(i,j); // matrix also stored in column order
      k++;
    }
  }
  

  if( wr == NULL || wi == NULL || vl == NULL || vr == NULL
      || work == NULL || a == NULL ) {
          printf("calcStiffnessEigens, memory allocation failed\n");
    return -1;
  }

 
  // debug
  printf("The original matrix!!!\n");
  opserr << theMatrix;

  /*printf("a, the array representing the matrix (columnwise)\n\n|");
  for (i = 0; i < lda*lda; i++ ){
    printf("%f ",a[i]);
    if (i %lda == 0) printf("|\n|");
    }*/
  printf("\n");

  // ok, let's run it
  dgeev_( &jobvl, &jobvr, &n, a, 
            &lda, wr, wi, vl, &ldvl, 
		       vr, &ldvr, work, &lwork, &info);


  printf("run returned %d\n", info);
  printf("real eigen vals: %f %f %f %f %f %f\n", 
	 wr[0], wr[1], wr[2], wr[3],wr[4],wr[5]);

  printf("imaginary %f %f %f %f %f %f \n",wi[0],wi[1],wi[2], wi[3], wi[4],
   wi[5]);
  printf("right eigen vectors\n");
  for(i = 0; i < 36; i++ ) {
    if (i %6 ==0 )
      printf("\n[%d: ", i/6);
    printf(" %f ", vr[i]);
    

  }
  printf("\n");
  //print some answers

 

  // delete wr;
  // delete wi;
  delete[] vl;
  // delete vr;
  delete[] work;
  delete[] a;

  // at this point we probably want to do setFlower to get 
  // the eigenvectors & also print the numbers to look at.
  // also, change setFlower to write to a local variable.


  return info;
}

// ok, now for symmetric tensors (stress, strain, some stiffness)

extern "C" int dsyev_( char * JOBZ, char * UPLO, int * N, 
		       double * A, int * LDA, double * W, double * WORK, 
		       int * LWORK, int * INFO );
/**
 * Calculate eignevectors & eigenvalues of symmetric square matrix
 * Caller must allocate eigenvector & eigenvalue arrays
 * @param leadingDim, the dimension of the square matrix
 * @param w array for eigenvalues
 * @param a array for eigenvectors
 * @returns 0 for success, i < 0, ith argument to dsyev call illegal, 
 * i > 0 for fail to converge  
 */
int  GlyphViewer::calcSymmEigens(  Matrix  theMatrix, int leadingDim, 
				   double * a, double * w) {
  char jobz = 'V'; //eigenvectors and eigenvalues
  char uplo = 'U'; // in reality we will give the whole matrix
  int n, lda, lwork, info, i, j,k;
  double * work;

  n = lda = leadingDim;
  lwork = 8*n; // lwork >= 3*n-1
  info = -1; 

  work = new double[lwork];

  // copy the matrix into 'a'
  k = 0;
  // copy the Matrix into a
  for( i = 0; i < lda; i++) {
    for(j = 0; j < lda; j++) {
      a[k] =  theMatrix(i,j); // matrix also stored in column order
      k++;
    }
  }

  dsyev_( &jobz, &uplo, &n, a, &lda, w, work, 
	  &lwork, &info );

  // do print test
  printf("\n^^^^^eigenvalues^^^^^^^\n[");
  for (k = 0; k < lda; k++)
    printf("%f ", w[k]);
  printf("]\n");
  //  printf("eigen vector rows\n");
  // k = 0;
  // for( i = 0; i < lda; i++) {
  //for(j = 0; j < lda; j++) {
  //  printf("%f ", a[k]); 
  //   k++;
  // }
  //   printf("\n");
  //  }
  printf("^^^^^^^^^^^^^^^\n");
  fflush(stdout);

  delete []work;
  return info;
}

/*
 SUBROUTINE DGETRF( M,	N, A, LDA, IPIV, INFO )

      INTEGER	     INFO, LDA,	M, N

      INTEGER	     IPIV( * )

      DOUBLE	     PRECISION A( LDA, * )

PURPOSE
  DGETRF computes an LU	factorization of a general M-by-N matrix A using par-
  tial pivoting	with row interchanges.

  The factorization has	the form
     A = P * L * U
  where	P is a permutation matrix, L is	lower triangular with unit diagonal
  elements (lower trapezoidal if m > n), and U is upper	triangular (upper
  trapezoidal if m < n).
*/
extern "C" int dgetrf_(int *M, int *N, double *A, int *LDA, 
		       int * IPIV, int * INFO);

int GlyphViewer::determinant( Matrix theMatrix, double & det ) {
  int m, n, lda,info;
  int * ipiv;
  double * a;
  int i,j,k;
  
  m = theMatrix.noRows();
  n = theMatrix.noCols();

  lda = m;
  info = -1;
  det = 0.0;

  try {
    if (n < m)
      ipiv = new int[n];
    else 
      ipiv = new int[m];

    a = new double[m*n];
  }
  catch( std::bad_alloc) {
    printf("GlyphViewer::determinant, unable to allocate work arrays\n");
    exit(0);
  }
 
 // copy the matrix into 'a'
  k = 0;
  // copy the Matrix into a
  // printf("calculating determinant for matrix:\n");
  for( i = 0; i < lda; i++) {
    for(j = 0; j < lda; j++) {
      a[k] =  theMatrix(i,j); // matrix also stored in column order
      k++;
      //  printf("%f ",  theMatrix(i,j) );
    }
    //    printf("\n");
  }
  

  dgetrf_(&m, &n, a, &lda, ipiv, &info);

  k = 0;
  /*  printf("LU");
  for( i = 0; i < lda; i++) {
    for(j = 0; j < lda; j++) {
      printf("%f ", a[k]);
      k++;
    }printf("\n");
    }*/

  //The LU decomposition also makes it possible to calculate the 
  //determinant of A, which is equal to the product of the diagonal 
  //elements of the matrix U
  if (info == 0) {
    det = 1;
    for (i = 0; i < n; i++) {
      j = i*n + i;
      det = det*a[j];
    }
  }

 //do  i = 1, n
  //        if (ipvt(i) .ne. i) det(1) = -det(1)
  // printf("det before permute %f\n", det); 
  for (i = 0; i < n; i++) {
    if (ipiv[i] != (i+1)) { // fortran starts counting at 1, thus i+1
      det = -det;
      // printf("swapsign %d %d\n", ipiv[i],i+1);
    }
  }

  // printf("Determinant success %d, value %f\n", info, det);
  // fflush(stdout);
  delete []ipiv;
  delete [] a;
  return info;
}
/*
   SUBROUTINE DSYEVX(JOBZ, RANGE, UPLO, N, A, LDA, VL, VU, IL, IU,
           ABTOL, NFOUND, W, Z, LDZ, WORK, LDWORK, IWORK2, IFAIL, INFO)

*/
extern "C"  void dsyevx_(char *jobz, char *range, char* uplo, int *n,  
			 double *a, int *lda, double *vl, double *vu, 
			 double *il, double *iu,
			 double *abtol, int *m, double *w,  double  *z,
			 int *ldz, double * work, int * lwork,
			 int * iwork, int *ifail, int *info);

int  GlyphViewer::calcSymmEigensDSYEVX(  Matrix  theMatrix, int leadingDim, 
				   double * a, double * w) 
{
  char jobz = 'V'; //eigenvectors and eigenvalues
  char range = 'A'; //find all eigenvalues;
  char uplo = 'U'; // in reality we will give the whole matrix
  int lda,n,m,ldz; // dimension of matrix
  double * copy; // copy of the matrix that can get destroyed
  double * work;
  int lwork;
  double vl, vu;
  double il, iu; // limits on values to calculate, not used
  double abstol = 0.0; // let algorithm determine precision of eigenvalues

  int * iwork;
  int liwork;
  int info;
  int ifail;

  int i,k,j;

  vl = vu = il = iu = 0.0;
  
  m = lda = n = ldz = leadingDim;
  lwork = 8*n;
  liwork = 5*n;

  try {
    copy   = new double[lda*n];
    work   = new double[lwork];
    iwork  = new int[liwork];
  }
  catch( std::bad_alloc) {
  
    printf(" GlyphViewer::calcSymmEigensDSYEVR, unable to allocate arrays\n");
    return -1;
  } 


 // copy the matrix into 'a'
  k = 0;
  // copy the Matrix into a
  // printf("calculating determinant for matrix:\n");
  for( i = 0; i < lda; i++) {
    for(j = 0; j < lda; j++) {
      copy[k] =  theMatrix(i,j); // matrix also stored in column order
      a[k] = 0.0;
      k++;
    }
  }

  // JOBZ, RANGE, UPLO, N, A, LDA, VL, VU, IL, IU,
  // $                   ABSTOL, M, W, Z, LDZ, WORK, LWORK, IWORK,
  // $                   IFAIL, INFO )

  dsyevx_(&jobz, &range, & uplo, &n,  copy, 
	  &lda, &vl, &vu, &il, &iu,
               &abstol, &m, w, a,
               &ldz, work, &lwork,iwork, &ifail, &info);


  delete [] copy;
  delete [] work;
  delete [] iwork;

  return info;
}

/*
    SUBROUTINE DSYEVD( JOBZ, UPLO, N, A, LDA, W, WORK, LWORK, IWORK,
     $                   LIWORK, INFO )
*  DSYEVD computes all eigenvalues and, optionally, eigenvectors of a
*  real symmetric matrix A. If eigenvectors are desired, it uses a
*  divide and conquer algorithm.
*/
extern "C" void dsyevd_(char * JOBZ, char * UPLO, int * N, double * A, 
			int * LDA, double * W, double * WORK, int * LWORK, 
			int * IWORK, int * LIWORK, int * INFO );
int  GlyphViewer::calcSymmEigensDSYEVD(  Matrix  theMatrix, int leadingDim, 
				   double * a, double * w) 
{
  char jobz = 'V'; //eigenvectors and eigenvalues
  char uplo = 'U'; // in reality we will give the whole matrix
  int n,lda,lwork,liwork,info;
  int * iwork;
  double * work;
  int i,j,k;

  n = lda = leadingDim;
  lwork = 1 + (6*n) + (2*n*n);
  liwork = 3 + 5*n;

  try {
    work = new double[lwork];
    iwork = new int[liwork];
  }
  catch( std::bad_alloc) {
  
    printf(" GlyphViewer::calcSymmEigensDSYEVD, unable to allocate arrays\n");
    return -1;
  } 

  // copy the matrix into 'a'
  k = 0;
  // copy the Matrix into a
  for( i = 0; i < lda; i++) {
    for(j = 0; j < lda; j++) {
      a[k] =  theMatrix(i,j); // matrix also stored in column order
      k++;
    }
  }

  dsyevd_(&jobz, &uplo, &n, a, 
	      &lda, w, work, &lwork, 
			iwork, &liwork, &info );

  delete [] work;
  delete [] iwork;
 
  return info;

}


//Zhou's code
int GlyphViewer::Matrix2TensorSysR4(const Matrix& M, tensor& T)
{
  int rank = T.rank();
  if (rank != 4) {
    opserr << "NewTemplate3Dep::Matrix2TensorSysR4 - tensor must be of rank 4" << endln;
    return 1;
  }

  int nr = M.noRows();
  int nc = M.noCols();
  if (nr < 6 || nc < 6) {
    opserr << "NewTemplate3Dep::Matrix2TensorSysR4 - matrix must be no less than (6, 6)" << endln;
    return 1;
  }

  double sqrthalf = sqrt(0.5);
  double half = 0.5;
  
  // Adopt method from Helnwein (2001):  
  
  T.val(1,1,1,1) = M(0,0);
  T.val(1,1,2,2) = M(0,1);
  T.val(1,1,3,3) = M(0,2);      
  T.val(1,1,1,2) = T.val(1,1,2,1) = M(0,3) *sqrthalf;
  T.val(1,1,2,3) = T.val(1,1,3,2) = M(0,4) *sqrthalf;
  T.val(1,1,1,3) = T.val(1,1,3,1) = M(0,5) *sqrthalf;      

  T.val(2,2,1,1) = M(1,0);
  T.val(2,2,2,2) = M(1,1);
  T.val(2,2,3,3) = M(1,2);      
  T.val(2,2,1,2) = T.val(2,2,2,1) = M(1,3) *sqrthalf;
  T.val(2,2,2,3) = T.val(2,2,3,2) = M(1,4) *sqrthalf;
  T.val(2,2,1,3) = T.val(2,2,3,1) = M(1,5) *sqrthalf; 

  T.val(3,3,1,1) = M(2,0);
  T.val(3,3,2,2) = M(2,1);
  T.val(3,3,3,3) = M(2,2);      
  T.val(3,3,1,2) = T.val(3,3,2,1) = M(2,3) *sqrthalf;
  T.val(3,3,2,3) = T.val(3,3,3,2) = M(2,4) *sqrthalf;
  T.val(3,3,1,3) = T.val(3,3,3,1) = M(2,5) *sqrthalf;

  T.val(1,2,1,1) = T.val(2,1,1,1) = M(3,0) *sqrthalf;
  T.val(1,2,2,2) = T.val(2,1,2,2) = M(3,1) *sqrthalf;
  T.val(1,2,3,3) = T.val(2,1,3,3) = M(3,2) *sqrthalf;      
  T.val(1,2,1,2) = T.val(2,1,1,2) = T.val(1,2,2,1) = T.val(2,1,2,1) = M(3,3) *half;
  T.val(1,2,2,3) = T.val(2,1,2,3) = T.val(1,2,3,2) = T.val(2,1,3,2) = M(3,4) *half;
  T.val(1,2,1,3) = T.val(2,1,1,3) = T.val(1,2,3,1) = T.val(2,1,3,1) = M(3,5) *half;

  T.val(2,3,1,1) = T.val(3,2,1,1) = M(4,0) *sqrthalf;
  T.val(2,3,2,2) = T.val(3,2,2,2) = M(4,1) *sqrthalf;
  T.val(2,3,3,3) = T.val(3,2,3,3) = M(4,2) *sqrthalf;      
  T.val(2,3,1,2) = T.val(3,2,1,2) = T.val(2,3,2,1) = T.val(3,2,2,1) = M(4,3) *half;
  T.val(2,3,2,3) = T.val(3,2,2,3) = T.val(2,3,3,2) = T.val(3,2,3,2) = M(4,4) *half;
  T.val(2,3,1,3) = T.val(3,2,1,3) = T.val(2,3,3,1) = T.val(3,2,3,1) = M(4,5) *half;

  T.val(1,3,1,1) = T.val(3,1,1,1) = M(5,0) *sqrthalf;
  T.val(1,3,2,2) = T.val(3,1,2,2) = M(5,1) *sqrthalf;
  T.val(1,3,3,3) = T.val(3,1,3,3) = M(5,2) *sqrthalf;      
  T.val(1,3,1,2) = T.val(3,1,1,2) = T.val(1,3,2,1) = T.val(3,1,2,1) = M(5,3) *half;
  T.val(1,3,2,3) = T.val(3,1,2,3) = T.val(1,3,3,2) = T.val(3,1,3,2) = M(5,4) *half;
  T.val(1,3,1,3) = T.val(3,1,1,3) = T.val(1,3,3,1) = T.val(3,1,3,1) = M(5,5) *half;

  return 0;
} 



/** 
 * copied from Zhou's NewTemplate3Dep 
 * DOES NOT WORK BECAUSE OF POINTERS
*/
int GlyphViewer::Tensor2MatrixSysR4(const tensor& T, Matrix& M)
{
  int rank = T.rank();
  if (rank != 4) {
    opserr << "NewTemplate3Dep::Tensor2MatrixSysR4 - tensor must be of rank 4" << endln;
    return 1;
  }

  int nr = M.noRows();
  int nc = M.noCols();
  if (nr != 6 || nc != 6) {
    opserr << "NewTemplate3Dep::Tensor2MatrixSysR4 - matrix must be of (6, 6)" << endln;;
    return 1;
  }
  
  double sqrt2 = sqrt(2.0);
  double two = 2.0;
  
  // Adopt method from Helnwein (2001):

  M(0,0) = T.cval(1,1,1,1);
  M(0,1) = T.cval(1,1,2,2);
  M(0,2) = T.cval(1,1,3,3);      
  M(0,3) = T.cval(1,1,1,2) *sqrt2;
  M(0,4) = T.cval(1,1,2,3) *sqrt2;
  M(0,5) = T.cval(1,1,1,3) *sqrt2;      
    
  M(1,0) = T.cval(2,2,1,1);
  M(1,1) = T.cval(2,2,2,2);
  M(1,2) = T.cval(2,2,3,3);      
  M(1,3) = T.cval(2,2,1,2) *sqrt2;
  M(1,4) = T.cval(2,2,2,3) *sqrt2;
  M(1,5) = T.cval(2,2,1,3) *sqrt2;            
    
  M(2,0) = T.cval(3,3,1,1);
  M(2,1) = T.cval(3,3,2,2);
  M(2,2) = T.cval(3,3,3,3);      
  M(2,3) = T.cval(3,3,1,2) *sqrt2;
  M(2,4) = T.cval(3,3,2,3) *sqrt2;
  M(2,5) = T.cval(3,3,1,3) *sqrt2;                  
    
  M(3,0) = T.cval(1,2,1,1) *sqrt2;
  M(3,1) = T.cval(1,2,2,2) *sqrt2;
  M(3,2) = T.cval(1,2,3,3) *sqrt2;      
  M(3,3) = T.cval(1,2,1,2) *two;
  M(3,4) = T.cval(1,2,2,3) *two;
  M(3,5) = T.cval(1,2,1,3) *two;                        
    
  M(4,0) = T.cval(2,3,1,1) *sqrt2;
  M(4,1) = T.cval(2,3,2,2) *sqrt2;
  M(4,2) = T.cval(2,3,3,3) *sqrt2;      
  M(4,3) = T.cval(2,3,1,2) *two;
  M(4,4) = T.cval(2,3,2,3) *two;
  M(4,5) = T.cval(2,3,1,3) *two;                              
    
  M(5,0) = T.cval(1,3,1,1) *sqrt2;
  M(5,1) = T.cval(1,3,2,2) *sqrt2;
  M(5,2) = T.cval(1,3,3,3) *sqrt2;      
  M(5,3) = T.cval(1,3,1,2) *two;
  M(5,4) = T.cval(1,3,2,3) *two;
  M(5,5) = T.cval(1,3,1,3) *two;

    return 0;
} 

//Zhou's code - invert a tensor DOES NOT WORK BECAUSE OF TENSOR POINTERS
int GlyphViewer::Stiffness2Compliance(const tensor& S, tensor& C)
{
  int rank = 0;
  int err = 0;

  rank = S.rank();
  if (rank != 4) {
    opserr << "NewTemplate3Dep::Stiffness2Compliance - tensor must be of rank 4" << endln;
    return 1;
  }

  rank = C.rank();
  if (rank != 4) {
    opserr << "NewTemplate3Dep::Stiffness2Compliance - tensor must be of rank 4" << endln;
    return 1;
  }

  Matrix S66(6,6);
  Matrix C66(6,6);

  err += Tensor2MatrixSysR4(S, S66);
  err += S66.Invert(C66);
  err += Matrix2TensorSysR4(C66, C);

  return err;
}



/**
 * Converts eigentensors coming out of the solver into a 
 * form usable by current methods
 * Assumption: each eig eigenTensor is [ xx yy zz xy yz xz] and
 * the 9 is [xx xy xz yx yy yz zx zy zz ]. So the indeces
 * are hardcoded
 */
void GlyphViewer::sixVectorToSymMatrix( double eigTensors[36] ) {

  int i;
  int base;// base index for eigTensor array
  printf("GlyphViewer::sixVectorToSymMatrix\n");
  for (i = 0; i < 6; i++){
    base = i*6;
    f_matrices[i][0] = eigTensors[base];
    f_matrices[i][4] = eigTensors[base + 1];
    f_matrices[i][8] = eigTensors[base + 2];
    f_matrices[i][1] =  f_matrices[i][3] = eigTensors[base + 3];
    f_matrices[i][5] =  f_matrices[i][7] = eigTensors[base + 4];
    f_matrices[i][2] =  f_matrices[i][6] = eigTensors[base + 5];
  }

  /*  for(i = 0; i < 36; i++ ) {
   
      printf("%f ", eigTensors[i]);
    
 
      }*/
}




/**
 * Rebecca Brannon's method to get stretch eigentensors 
 * Start with 6 x 6 Helnwein matrix 
 * @param diffStiff: difference between current stiffness and elastic stiffness
 * @returns 0 for correct execution 
 */
int GlyphViewer::highamPolar( Matrix diffStiff, Matrix & stretch, Matrix & rotate ) {
  // iterative: D = RU, R is rotation, U is % stretch
  // R_0 = D, D = diffstiff
  int i,j,k,l;
  Matrix Rcurr(6,6);
  Matrix RcurrInvert(6,6);
  Matrix RcurrInverseTranspose(6,6);
  Matrix Rnext(6,6);
  Matrix testDiff(6,6); // how much values are changing on this iteration
  double frobenius = 0;
  
  Rcurr = diffStiff;
  transposeMatrix( Rcurr, testDiff );
  testDiff = testDiff*Rcurr;
  // printf("Start Identity Test WITH ORIGINAL DIFFERENCE OF STIFF\n");
  // for (k = 0; k < 6; k++ ) {
  // for (l = 0; l < 6; l++ ) {
  //  printf("%f ", testDiff(k,l) );
  // }
  // printf("\n");
  //  }
  // printf("End Identity TestWITH ORIGINAL DIFFERENCE OF STIFF\n");
  
  
  for (i = 0; i < 10; i++) { // Newton Iterations

    if ( diffStiff.Invert(RcurrInvert ) != 0) { // already done by Frank.!
      printf("GlyphViewer::highamPolar, unable to invert matrix\n");
      return -1;
    }
    
    transposeMatrix( RcurrInvert,  RcurrInverseTranspose );
    
    // R_k+1 = 1/2(R_k + R_k^(-1T) ), where ^(-1T) means inverse transpose
    Rnext = (Rcurr + RcurrInverseTranspose);
    Rnext = 0.5*Rnext;
    
    //------------ DEBUG TEST ---------------------
    testDiff = diffStiff - Rnext;
    
    
    double test = 0;
    // check how we're doing by the Frobenius norm of the diff from the new iter and
    // original
    for (k = 0; k < 6; k++ )
      for (l = 0; l < 6; l++ )
	test = test + (fabs(testDiff(k,l)));

    printf("highamPolar, iteration %d change %f,\n", i, test);
    // for (k = 0; k < 6; k++ ) {
    // for (l = 0; l < 6; l++ ) {
    //	printf("%f ", Rnext(k,l) );
    // }
    // printf("\n");
    //}
  //-------------- END TEST

    Rcurr = Rnext;

  } // Newton iterations

  rotate = Rcurr;

  // Test: does Rcurr*Rcurr^T = I?
 
  transposeMatrix( Rcurr, testDiff );
  testDiff = testDiff*Rcurr;
  // printf("Start Identity Test\n");
  //for (k = 0; k < 6; k++ ) {
  // for (l = 0; l < 6; l++ ) {
  //  printf("%f ", testDiff(k,l) );
  // }
  // printf("\n");
  // }

  // printf("End Identity Test\n");

  // U = R^-1 D
  if ( Rcurr.Invert(RcurrInvert) != 0 )  {
    printf("GlyphViewer::highamPolar, unable to invert matrix\n");
      return -1;
  } 
  
  stretch = RcurrInvert*diffStiff;
  return 0;
 
}



/*
 * @pre: both matrices have the same dimensions */
void GlyphViewer::transposeMatrix( Matrix source, Matrix & transpose ) {
  int i,j,rows,cols;
  rows = source.noRows();
  cols = source.noCols();

  for (i = 0; i < rows; i++ )
    for (j = 0; j < cols; j++)
      transpose(i,j) = source(j,i);

}

/**
 * Symmetric part of a matrix is 1/2(T + T^T)
 */
void GlyphViewer::symmetrizeMatrix(Matrix source, Matrix & symm ) {
  int i,j,rows,cols;


  rows = source.noRows();
  cols = source.noCols();

  for (i = 0; i < rows; i++ )
    for (j = 0; j < cols; j++)
      symm(i,j) = 0.5*( source(i,j) + source(j,i) );

}


/**
 * Find the difference of the Frobenius norm of 2 matrices */
double GlyphViewer::differenceOfNorm(Matrix a, Matrix b) {
  double norm;
  int i,j;
  norm = 0;
  Matrix diff(6,6);

  /* for ( i = 0; i < a.noRows(); i++)
    for (j = 0; j < a.noCols(); j++)
    diff(i,j) = a(i,j) - b(i,j);*/
  diff = a-b;
  // printf("a - b\n");
  //opserr << diff;
  // now, the frobenius norm of diff

  for ( i = 0; i < a.noRows(); i++) {
    for (j = 0; j < a.noCols(); j++) {
      if (diff(i,j) != 0 )  // try to avoid overflow garbage
	norm += (diff(i,j)*diff(i,j));
    }
  }
  printf("norm %f\n", norm);
  norm = sqrt( norm );
  printf("sqrt %f\n", norm);

  return norm;
}




/**
 * contract the current stiffness aigainst the
 * eigentensor with the smallest eigenvalue
 * @pre! eigentensors have been calculated
 */
void GlyphViewer::testSmallest( Tensor stiffness) {
  int i;
  int pdim[2] = {3,3};
  int pdim4[4] = {3,3,3,3};
  
  int smallest = 0;
  double smallVal =  f_real_eigens[0];
  for (i = 1; i < 6; i++) {
    if (f_real_eigens[i] < smallVal) {
      smallest = i;
      smallVal = f_real_eigens[i];
    }
  }
  
  // ok, now we can use the eigenvector row-wise because 
  // it is symmetric and equals its own transpose
  //BJtensor (int rank_of_BJtensor, const int *pdim, double *values)
  Tensor eigen(2, pdim, f_matrices[smallest]);
  printf("eigentensor\n");
  eigen.print(); 

  Matrix m_stiff(6,6);
  Matrix m_compliance(6,6);
  Tensor compliance(4, pdim4, 0.0 );
  
  Tensor2MatrixSysR4( stiffness, m_stiff );
  m_stiff.Invert( m_compliance );

  Matrix2TensorSysR4( m_compliance, compliance );

  Tensor result = compliance("ijkl")*eigen("kl");

  printf("Frobenius norm of contracted eigentensor: %f\n", result.Frobenius_norm());

 
  // try testing against standard orthonormal directions.
  for (i = 0; i < 6; i++) {
    Tensor  x( 2, pdim, orthoEigens[i] );
    result = compliance("ijkl")*x("kl");
    printf("Frobenius norm of contracted orthoEigen %d: %f\n", i, result.Frobenius_norm());
  }
}
  
/**
 * check two doubles for equality */ 
int GlyphViewer::sameDouble( double t1, double t2) {
  
  double bigger;
  double smaller;
  double nearZero = 0.000001;

  if ( (fabs(t1) < nearZero) && (fabs(t2) < nearZero) )
    return 1; // they are both near zero

  if ( t1 > 0 ) { // positive value: establish squeeze
     bigger =  t1*1.01;
     smaller = t1*0.99;	      
  }
  else {
      bigger =  t1*0.99;
      smaller =  t1*1.01;
   }
	  
   if ( t2 > bigger || t2 < smaller ){
     return 0;
   }

	
  return 1; // they are the same within  0.01
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
double  GlyphViewer::orient3D( double a[3], double b[3], double c[3], double d[3] )

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

/*
* function [R,U,T,J,SQRTJ,U0,F,Gorig,G,C,zeroEigs0] = polarDecompA(D, toler)
* function [R,U,T,J,SQRTJ,U0,F,Gorig,G,C,zeroEigs0] = polarDecompA(D, toler)
* D is a square matrix with det(D) >= 0 (not checked, see below).
* toler controls near-zero threshold.
*   Plausible values of toler range from 1 to 10000, but 0 and sqrt(1/eps)
*   are also worth experimenting with.
*   The actual cutoff for an eigenvalue of J = D^T D is toler * eps * J(n,n),
*   where J(n,n) is the maximum eigenvalue.
*
* Returns
*    R orthogonal rotation, with det(R) = +1.
*    U symmetric positive semidefinite stretch, D = R U is the polar decomp
*    T eigenvectors of U (actually of U0); column-1 is eigenvector of 0 for
*        D IF D IS SINGULAR.
*    Other diagnostic values. zeroEigs0 is how many eigenvalues were at
*    or below the cutoff.
*
* Adjustment for det(D) < 0:  Define T1 = T with column-1 negated.
* Define R1 = R * T * T1'.  Now det(R1) < 0 and D ~= R1 * U.
* In theory T1 can be T with any odd number of columns negated.
* Is it true that R1 comes out the same in all cases ????
* The answer is unclear, so this is not implemented, and we require
* det(D) >= 0 to guarantee a correct decomposition.
* A good value for toler 10000
* @pre! we have already checked for bigger negative eigenvalues (at best,
* only very small negative eigenvalues)
* @param  D:curr stiffness , rotate, stretch, tolerance for zero eigenvalue
*/

int GlyphViewer::vanGelderPolar( Matrix D, Matrix &R, Matrix &U, double toler)
{
  double * T; // eigenvectors
  double * J; // jacobean like eigenvalues
  double * wi;
  int smallest,i,j,k, success, zeroEigs;
  double det;
  double eps =  2.2204e-16; // matlab's epsilon
  int n = D.noRows();

  Matrix DTD(n,n);
  Matrix T_m(n,n);
  Matrix sqrtJ(n,n);
  Matrix trans(n,n);
  Matrix F(n,n);
  Matrix G(n,n);
  Vector G_k(n);
  Vector temp(n);
  Vector proj(n);
  success = 0;

 if (n != D.noCols()) {
   printf("GlyphViewer::vanGelderPolar. Non-square matrix, cannot compute\n");
   return -1;
 }


 //DTD = (D')*D;
  transposeMatrix( D, DTD );

  DTD = DTD*D; // = stretch^2 non-symmetric if D nonsymmetric


 try {
   T = new double[n*n];
   J = new double[n];
   wi = new double[n];
 }
 catch( std::bad_alloc) {
    printf("GlyphViewer::vanGelderPolar, Unable to allocate T &J arrays\n");
    exit(0);
 }

 //[T,J] = eig(DTD); stretch = T J T'
 success = calcSymmEigens(  DTD, n, T, J ); 


 if ( success != 0 ) {
   printf("GlyphViewer::vanGelderPolar: eigenvector calculation unsuccessful,\n");
   printf("Lapack  dsyev returns %d\n", success);
   delete [] T;
   delete [] J;
   return success;
 }

 //if J(n,n) == 0, ALL_ZERO_EIGS
 if (J[n-1] == 0.0) {
   printf("GlyphViewer::vanGelderPolar: WARNING, ALL_ZERO_EIGS");
   delete [] T;
   delete [] J;
   return -1;
 }




  k = 0;
  // copy the array into the matrix
  for( i = 0; i < n; i++) {
    for(j = 0; j < n; j++) {
      T_m(j,i) = T[k]; // matrix also stored in column order
      k++;
    }
  }
 
  printf("********** T ******\n");
  opserr << T_m;

  printf("****** J *************\n");
  
  for( i = 0; i < n; i++) 
    printf("%f ", J[i]);
  printf("\n*********************\n");

  j = determinant( T_m, det );

  if (j != 0) {
    printf("GlyphViewer::vanGelderPolar, 1st determinant calculation failed\n");   delete [] T;
   delete [] J;
   return -1;

  }
  
  //if det(T) < 0
  //  T(:,1) = -T(:,1);
  if (j == 0 && det < 0 ) {
    printf("Swapping sign of first, negative detirminant\n");
    for (i = 0; i < n; i++) {
      T[i] = -T[i];
      T_m(i,0) = -T_m(i,0); 
    }
  }

 printf("********** T ******\n");
  opserr << T_m;

  // precondition: sorted eigenvalues, which lapack gives us
  zeroEigs = 0;
  for (i = 0; i < n; i++ ) {
    if (J[i] < 0) {
      J[i]= 0;
    }
    printf("J[%d]: %f, ( toler * eps * J[n-1]):%f\n",
	   i,J[i],( toler * eps * J[n-1])  );
    if (J[i] <= ( toler * eps * J[n-1]) ){
      zeroEigs = zeroEigs + 1;
      J[i]= 0;
    }
  }

  printf("ZeroEigs: %d, should be 1\n", zeroEigs); 
  // Later code cannot handle zeroEigs > 1.
  if  (zeroEigs > 1) {
      return -zeroEigs;
  }
    // set up sqrtJ matrix
    for(i =0; i < n; i++) {
      for ( j =0; j < n; j++) {
	if (i == j) {
	  sqrtJ(i,j) = sqrt(J[i]);
	}
	else 
	  sqrtJ(i,j) = 0.0;
      }
    } 
    
    printf("SQRTJ = \n");
     for (i = 0; i < n; i++) {
	for (j = 0; j < n; j++) {
	  if (i == j)
	  printf("%f ", sqrtJ(i,j));
	}printf("\n");
	}

    if (zeroEigs == 0 ) { // non-singular U, R = DU^(-1)


      //U0 = (T)*SQRTJ*(T');
   
      transposeMatrix( T_m, trans );
      //  printf("*********T**********\n");
      //	opserr << T_m;

      U = T_m*sqrtJ*trans;

      printf( "U = T*sqrtJ*T_trans\n");
      opserr << U;

      // Ensure U0 is symmetric
      //U0 = (U0 + U0') / 2;    
      transposeMatrix(U,trans);
      U = 0.5 * (U + trans);

      printf("U tuned up\n");
      opserr << U;
      //R = DU^(-1)

     for (i = 0; i < n; i++) {
	for (j = 0; j < n; j++) {
	  DTD(i,j) = 0.0;
	}
     }

     j =  invertMatrix( U, DTD );
     //      j = U.Invert(DTD);
      printf("U inverse\n");
      opserr << DTD;

      

      if ( j != 0 ) {
	printf("GlyphViewer::vanGelderPolar, error inverting U\n");
	 success -1;
      }
      else {
	R = D*DTD;  //  D*U^(-1)
	success = 0;
      }
      /*
      printf("Simple Vangelder stretch (U)\n");
      for (i = 0; i < n; i++) {
	for (j = 0; j < n; j++) {
	  printf("%f ", U(i,j));
	}printf("\n");
      }
      */
     
     printf("Simple Vangelder rotate (R)\n");
      for (i = 0; i < n; i++) {
	for (j = 0; j < n; j++) {
	  printf("%f ", R(i,j));
	}printf("\n");	}

      //  calcSymmEigens(  U, n, T, J ); 
	printf("++++++++++++\nRotation matrix, no near-zero eigen\n");
	opserr << R;
	printf("++++++++++++\n");
	printf("STRETCH \n");
	  opserr << U;
    }
    else { // case 1 zero eigenvalue
      // C = RT, C unknown, ultimately we will solve for R orthogonal
      //calculate F, the fake inverse of sqrt(J)
 
   
      // case of one or more zero eigenvalue
      for (j = 0; j < zeroEigs; j++) {
	if (sqrtJ(j,j) == 0) 
	  sqrtJ(j,j) = 1.0; 
      }
    
      j = sqrtJ.Invert(F); // check return value!
      if ( j != 0 ) {
	printf("GlyphViewer::vanGelderPolar, error inverting F\n");
	delete [] T;
	delete [] J;
	return -1;
      }
      /*
      printf("F after inverting\n");
      for (i = 0; i < n; i++) {
	for (j = 0; j < n; j++){
	  printf("%f ", F(i,j) );
	}
	}*/
     // DTF = CF
      G = (D*T_m)*F;
      /*
      printf("Gorig\n");
      for (i = 0; i < n; i++) {
	for (j = 0; j < n; j++){
	  printf("%f ", G(i,j) );
	}printf("\n");
	}*/
// Tune up G.  All columns of G not corresponding to a zero eigenvalue
// should be orthogonal.
// All columns of G not corresponding to a zero eigenvalue
// should be unit length.  Force this to be true.
//
//for j = n : -1 : zeroEigs+1
//    temp = G(:,j);
//    for k = j+1:n
//        temp = temp - orthoprojection(temp, G(:,k));
//        end;
//    G(:,j) = normc(temp);
//    end;
      
      
      for( i = n-1; i > (zeroEigs -1); i-- ) { // columns to orthonorm

	// get current column of G
	for ( j = 0; j < n; j++ ) {
	  temp(j) = G(j,i); // G is organized column-wise, use J for temp
	}
	
	// now iterate through previous columns
	for ( j = i+1; j < n; j++ ) {

	  for ( k = 0; k < n; k++ ) { // copy the column from G
	    G_k(k) = G(k,j);
	  }

	  // now project
	  orthoprojection(temp, G_k, proj);
	  temp -= proj;
	}
	// normalize and copy back into G
	temp.Normalize();
	for ( j = 0; j < n; j++ ) {
	  G(j,i) = temp(j); 
	}
	
      }// non-zero columns to orthonorm

  
// find ROW with smallest magnitude in G for 1
//X = zeros(n,1);
//for i = 1:n
//    X(i) = norm(G(i,:));
//    end;
      for (i = 0; i < n; i++) {
	temp(i) = 0.0;
      }
      for (i = 0; i < n; i++) {
	for (j = 0; j < n; j++){
	  temp(i) += G(i,j);
	}
      }
      
      

// sum the first row of G - this is smallest
      det = temp(0);
      j = 0;
      // identify smallest row of G
      for( i = 1; i < n; i++ ) { 
	if ( temp(i) < det ) {
	  det = temp(i);
	  j = i;
	}
      }
      
      // create "guess" vector for 1st column of G
      for( i = 0; i < n; i++ ) { 
	if ( i == j )
	  temp(i) = 1.0;
	else 
	  temp(i) = 0.0;
      }
      
      // now see if the determinant of G with this new vector 
      // will be negative.
       for( i = 0; i < n; i++ ) { 
	 G(i,0) = temp(i);
       }
       /*   printf("G_Later_later\n");
      for (i = 0; i < n; i++) {
	for (k = 0; k < n; k++){
	  printf("%f ", G(i,k) );
	}printf("\n");
	}*/    

      i = determinant( G, det );
      if ( det < 0 ) {
	temp(j) = -temp(j); // j is still guess
      }

      // Use Gram-Schmidt
      //Only coded for zeroEigs = 0 or 1, but generalizes.
      //P1 = zeros(n,1);
      //for j = 2:n
      //P1 = P1 - orthoprojection(G1, G(:,j));
      // end;
     
      for( i = n-1; i > (zeroEigs -1); i-- ) { // columns to orthonorm

	// get current column of G
	for ( j = 0; j < n; j++ ) {
	  G_k(j) = G(j,i); // G is organized column-wise
	}
	// now project
	orthoprojection(temp, G_k, proj);
	temp -= proj;
	}
	// normalize and copy back into G
	temp.Normalize();

	for ( j = 0; j < n; j++ ) {
	  G(j,0) = temp(j);

	}
	printf("\nC =\n");
	for (i = 0; i <n; i++) {
	  for ( j = 0; j < n; j++ ) {
	    printf("%f ", G(i,j));
	  }printf("\n");
	}
	//C = [C1 G(:, 2:n)];
	//R = C*(T');
	//U2 = inv(R) * D;   % more accurate than assuming inv(R) = R'.
	//U = (U2 + U2') / 2;	  
	transposeMatrix( T_m, trans );
	  
	R = G*trans;
	R.Invert(trans);
	
	printf("++++++++++++\nRotation matrix\n");
	opserr << R;
	printf("++++++++++++\n");

	U = trans*D;
	transposeMatrix( U, trans );
	U = ( U + trans)*0.5;

	printf("STRETCH \n");
	  opserr << U;

	/*	printf("\nR =\n");
	for (i = 0; i <n; i++) {
	  for ( j = 0; j < n; j++ ) {
	    printf("%f ", R(i,j));
	  }printf("\n");
	}

	printf("\nU =\n");
	for (i = 0; i <n; i++) {
	  for ( j = 0; j < n; j++ ) {
	    printf("%f ", U(i,j));
	  }printf("\n");
	  }*/

    } // end of case 1 zero eigenvalue

    

    delete [] T;
    delete [] J;
    delete [] wi;
    return success;
}

/*function [T] = orthoprojection(V,W)
* function [T] = orthoprojection(V,W) V is column vector, 
* W is a column vector 
* T is the result of projecting V against W ;
* you would subtract this whole thing from V to get the actual projection
*/

void GlyphViewer::orthoprojection(Vector V, Vector W, Vector& T ) 
{ 
  int i; 
  double magnitude, denom; 
 
  
 //T =   (( (W') * V ) / ( (W')*W ))*W;

  magnitude = denom = 0.0;
  for( i = 0; i < V.Size(); i++ ) {
    magnitude = magnitude + (W[i]*V[i]);
    denom = denom + (W[i]*W[i]);
  }
  
  magnitude = magnitude/denom;

  for( i = 0; i < V.Size(); i++ )
    T[i] = W[i]*magnitude;
}

/**
 * check the symmetric part of the tensor for zero or
 * negative eigenvalues
 */
void GlyphViewer::testBifurcation( Matrix D ){

 int i,j,rows,cols, count;
 double * a;
 double * w;
 double * wi;
 double eps =  2.2204e-16*10000; //matlab's epsilon plus tolerance

  rows = D.noRows();
  cols = D.noCols();

  /*Matrix  symm(rows, cols);
  for (i = 0; i < rows; i++ ) {
    for (j = 0; j < cols; j++) {
      symm(i,j) = 0.5*( D(i,j) + D(j,i) );
    }
    }*/

  try{ 
    a = new double[rows*cols];
    w = new double[rows];
    wi = new double[rows];
  }
  catch( std::bad_alloc) {
    printf("ReynoldsGlyph::testBifurcation, Unable to allocate eigens\n");
    exit(0);
  }
  

  j = calcAsymmEigens(  D, rows, a, w, wi);
  // NOTE: eigenvalues are sorted

  if (j != 0 )
    printf("ReynoldsGlyph::testBifurcation, eigen test failed\n");

  nonPositiveDefinite = 0;
  indeterminate = 0;
  count = 0;
  printf("\n ReynoldsGlyph::testbifurcation\n");
  for (i = 0; i < 6; i++ ) {
    printf("%f ", w[i]);
    if (fabs(w[i]) < eps ) // we consider this close to zero
      count++;

    if (w[i] < -0.5) 
      indeterminate = 1;

  }
  if (count > 0 )
    nonPositiveDefinite = 1;
  if (count > 1)
    indeterminate = 1;
  printf("\n indeterminate %d count %d\n\n",
	 indeterminate,count);

  delete [] a;
  delete [] w;
  delete [] wi;
} 


/**
 * Eigentensors are not unique - they can be multiplied by
 * any constant (including -1) and still fulfill E.x = lambda.x
 * Test whether the first eigentensor or its negative
 * produces greater strain
 * @param  Tensor curr is current stiffness 
 */
/*
void GlyphViewer::setEigenDirection(Tensor curr){
  double * data;
  int pdim[] = {3,3};
  int i;
  try {
    data = new double[9];
  }
  catch( bad_alloc &ba) {
    printf("GlyphViewer::setEigenDirection, Unable to allocate data array\n");
  }

  for(i = 0; i < 8; i++)
    data[i] = f_matrices[0][i];

  //Tensor curr is current stiffness 
  // contract with first eigen tensor, and 
  // get frob. norm of strain increment
  // create BJtensor as copy of f_matrices[0]
  tensor firstEigen = new tensor(2, pdim, data);

  tensor result;
  result = curr("ijkl")*firstEigen("kl");
  // then negate first eigentensor, and get frob. norm. 

  // compare norms. If negated is bigger, negate
  // the first eigentensor in f_matrices. 
}
*/

void GlyphViewer::eigenMode() {
  double * a;
  double * w;
  Matrix theMatrix(3,3);
  int i,j, count;
  double test, iso, aniso, red, green, blue;
  int min, max;
  double color[3];
  double eps =  2.2204e-16*10000; // matlab's epsilon*10000

  double w_abs[3];

  // check f_real_eigens for 5 alike -- isotropy!!
  test= f_real_eigens[0];
  count = 1;

  for ( i = 1; i < 6; i++ ) {
    if (sameDouble(test, f_real_eigens[i]))
	count++;
  }
    if ( count == 5 ) {
          printf("######## Eigen mode ISOTROPY!!! check min eigen for hardening#######\n");
      color[0] = 1;
      color[1] = color[2] = 0;
      return;
    } 
    
    count = 1;
    test = f_real_eigens[5];
    for ( i = 0; i < 5; i++ ) {
      if (sameDouble(test, f_real_eigens[i]))
	count++;
    }
    if ( count == 5 ) {
        printf("######## 2nd case Eigen mode ISOTROPY!!! check min eigen for hardening#######\n");
      color[0] = 1;
      color[1] = color[2] = 0;
      return;
    } 

  // otherwise do eigendecomposition on first matrix
    // copy the array into the matrix
  try{ 
    a = new double[9];
    w = new double[3];
  
  }
  catch( std::bad_alloc) {
    printf("ReynoldsGlyph::eigenMode, Unable to allocate eigens\n");
    color[0] = color[1] = color[2] = 0;
    return;
  }
  
  //[xx xy xz yx yy yz zx zy zz ], symmetric
  for ( i = 0; i < 3; i++) {
    for( j = 0; j < 3; j++ ) {
      theMatrix(i,j) = f_matrices[0][i*3+j];
    }
  }

  count = calcSymmEigens(  theMatrix, 3, a, w);   

  if (count != 0) {
      printf("ReynoldsGlyph::eigenMode, Unable to calculate eigens\n");
    color[0] = color[1] = color[2] = 0;
  }
  else {
    printf("Eigentensor:\n");
    opserr << theMatrix;
    printf("Low eigentensor egenvalues %f %f %f \n", w[0], w[1], w[2] );
    // check for isotropic part
    iso = (w[0] + w[1] + w[2])/3;
    if ( fabs(iso) > eps*10000 ) {
      printf("####### EIGENMODE isotropic part of %f \n", iso );
      
      // remove isotropy
      for (i = 0; i < 3; i++ )
	w[i] = w[i] - iso;
      
    }
    // now, let's get the anisotropic part
    // find abs(min) and abs(max), divide min by max (range will be 0 .. 0.5)
    
    for (i = 0; i < 3; i++ )
      w_abs[i] = fabs(w[i]);
    
    
    min = max = 0;
    test = w_abs[0];
    for (i = 1; i < 3; i++ ) {
      if (w_abs[i] < test) {
	test = w_abs[i];
	min = i;
      }
    }
    
    test = w_abs[0];
    for (i = 1; i < 3; i++ ) {
      if (w_abs[i] > test) {
	test = w_abs[i];
	max = i;
      }
    }
   
    // 0 - 0.5 mapped to 0 - 1
    aniso = (w_abs[min]/w_abs[max])*2;

    red = fabs(iso);
    
    // aniso is between 0 and 1
    blue = aniso;
    green = (1.0 - aniso); // if aniso is 0, shear is 1

    // remove isotropy factor
    if ( red > eps ) {
      blue = blue * ( 1.0 - red );
      green = green * ( 1.0 - red );
    } 
    
    printf("iso %f min %f max %f min/max %f\n", iso, w_abs[min], w_abs[max], 
 (w_abs[min]/w_abs[max]));
    printf("########EIGENMODE Shear vs. Axisymmetric: %f , where 0 is shear and 0.5 is axis\n",
	 (w_abs[min]/w_abs[max]) );
    printf("red(iso) %f green(shear) %f blue(axi) %f\n", red, green, blue);

    // 1.0 - iso = aniso
    color[0] = red;
    color[1] = green;
    color[2] = blue;
  
    lodeAngle( theMatrix, red, blue);

    printf("Lode coords: radius %f theta %f (-30 TXC, 0 Shear, +30 TXE)\n",
	   red,radiansToDegrees( blue) );
  } // able to calculate eigens
  
  
  delete[] a;
  delete[] w;
  
}


/**
 * get yield surface normal & plastic flow normal for this integration point */
/*void GlyphViewer::getMaterialNormals(Matrix & yieldNormal, Matrix & flowNormal)
{
  // NewTemplate3Dep

  // Template3Dep  

}*/

/** go through f_real_eigens, count similar eigenvalues within 0.5
 * if 5 the same, it's spherical
 * @returns true for spherical */
bool GlyphViewer::isSpherical() {
  int i,count1;

  count1 = 1;
  for (i = 1; i < 6; i++) {
    if (sameDouble(f_real_eigens[0],f_real_eigens[i])){

	count1++;
    }

  }
  if (count1 == 5)
    return true;
  
  count1 = 1;
  for (i = 0; i < 5; i++) {
    if (sameDouble(f_real_eigens[5],f_real_eigens[i])){

	count1++;
    }
  }
  if (count1 == 5)
    return true;

  return false;
}

void GlyphViewer::eigenMode2( double iso, double txc, double txe, double shear)
{
  double * a;
  double * w;

  Matrix theMatrix(3,3);
  int i,j, count;
  double test, aniso;
  int min, max;
  double w_abs[3]; // absolute value of eigenvalues from eigentensor

  int dim[4] = {3,3,3,3};
  double * eigCopy;

  double eps =  2.2204e-16*10000; // matlab's epsilon*10000


  // check f_real_eigens for 5 alike -- isotropy!!
  txc = txe = shear = iso = 0;
  if ( isSpherical() ) {
    iso = 1;
    return;
  } 
  


  // otherwise do eigendecomposition on first eigentensor
  try{ 
    a = new double[9];
    w = new double[3];
    
  }
  catch( std::bad_alloc) {
    printf("ReynoldsGlyph::eigenMode2, Unable to allocate eigens\n");
    exit(-1);
  }
  
  //[xx xy xz yx yy yz zx zy zz ], symmetric, copy into matrix
  for ( i = 0; i < 3; i++) {
    for( j = 0; j < 3; j++ ) {
      theMatrix(i,j) = f_matrices[0][i*3+j];
    }
  }

  // Test for match with Lode angle
  // YieldSurfNormal_ij Elastic_ijkl Eigentensor_kl > 0
 


  if (  theElement != NULL && 
	theElement->getClassTag() == ELE_TAG_EightNodeBrick 
	&& 	intgrPoint > -1 
	&& 	intgrPoint < 8	) {
    
    Tensor elastic( 4,dim, 0.0 );
    Tensor yieldNorm(2,dim,0.0);
    

    EightNodeBrick * el = (EightNodeBrick *) ( theElement );
      
    elastic =  el->ElasticStiffnessTensor( intgrPoint );
    
    yieldNorm = el-> dFods(intgrPoint );

    // finally, put eigentensor into a BJtensor
    // allocate and 
    try {
      eigCopy = new double[9];
    }
    catch( std::bad_alloc) {
      printf("ReynoldsGlyph::eigenMode2, unable to allocate data\n");
      exit(-1);
    }
    for (i = 0; i < 9; i++) {
      eigCopy[i] =  f_matrices[0][i];
    }
    stresstensor t_eig(eigCopy); // destructor will de-allocate
   
    t_eig.print("EIG", "low eigentensorm");
    test = testOrientation(yieldNorm, elastic, t_eig);
    printf(" testOrientation returns %f\n", test);
    printf("%%%%%%%%%%%%%%%%%%%%%%\n");
  
    if ( test < 0 ) {
      for ( i = 0; i < 3; i++) {
	for( j = 0; j < 3; j++ ) {
	  theMatrix(i,j) = theMatrix(i,j)*(-1);
	}
      }
      for (i = 0; i < 9; i++) {
       f_matrices[0][i] =  (-1)*f_matrices[0][i];
      }
    }
    count = calcSymmEigens(  theMatrix, 3, a, w); 
    if (count != 0) {
      printf("ReynoldsGlyph::eigenMode2, Unable to calculate eigens\n");
      return;
    }
   

    printf("eigenvalues %f %f %f\n", w[0], w[1], w[2]);
    // check for isotropic part
    iso = (w[0] + w[1] + w[2])/3; 
    if ( fabs(iso) > eps ) {
      
      for (i = 0; i < 3; i++ )
	w[i] = w[i] - iso;
    } // remove isotropy
    
    
    for (i = 0; i < 3; i++ )
      w_abs[i] = fabs(w[i]);
    
    //get shear percent
    min = max = 0;
    test = w_abs[0];
    for (i = 1; i < 3; i++ ) {
      if (w_abs[i] < test) {
	test = w_abs[i];
	min = i;
      }
    }
    
    test = w_abs[0];
    for (i = 1; i < 3; i++ ) {
      if (w_abs[i] > test) {
	test = w_abs[i];
	max = i;
      }
    }
   
    // 0 .. 0.5 mapped to 0 .. 1
    aniso = (w_abs[min]/w_abs[max])*2;
    
    shear = (1.0 - aniso); // if aniso is 0, shear is 1
    
    //test txe vs txc if there is any.
    if (aniso >  eps ) {
      // if max is positive, triaxial extension
      if ( w[max] > 0 )
	txe = aniso;
      //else triaxial compression
      else
	txc = aniso;
      
    }
    printf("iso:%f shear:%f txc:%f txe:%f \n%%%%%%%%%%%%%%%%%\n", iso, shear, txc, txe);
    
  } // if EightNodeBrick
  /////////////

 

 
  
  delete[] a;
  delete[] w;
  

}

/**
 * Test an eigentensor's orientation (both t and -1*t are valid eigentensors)
 ** -- match with Lode angle
 * YieldSurfNormal_ij Elastic_ijkl Eigentensor_kl > 0
 */
double GlyphViewer::testOrientation(BJtensor  yieldNorm,  BJtensor elastic, BJtensor eig) {
  int i,j,k,m;
  double sum = 0;
  double curr = 0;
  for (i = 1; i < 4; i++ ) {
    for (j = 1; j < 4; j++) {
      for (k = 1; k < 4; k++ ) {
	for (m = 1; m < 4; m++) {
	
	  curr = ((nDarray)yieldNorm).cval(i,j)*((nDarray)elastic).val4(i,j,k,m)*((nDarray)eig).cval(k,m);
	 
	  sum +=curr;
 

	}
      }
    }
  }
  return sum;

}


extern "C" int dgetrf_( int * M, int * N, double * A, 
int * LDA, int * IPIV, int * INFO );

extern "C" int dgetri_(  int * N,  double * A, int * LDA, int * IPIV, double * WORK, int * LWORK, int * INFO );
/**
 * @returns 0 for correct execution */
int GlyphViewer::invertMatrix( Matrix source, Matrix & inverse ) {

  int n, lda, lwork, info,i,j,k;
  double * a;
  int * ipiv;
  double * work;
  n = lda = source.noCols();

  a = new double[n*n];
  ipiv = new int[n];

  lwork = 10*n; //no idea what the optimal block size is..
  work = new double[lwork];
  info = -1;

  if ( a == NULL || work == NULL || ipiv == NULL ) {
    printf(" GlyphViewer::invertMatrix, unable to allocate arrays\n");
    return info;
  } 

 // copy the Matrix into 'a'
  info = k = 0;
  for( i = 0; i < lda; i++) {
    for(j = 0; j < lda; j++) {
      a[k] = source(i,j); // matrix also stored in column order
      k++;
    }
  }
  
  // ok, let's make the call
  dgetrf_( &n, &n, a, &n, ipiv, &info);  // get A into L-U factorized form
 
  if ( info != 0 ) {
    printf("GlyphViewer::invertMatrix, LU factorization did not succeed\n");
    return info;
  }
  dgetri_(  &n,  a, &lda, ipiv, work, &lwork, &info ); // then get inverse

  if ( info != 0 ) {
    printf("GlyphViewer::invertMatrix, inversion did not succeed\n");
    return info;
  }
  // copy into matrix
  Matrix m(a,6,6);
  opserr << "m\n" << m;
  for (i = 0; i < 6; i++)
    for (j = 0; j < 6; j++)
   inverse(i,j) = m(i,j); // deep copy so we can clean up? no, not working. do a copy.
   
   //delete []a;
   delete []ipiv;
   delete []work;
   return info;

}


/**
 * color selected eigentensor by its lode angle 
 * @pre if first eiegntensor, it has correct direction
 */
void GlyphViewer::eigenColor(int petal, float color[4]) {
  Matrix theMatrix(3,3);
  int i,j;
  double test, radius,theta;
  //  int min, max;
  // double w_abs[3]; // absolute value of eigenvalues from eigentensor

  //  int dim[4] = {3,3,3,3};
  // double * eigCopy;

  double eps =  2.2204e-16*10000; // matlab's epsilon*10000
  //  double piOverTwo = 1.57079633;
  // double pi = 3.14159265;

  // check f_real_eigens for 5 alike -- isotropy!!
  test = radius = theta  = 0.0;

  color[0] = color[1] = color[2] = 0.8; // pale grey
  color[3] = 1;
  
  
  //[xx xy xz yx yy yz zx zy zz ], symmetric, copy into matrix
  for ( i = 0; i < 3; i++) {
    for( j = 0; j < 3; j++ ) {
      theMatrix(i,j) = f_matrices[petal][i*3+j];
    }
  }
  
  
  lodeAngle( theMatrix,  radius,  theta);
  
  if (radius < eps ) { // no significant anisotropic component+

    color[0] = 0.9;
    color[1] = color[2] = 0.0; // dark red
    color[3] = 1.0;
    return;
  }
  theta = theta + 30; //degrees, rangle -30 to +30
  theta = theta/60; // -30 txc 0 shear +30 txe
  get_color_by_param ( color, theta );
 
}





/**
 * a way to classify deviatoric stress. r is like magnitude,
 * theta is between triaxial extension, shear, and triaxial 
 * compression
 * Source: Notes on Quantifying modes of a second order tensor
 * By Rebecca Brannon
 */ 
void GlyphViewer::lodeAngle( Matrix stress, double & radius, double & theta) 
{
  Matrix S(3,3);
  double mean;
  int i,j;

  double J_two,J_three; // deviatoric invariants

  double eps =  2.2204e-16; // matlab's epsilon

  J_two = J_three = 0.0;

  mean = 0.0;

  for ( i = 0; i < 3; i++ ) {
    mean = mean + stress(i,i);
  }
  mean = mean/3.0;


 // S = stress - 1/3(stress_kk) kronecker_ij
  for ( i = 0; i < 3; i++ ) {
    for (j = 0; j < 3; j++ ) {
      S(i,j) = stress(i,j);
      if (i == j)
	S(i,j) = S(i,j) - mean;
    }
  }

 
  J_two = 0.0;
  // magnitude of S_ij
  for ( i = 0; i < 3; i++ ) {
    for (j = 0; j < 3; j++ ) {
      J_two = J_two + S(i,j)*S(j,i);
     
    }
  }
  radius = J_two; // ||S||^2

  J_two = 0.5*J_two; // J_two = 1/2 ||S||^2


  // r = sqrt(|| S ||^2)
  if ( radius > eps*1000 ) // watch for underflow
    radius = sqrt(radius);
  else radius = 0.0;

  
  // J_3 = det(S)
  J_three = J_three + ( S(0,0)*S(1,1)*S(2,2) );
  J_three = J_three + ( S(0,1)*S(1,2)*S(2,0) );
  J_three = J_three + ( S(0,2)*S(1,0)*S(2,1) );
  J_three = J_three - ( S(0,2)*S(1,1)*S(2,0) );
  J_three = J_three - ( S(0,0)*S(1,2)*S(2,1) );
  J_three = J_three - ( S(0,1)*S(1,0)*S(2,2) );

  if ( radius == 0.0 )
    theta = -31; // impossible value -- theta must be +90 to -90 degrees 

  // theta = 1/3 asin[(J_3/2) * (3/J_2)^(3/2)]
  else { 
    theta = (3/J_two);

    theta = theta*theta*theta;

    theta = sqrt(theta); // theta ^(3/2)
 
    theta = theta*(J_three/2);

    theta = asin(theta);
   
    theta = theta/3.0;

    // change to degrees
    theta = radiansToDegrees(theta);
   
  }

  printf("Theta %f radius %f, mean stress %f\n", theta, radius, mean);
}

void GlyphViewer::setOrientation(){
  double * eigCopy;
  int dim[4] = {3,3,3,3};
  int i;
  double test;

  if (  theElement != NULL && 
	theElement->getClassTag() == ELE_TAG_EightNodeBrick 
	&& 	intgrPoint > -1 
	&& 	intgrPoint < 8	) {
    
    Tensor elastic( 4,dim, 0.0 );
    Tensor yieldNorm(2,dim,0.0);
    

    EightNodeBrick * el = (EightNodeBrick *) ( theElement );
      
    elastic =  el->ElasticStiffnessTensor( intgrPoint );
    
    yieldNorm = el-> dFods(intgrPoint );

    // finally, put eigentensor into a BJtensor
    // allocate and 
    try {
      eigCopy = new double[9];
    }
    catch( std::bad_alloc) {
      printf("ReynoldsGlyph::eigenMode2, unable to allocate data\n");
      exit(-1);
    }
    for (i = 0; i < 9; i++) {
      eigCopy[i] =  f_matrices[0][i];
    }
    stresstensor t_eig(eigCopy); // destructor will de-allocate
   
    test = testOrientation(yieldNorm, elastic, t_eig);
    if (test < 0) {
       for (i = 0; i < 9; i++) {
	 f_matrices[0][i] = -f_matrices[0][i];
       }
    }
    delete[] eigCopy; // clean up
  }
}
