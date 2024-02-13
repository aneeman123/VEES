
#include "ReynoldsGlyph.h"
#include "auxiliary.h"
#include "materials.h"
#include "Camera.h"
#include <Brick.h>




ReynoldsGlyph::ReynoldsGlyph() {
  int i;

  saveImage = 0;
  theElement = NULL;
  intgrPoint = -1;
  eleNum = -1;
  nonPositiveDefinite = 0;
  indeterminate = 0;
  rotate = new Matrix(6,6);
 

  try {
    sizeRange = new ColorRange();
    eigenModeColor = new double[3];
  }
  catch( std::bad_alloc) {
    printf("ReynoldsGlyph, Unable to allocate ColorRange\n");
    exit(0);
  }

  sizeRange->setDataName(1);
  
  drawMethod = REYNOLDS_GLYPH4R;
 

  // initialize eigentensors eigenvalues
  for ( i = 0; i < 6; i++ ) {
    f_real_eigens[i]  = 0.0;
  }


  drawStretch = 1;
  changedStiffness = 0; //no change
  smallestEigen = 0;

}

ReynoldsGlyph::~ReynoldsGlyph() {
  //delete ext;
  delete sizeRange;
  delete rotate;
}

void ReynoldsGlyph::setDrawMethod(int method) {
  if ( method == STRESS ) {
    drawMethod = method;
    setStress( intgrPoint );
    return;
  }
  if (method >= PETAL0R && method <= REYNOLDS_GLYPH4R) {
    drawMethod = method;
    setIntegrPoint( intgrPoint );
    return;
  }
 
  if ( method ==  EIGEN_STIFFNESSR ) {
      drawStretch = 0;
      setIntegrPoint( intgrPoint );
      return;
  }

  if ( method == EIGEN_STRETCHR ) {
    drawStretch = 1;
    setIntegrPoint( intgrPoint );
    return;
  }   
  
}




/**
 * precondition: color noormalized based on 4th order tensor
 */

void ReynoldsGlyph::setPetal( double eigSet[9] ) {
  int j, k;
  double max = -100000000;
  double min =  100000000;
  
  bool positiveDef = true;
	


  // create matrix
  Matrix theMatrix( eigSet, 3,3);
  
  //calculate morph constants
  for( j = 0; j < NUM_PTS_R; j++ ) { //for each exterior point
    p_ext[j] = einsteinSum2( theMatrix, ptLocR[j] );

   
    if ( p_ext[j] > max )
      max =  p_ext[j];
    if( p_ext[j] < min )
      min =  p_ext[j];
  }		

  //morph along ray
  // Since the tensors are symmetric, 
  // we use the absolute value and avoid 
  // turning the polygons inside out

  for ( j = 0; j < NUM_PTS_R; j++ ) {
    for ( k = 0; k < 3; k++ ) {
      p_extPts[j][k] = fabs(p_ext[j])*ptLocR[j][k]; 
      if (  p_extPts[j][k] < 0.0 )
	positiveDef = false;
    } //k
  } //j



  //calculate triangle normals & point normals
  setNormals(p_extPts, p_ptNorm);//, p_norm); 
}





/**
 * calculate normals for each point on the glyph */
void ReynoldsGlyph::setNormals(double extPts[NUM_PTS_R][3],
		       double ptNorm[NUM_PTS_R][3] ) {
  int i, j, count;
  double norm[NUM_TRI_R][3];
  for(i = 0; i < NUM_TRI_R; i++ ) {
		calcNormal( extPts[triPtsR[i][0]],extPts[triPtsR[i][1]], 
			extPts[triPtsR[i][2]], norm[i]);
	  }
  /* calculate normal for each point based on average
   * of surrounding triangles */
  for( i = 0; i < NUM_PTS_R; i++ ) {
   	for( j = 0; j < 6; j++ ) {
	  ptNorm[i][j] = 0.0;
	}
  }

  for( i = 0; i < NUM_PTS_R; i++ ) {
	count = 0;
	for( j = 0; j < 6; j++ ) {
	  if ( neighborsR[i][j] != -1 ) {
		ptNorm[i][0] = 	ptNorm[i][0] + norm[neighborsR[i][j]][0];
		ptNorm[i][1] = 	ptNorm[i][0] + norm[neighborsR[i][j]][1];
		ptNorm[i][2] =  ptNorm[i][0] + norm[neighborsR[i][j]][2];
		count++;
	  }
	} // 6 neighborsR
 
	ptNorm[i][0] /= (double)count;
	ptNorm[i][1] /= (double)count;
	ptNorm[i][2] /= (double)count;
  } //each point
}

/** 
 * draw the triangles if the stiffness tensor changed
 * @pre: color and mesh vs. solid already set */
void  ReynoldsGlyph::draw(double color){
  //  printf(" ReynoldsGlyph::draw() smallestEigen %d\n", smallestEigen);
  //fflush(stdout);

  if (drawMethod == STRESS ) {
    	drawPetal(color);

  }
  else if (changedStiffness) {
      /*      if (smallestEigen < 0 )
	drawE_ijkl();

	else */
    	drawPetal(color);
   
  }
 
} 

/**
 * color by eigenmode type : isotropy (red), shear (green),
 * axisymmetric (blue) or some combination of the modes (grey for all 3)
 */
void ReynoldsGlyph::drawEigenMode() {
   int i;
   float color[4];
   color[0] = color[1] = color[2] = 0.5;
   color[3] = 1.0;
  glPolygonMode(GL_FRONT_AND_BACK, GL_FILL);
  glEnable(GL_CULL_FACE);
  glCullFace(GL_BACK);

  if (indeterminate) {
    //DRAW THE BLACK BOX
    set_saturation_by_param(0.05, 1.0);
    glutSolidCube(0.9);
   
  }
  else if (isSpherical(f_real_eigens)) {
    set_color_by_param(1.0, 1.0); //red
    glutSolidSphere(0.6, 10,10);
  } 
  else {

  eigenMode2( color ); // get right color
  glMaterialfv( GL_FRONT_AND_BACK, GL_AMBIENT, color );
  glMaterialfv( GL_FRONT_AND_BACK, GL_DIFFUSE, color );

  glBegin( GL_TRIANGLES );

	for(i = 0; i < NUM_TRI_R; i++){

	  glNormal3d(p_ptNorm[triPtsR[i][0]][0], 
		     p_ptNorm[triPtsR[i][0]][1],
			   p_ptNorm[triPtsR[i][0]][2]);
	  
	  
	  
	  glVertex3d(p_extPts[triPtsR[i][0]][0],p_extPts[triPtsR[i][0]][1],
		     p_extPts[triPtsR[i][0]][2] );
	  
	  glNormal3d(p_ptNorm[triPtsR[i][1]][0], 
		     p_ptNorm[triPtsR[i][1]][1],
		     p_ptNorm[triPtsR[i][1]][2]);
	  
	  // set_color_minus_plus( p_ext[triPtsR[i][1]], 1.0 );
	  // set_color_minus_plus( f_real_eigens[smallestEigen], 1.0 );
	  glVertex3d(p_extPts[triPtsR[i][1]][0],
		     p_extPts[triPtsR[i][1]][1],
		     p_extPts[triPtsR[i][1]][2] );
	  
	  glNormal3d(p_ptNorm[triPtsR[i][2]][0], 
		     p_ptNorm[triPtsR[i][2]][1],
		     p_ptNorm[triPtsR[i][2]][2]);
	  
	  //	  set_color_minus_plus( p_ext[triPtsR[i][2]], 1.0 );
	  // set_color_minus_plus( f_real_eigens[smallestEigen], 1.0 );
	  glVertex3d(p_extPts[triPtsR[i][2]][0],
		     p_extPts[triPtsR[i][2]][1],
		     p_extPts[triPtsR[i][2]][2] );	
	}
   glEnd();
  }

}

/**
 * Draw a single eigentensor (with smalles eigenvalue)
 * as a Reynolds Glyph 
 * @pre: smallest eigenvalue has been set, and petal calculated
 */
void ReynoldsGlyph::drawPetal( double color ) {
  int i;

  glPolygonMode(GL_FRONT_AND_BACK, GL_FILL);
  glEnable(GL_CULL_FACE);
  glCullFace(GL_BACK);

  if (indeterminate) {
    //DRAW THE BLACK BOX
    set_saturation_by_param(0.05, 1.0);
    glutSolidCube(0.9);
   
  }
  else if (isSpherical(f_real_eigens)) {
    set_color_minus_plus( color, 1.0);
    glutSolidSphere(0.6, 10,10);
  } 
  else {
  glBegin( GL_TRIANGLES );

	for(i = 0; i < NUM_TRI_R; i++){

	  glNormal3d(p_ptNorm[triPtsR[i][0]][0], 
		     p_ptNorm[triPtsR[i][0]][1],
			   p_ptNorm[triPtsR[i][0]][2]);
	  
	  
	  //	  set_color_minus_plus( p_ext[triPtsR[i][0]], 1.0 );
	  //  set_color_minus_plus( f_real_eigens[smallestEigen], 1.0 );
	  set_color_minus_plus( color, 1.0);
	  glVertex3d(p_extPts[triPtsR[i][0]][0],p_extPts[triPtsR[i][0]][1],
		     p_extPts[triPtsR[i][0]][2] );
	  
	  glNormal3d(p_ptNorm[triPtsR[i][1]][0], 
		     p_ptNorm[triPtsR[i][1]][1],
		     p_ptNorm[triPtsR[i][1]][2]);
	  
	  // set_color_minus_plus( p_ext[triPtsR[i][1]], 1.0 );
	  // set_color_minus_plus( f_real_eigens[smallestEigen], 1.0 );
	  glVertex3d(p_extPts[triPtsR[i][1]][0],
		     p_extPts[triPtsR[i][1]][1],
		     p_extPts[triPtsR[i][1]][2] );
	  
	  glNormal3d(p_ptNorm[triPtsR[i][2]][0], 
		     p_ptNorm[triPtsR[i][2]][1],
		     p_ptNorm[triPtsR[i][2]][2]);
	  
	  //	  set_color_minus_plus( p_ext[triPtsR[i][2]], 1.0 );
	  // set_color_minus_plus( f_real_eigens[smallestEigen], 1.0 );
	  glVertex3d(p_extPts[triPtsR[i][2]][0],
		     p_extPts[triPtsR[i][2]][1],
		     p_extPts[triPtsR[i][2]][2] );	
	}
   glEnd();
  }
}



/**
 * @pre: orthographic projection to avoid squeezing the glyph */
/*void ReynoldsGlyph::drawE_ijkl() {
  int i;


  glPolygonMode(GL_FRONT_AND_BACK, GL_FILL);
  glEnable(GL_CULL_FACE);
  glCullFace(GL_BACK);

  glBegin( GL_TRIANGLES );
  for(i = 0; i < (NUM_TRI_R); i++ ) {

	glNormal3d(ptNorm[triPtsR[i][0]][0], ptNorm[triPtsR[i][0]][1],
			   ptNorm[triPtsR[i][0]][2]);
   	set_color_by_param_shiny( ext[triPtsR[i][0]], 0.1 );
	glVertex3d(extPts[triPtsR[i][0]][0],extPts[triPtsR[i][0]][1],
			   extPts[triPtsR[i][0]][2] );

	glNormal3d(ptNorm[triPtsR[i][1]][0], ptNorm[triPtsR[i][1]][1],
		 ptNorm[triPtsR[i][1]][2]);
	set_color_by_param_shiny( ext[triPtsR[i][1]], 0.1 );
	glVertex3d(extPts[triPtsR[i][1]][0],extPts[triPtsR[i][1]][1],
			   extPts[triPtsR[i][1]][2] );

	glNormal3d(ptNorm[triPtsR[i][2]][0], ptNorm[triPtsR[i][2]][1],
			   ptNorm[triPtsR[i][2]][2]);
	set_color_by_param_shiny( ext[triPtsR[i][2]], 0.1);
	glVertex3d(extPts[triPtsR[i][2]][0],extPts[triPtsR[i][2]][1],
			   extPts[triPtsR[i][2]][2] );
  }
  glEnd();
  }*/


double  ReynoldsGlyph::getCenterAndMaxDim( double center[3] ){
  center[0] = center[1] = center[2] = 0;

  return 3.0; 
  
}

int  ReynoldsGlyph::getNumElements() {
  return 1;
}

void  ReynoldsGlyph::setElement( Element * ele ) {
 
  theElement = NULL;
  theElement =  ele;
  eleNum = ele->getTag();
}

/**
 * Get current and elastic stiffness, matrices must be 6x6
 * also yield surface normal (a 3x3 symmetric tensor)
 * and plastic flow normal (3x3 symmetric)
 */
void ReynoldsGlyph::getStiffData( Matrix & m_elastic, Matrix & m_curr,
Matrix & yieldNorm, Matrix & flowNorm , int & isIndeterminate , int & sameStiff) {
  int dim[4] = {3,3,3,3};

  int i,j;

  isIndeterminate = indeterminate;


  if (  theElement != NULL && 
	theElement->getClassTag() == ELE_TAG_EightNodeBrick 
	&& 	intgrPoint > -1 
	&& 	intgrPoint < 8	) {
    
    Tensor init( 4,dim, 0.0 );
    Tensor curr( 4,dim, 0.0 );
    Tensor secOrder(2,dim,0.0);
   Tensor secOrderB(2,dim,0.0);

    EightNodeBrick * el = (EightNodeBrick *) ( theElement );
    curr = (el->getTangentTensor( intgrPoint ));
    
    init =  el->ElasticStiffnessTensor( intgrPoint );

    sameStiff = sameTensor( curr, init );
    
    Tensor2MatrixSysR4( init, m_elastic );
    Tensor2MatrixSysR4( curr, m_curr );

    secOrder = el-> dFods(intgrPoint );
    for (i = 0; i < 3; i++ ) {
      for ( j = 0; j < 3; j++) {
	yieldNorm(i,j) = secOrder.val(i+1,j+1);
      }
    }

   secOrderB = el-> dQods(intgrPoint );
    for (i = 0; i < 3; i++ ) {
      for ( j = 0; j < 3; j++) {
	flowNorm(i,j) = secOrderB.val(i+1,j+1);
      }
    }

    

  }

}

/**
 * get current stiffness matrix
 * @pre: m_curr is allocated 6x6 matrix
*/
void ReynoldsGlyph::getStiffness(Matrix & m_curr ) {
  int dim[4] = {3,3,3,3};
  if ( theElement->getClassTag() != ELE_TAG_EightNodeBrick ) {
    printf( "ReynoldsGlyph::getStiffness, element is not 8 node brick\n");
    return;
  }

  Tensor curr( 4,dim, 0.0 );
  EightNodeBrick * el = (EightNodeBrick *) ( theElement );
  curr = (el->getTangentTensor( intgrPoint ));
  Tensor2MatrixSysR4( curr, m_curr );
}


/**
 * Set the integration point, get its stiffness tensor, 
 * calculate eigentensors and their eigenvalues, update glyph
 */
void  ReynoldsGlyph::setIntegrPoint( int theNum ) {
  double * vr;
  int i,j;



  if (theNum != -1 ) { // actual hit
	intgrPoint  = theNum;
	
	for ( i = 0; i < 6; i++ ) {
	  f_real_eigens[i]  = 0.0;
	  smallestEigen = 0;
	}
	//reset rotation to identity
	for ( i = 0; i < 6; i++ ) {
	  for (j = 0; j < 6; j++ ) {
	    if (i == j)
	      (*rotate)(i,j) = 1.0;
	    else 
	      (*rotate)(i,j) = 0.0;
	  }
	}

	// do the big cast
       if ( theElement->getClassTag() == ELE_TAG_EightNodeBrick ) {

	 int dim[4] = {3,3,3,3};
	 //	 Tensor init( 4,dim, 0.0 );
	 //	 Tensor inverse( 4,dim, 0.0 );
	 Tensor curr( 4,dim, 0.0 );;
	 
	 //	 Matrix m_init(6,6);
	 //	 Matrix m_diff(6,6);
	 Matrix m_curr(6,6);
	 //	 Matrix m_inverse(6,6);
	
	 EightNodeBrick * el = (EightNodeBrick *) ( theElement );
	 curr = (el->getTangentTensor( intgrPoint ));

	 //	 init =  el->getInitialTangent();
	  
	 //	 if ( sameTensor(curr, init) ) { 
	 // changedStiffness = 0;
	 // return;
	 // }
	 // else 
	 changedStiffness = 1; //always draw

       
	 // change from tensor to matrix
	 //	 Tensor2MatrixSysR4( init, m_init );
	 Tensor2MatrixSysR4( curr, m_curr );
	 
	 // test for negative eigenvalue

	 // try to invert
	 // if (m_init.Invert( m_inverse ) != 0 ) {
	 // printf("ReynoldsGlyph::setIntgr, cannot invert initial stiff\n");
	 // return; // bail out
	 // }
	
	 // m_diff = m_curr*m_inverse; // C.inverse(C_0); 

	 // TEST FOR NONPOSITIVE DEFINITE 
	 testBifurcation(m_curr,nonPositiveDefinite, indeterminate );
	 if (!indeterminate) {
	 // DONE TESTING, CONTINUE
	   try {
	     vr = new double[36];
	   }
	   catch( std::bad_alloc) {
	     printf("ReynoldsGlyph::setIntgr, Unable to allocate eigenvectors\n");
	     exit(0);
	   }
	   Matrix stretch;
	   // Matrix rotate;
	 
	 /* if (highamPolar( m_diff, stretch,  rotate )  == 0) {
	   printf("ReynoldsGlyph::setIntgr, polar decomposition failed\n");
	   smallestEigen = 0;
	   return;
	   }*/

	   // tolerance smallest eigen within 4 orders of magnitude or 
	   // force Gramm Schmidt
	   if (vanGelderPolar( m_curr, *rotate, stretch, 1.0e12 ) != 0 ) {
	     printf("ReynoldsGlyph::setIntgr, polar decomposition failed\n");
	     smallestEigen = 0;
	     return;
	   }
	 
	   // try copying stretch
	   /*
	   if ( calcAsymmEigens( stretch, 6, vr, f_real_eigens ) != 0 ) {
	     printf("ReynoldsGlyph::setIntgr, calculating eigenvalues failed\n");
	     delete[] vr;
	     smallestEigen = 0;
	     return;
	     }*/
	   if( calcSymmEigens(  stretch, 6, vr, f_real_eigens) != 0 ) {
	     printf("ReynoldsGlyph::setIntgr, calculating eigenvalues failed\n");
	     delete[] vr;
	     smallestEigen = 0;
	     return;
	   }
	   
	   
		
	   sixVectorToSymMatrix( vr ); // turns them into 3x3 matrices
	  
	   setOrientation(); // in case we got (-1)*eigenvector
	   smallestEigen = 0; // they come out sorted
	   setPetal(f_matrices[smallestEigen]);
	  
	     

	   delete[] vr;
	  
	 } // !indeterminate, calculated polar 
	  
	 else {
	   // signal that rotation is non-meaningful
           // stretch is non meaningful
	 }
       } // is an 8 node brick
  
      
  } //hit
       
  else { 
    //  printf("miss\n:");
  }

}




/**
 * C_ij n_i n_j for appropriate number of indeces.
 * Constract the tensor to a scalar with each normalized
 * vector on the sphere. 
 */
/*void ReynoldsGlyph::morphGlyphPoints4(){
  int i, j;

 
  for( i = 0; i < NUM_PTS_R; i++ ) {
    ext[i] = einsteinSum4( ptLocR[i] ); 
  }

 
  sizeRange->setData( ext, NUM_PTS_R);
 

  if ( sizeRange->getBaseValue() == 1.0 ){ //default value
                               // scalar, length, posDef, logScaling
    sizeRange->scaleAndSetRange( ext, NUM_PTS_R, true, false ); 

    }
  else {
    sizeRange->scaleDataSet();
    }
  
  
 
  // normalizeMinusPlus( ext, NUM_PTS_R ); // worry about   //singles laterl;l 

  for ( i = 0; i < NUM_PTS_R; i++ ) {
	for ( j = 0; j < 3; j++ )
	  extPts[i][j] = ext[i]*ptLocR[i][j];
  } 
 
  } */ 

/**
 * Use for single Reynolds Stress/Strain glyph, not flower
 */
/*void ReynoldsGlyph::morphGlyphPoints2(  Matrix theMatrix ){

  int i, j;

  for( i = 0; i < NUM_PTS_R; i++ ) {
	ext[i] = einsteinSum2( theMatrix, ptLocR[i] ); 
	// ext[i] = HWYglyph( theMatrix, ptLocR[i] ); 
  }
  normalizeColorScale( ext, NUM_PTS_R, true ); // worry about   //singles later
  
  for ( i = 0; i < NUM_PTS_R; i++ ) {
	for ( j = 0; j < 3; j++ )
	  extPts[i][j] = ext[i]*ptLocR[i][j];
  } 

  }*/

/**
 * Tarantalo's unrolling of the 3x3x3x3 tensor into 9x9 
 * unrolls current tangent tensor
 */
/*
void ReynoldsGlyph::calcTarantalo() {
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
 

  }*/

/**
 * Stress or strain tensor. Hashah's method to show shear
 * sigma_S = sqrt{T^2 - sigma^2_N}  = sqrt{T_iTi - sigma^2_N} 
 * T_i = sigma_{ij}n_j
 * @returns the result of the formula or 0 if not a 2nd order tensor
 */
double ReynoldsGlyph::HWYglyph( Matrix theMatrix,double n[3] ) {
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
double ReynoldsGlyph::einsteinSum2( Matrix theMatrix, double n[3] ) {
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
/*double ReynoldsGlyph::einsteinSum4(  double n[3] ) {
  int  i, j, k,l;
  double sum;

  sum = 0.0;
  

  for( i = 0; i < 3; i++ ) {
	for ( j = 0; j < 3; j++ ) {
	  for ( k = 0; k < 3; k++ ) {
		for ( l = 0; l < 3; l++ ) {	
	   
		  sum = sum + tangent[i][j][k][l]*n[i]*n[j]*n[k]*n[l];
		} // l
	  } // k 
	} // j
  } // i

  return sum;
}  
*/


/**
 * copy from BJtensor into standard 4D array 
 * ( 3 x 3 x 3 x 3 ) 
 */
/*void ReynoldsGlyph::setTangentT( tensor theTensor )
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

}*/
 

/** map 6 x 6 Matrix => 3 x 3 x 3 x 3 tensor with minor symmetry */
/*void ReynoldsGlyph::setTangentM( Matrix theMatrix ) {
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
  
 


   }*/

//Matlab 9 x 9 matrix
/*
void ReynoldsGlyph::printTangent() {
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
  }*/

//Matlab 6 x 6 matrix
void ReynoldsGlyph::printMatrix(Matrix theMatrix) {
  int i,j;

  printf("M = [");
  for (i = 0; i < 6; i++) {
	for ( j = 0; j < 6; j++ ) {
	  printf("%16g ", theMatrix(j,i));//Matrix class indexing reversed
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
void ReynoldsGlyph::indexMap( int mIndex, int & i, int & j ) {

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
	printf ("Error ReynoldsGlyph::index map: called with %d, should be 0 to 5\n",
			mIndex);
			
  }
}

/**
 * take data that has both positive & negative values, 
 * normalize and center over 0-1. Range may be preset or calculated
 * on the fly
 */
void ReynoldsGlyph::normalizeColorScale( double * scalar, int length, bool setRange ) {
  int i;
  double glyphMax, glyphMin, range;
  bool positiveDefinite;
  positiveDefinite = false;
  range = 0;
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

void ReynoldsGlyph::normalizeMinusPlus( double * scalar, int length) {
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
int  ReynoldsGlyph::calcAsymmEigens(  Matrix  theMatrix, int leadingDim, double * vr, double * wr )
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
  double * wi; // imaginary part of eigenvalues
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
  wi = new double[lda];
  vl = new double[lda*lda];
  // vr = new double[lda*lda];
  work = new double[lwork];
  a = new double[lda*lda];

 /* LAPACK uses column first ordering, e.g.
   * |  5  7 |
   * | -2 -4 | == [ 5 -2 7 -4 ] */

  //  a = new double[lda*lda]; test as static array

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
  //  printf("The original matrix!!!\n");
  //opserr << theMatrix;

 
  // ok, let's run it
  dgeev_( &jobvl, &jobvr, &n, a, 
            &lda, wr, wi, vl, &ldvl, 
		       vr, &ldvr, work, &lwork, &info);


 
  //vr real eigenvectors
  //wr real eigenvalues

  // delete wr;
  delete wi;
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
/** SUBROUTINE DSYEVR(
    JOBZ, RANGE, UPLO, N, A, LDA, VL, VU, IL, IU, ABSTOL, M, W, Z, LDZ, ISUPPZ, WORK, LWORK, IWORK, LIWORK, INFO )

*/

/**
 * Zhou's code, returns 1 for failure*/
int ReynoldsGlyph::Matrix2TensorSysR4(const Matrix& M, tensor& T)
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
  
  T.val4(1,1,1,1) = M(0,0);
  T.val4(1,1,2,2) = M(0,1);
  T.val4(1,1,3,3) = M(0,2);      
  T.val4(1,1,1,2) = T.val4(1,1,2,1) = M(0,3) *sqrthalf;
  T.val4(1,1,2,3) = T.val4(1,1,3,2) = M(0,4) *sqrthalf;
  T.val4(1,1,1,3) = T.val4(1,1,3,1) = M(0,5) *sqrthalf;      

  T.val4(2,2,1,1) = M(1,0);
  T.val4(2,2,2,2) = M(1,1);
  T.val4(2,2,3,3) = M(1,2);      
  T.val4(2,2,1,2) = T.val4(2,2,2,1) = M(1,3) *sqrthalf;
  T.val4(2,2,2,3) = T.val4(2,2,3,2) = M(1,4) *sqrthalf;
  T.val4(2,2,1,3) = T.val4(2,2,3,1) = M(1,5) *sqrthalf; 

  T.val4(3,3,1,1) = M(2,0);
  T.val4(3,3,2,2) = M(2,1);
  T.val4(3,3,3,3) = M(2,2);      
  T.val4(3,3,1,2) = T.val4(3,3,2,1) = M(2,3) *sqrthalf;
  T.val4(3,3,2,3) = T.val4(3,3,3,2) = M(2,4) *sqrthalf;
  T.val4(3,3,1,3) = T.val4(3,3,3,1) = M(2,5) *sqrthalf;

  T.val4(1,2,1,1) = T.val4(2,1,1,1) = M(3,0) *sqrthalf;
  T.val4(1,2,2,2) = T.val4(2,1,2,2) = M(3,1) *sqrthalf;
  T.val4(1,2,3,3) = T.val4(2,1,3,3) = M(3,2) *sqrthalf;      
  T.val4(1,2,1,2) = T.val4(2,1,1,2) = T.val4(1,2,2,1) = T.val4(2,1,2,1) = M(3,3) *half;
  T.val4(1,2,2,3) = T.val4(2,1,2,3) = T.val4(1,2,3,2) = T.val4(2,1,3,2) = M(3,4) *half;
  T.val4(1,2,1,3) = T.val4(2,1,1,3) = T.val4(1,2,3,1) = T.val4(2,1,3,1) = M(3,5) *half;

  T.val4(2,3,1,1) = T.val4(3,2,1,1) = M(4,0) *sqrthalf;
  T.val4(2,3,2,2) = T.val4(3,2,2,2) = M(4,1) *sqrthalf;
  T.val4(2,3,3,3) = T.val4(3,2,3,3) = M(4,2) *sqrthalf;      
  T.val4(2,3,1,2) = T.val4(3,2,1,2) = T.val4(2,3,2,1) = T.val4(3,2,2,1) = M(4,3) *half;
  T.val4(2,3,2,3) = T.val4(3,2,2,3) = T.val4(2,3,3,2) = T.val4(3,2,3,2) = M(4,4) *half;
  T.val4(2,3,1,3) = T.val4(3,2,1,3) = T.val4(2,3,3,1) = T.val4(3,2,3,1) = M(4,5) *half;

  T.val4(1,3,1,1) = T.val4(3,1,1,1) = M(5,0) *sqrthalf;
  T.val4(1,3,2,2) = T.val4(3,1,2,2) = M(5,1) *sqrthalf;
  T.val4(1,3,3,3) = T.val4(3,1,3,3) = M(5,2) *sqrthalf;      
  T.val4(1,3,1,2) = T.val4(3,1,1,2) = T.val4(1,3,2,1) = T.val4(3,1,2,1) = M(5,3) *half;
  T.val4(1,3,2,3) = T.val4(3,1,2,3) = T.val4(1,3,3,2) = T.val4(3,1,3,2) = M(5,4) *half;
  T.val4(1,3,1,3) = T.val4(3,1,1,3) = T.val4(1,3,3,1) = T.val4(3,1,3,1) = M(5,5) *half;

  return 0;
} 



/** 
 * copied from Zhou's NewTemplate3Dep */
int ReynoldsGlyph::Tensor2MatrixSysR4(const tensor& T, Matrix& M)
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
  M(0,3) = T.cval(1,1,1,2) *sqrt2; // 1/sqrt2(T.cval(1,1,1,2) + T.cval(1,1,2,1))
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

//Zhou's code - invert a tensor
int ReynoldsGlyph::Stiffness2Compliance(const tensor& S, tensor& C)
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
 * are hardcoded. The assumption is that the vector has Mandel components
 * so we need to reverse that.
 */
void ReynoldsGlyph::sixVectorToSymMatrix( double eigTensors[36] ) {

  int i;
  int base;// base index for eigTensor array
 
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
 * @returns 1 for correct execution, 0 for failure 
 */
int ReynoldsGlyph::highamPolar( Matrix diffStiff, Matrix & stretch, Matrix & rotate ) {
  // iterative: D = RU, R is rotation, U is % stretch
  // R_0 = D, D = diffstiff
  int i;
  Matrix Rcurr(6,6);
  Matrix RcurrInvert(6,6);
  Matrix RcurrInverseTranspose(6,6);
  Matrix Rnext(6,6);
  
  Rcurr = diffStiff;
  
  for (i = 0; i < 7; i++) {

    if ( diffStiff.Invert(RcurrInvert ) != 0) { // already done by Frank.!
       printf("ReynoldsGlyph::highamPolar, unable to invert DIFFERENCE matrix\n");
       printf("Iteration %d\n", i);   
      return 0;
    }

    transposeMatrix( RcurrInvert,  RcurrInverseTranspose );

    // R_k+1 = 1/2(R_k + R_k^(-1T) ), where ^(-1T) means inverse transpose
    Rnext = (Rcurr + RcurrInverseTranspose);
    Rnext = 0.5*Rnext;

    // frobenius = differenceOfNorm( diffStiff, Rnext);
    // printf("highamPolar, iteration %d, curr %f next %f\n", i,Rcurr, Rnext);
  // check how we're doing by the Frobenius norm of the diff from the new iter and
  // original
    Rcurr = Rnext;

  } // Newton iterations

  rotate = Rcurr;

  // U = R^-1 D
  if ( Rcurr.Invert(RcurrInvert) != 0 )  {
    printf("ReynoldsGlyph::highamPolar, unable to invert rotate matrix to get stretch\n");
      return 0;
  } 
  
  stretch = RcurrInvert*diffStiff;
  return 1;
 
}



/*
 * @pre: both matrices have the same dimensions */
void ReynoldsGlyph::transposeMatrix( Matrix source, Matrix & transpose ) {
  int i,j,rows,cols;
  rows = source.noRows();
  cols = source.noCols();

  for (i = 0; i < rows; i++ )
    for (j = 0; j < cols; j++)
      transpose(i,j) = source(j,i);

}


/**
 * Find the difference of the Frobenius norm of 2 matrices */
double ReynoldsGlyph::differenceOfNorm(Matrix a, Matrix b) {
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
  // printf("norm %f\n", norm);
  norm = sqrt( norm );
  // printf("sqrt %f\n", norm);

  return norm;
}

/**
 * @pre: eigenvalues and eigenvectors have been calculated
 */
void ReynoldsGlyph::setSmallest() {
  int i;
   
  smallestEigen = 0;
  double smallVal =  f_real_eigens[0];
  for (i = 1; i < 6; i++) {
    if (f_real_eigens[i] < smallVal) {
      smallestEigen = i;
      smallVal = f_real_eigens[i];
    }
  }
  if ( f_real_eigens[smallestEigen] < -0.01 ) {
    
  printf("NEGATIVE smallest eigen %d value %f\n", 
	 smallestEigen, f_real_eigens[smallestEigen]);

  }
  
}



/**
 * contract the current stiffness aigainst the
 * eigentensor with the smallest eigenvalue
 * @pre! eigentensors have been calculated
 */
void ReynoldsGlyph::testSmallest( Tensor stiffness) {
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
  //printf("eigentensor\n");
  //eigen.print(); 

  Matrix m_stiff(6,6);
  Matrix m_compliance(6,6);
  Tensor compliance(4, pdim4, 0.0 );
  
  Tensor2MatrixSysR4( stiffness, m_stiff );
  m_stiff.Invert( m_compliance );

  Matrix2TensorSysR4( m_compliance, compliance );

  Tensor result = compliance("ijkl")*eigen("kl");

  // printf("Frobenius norm of contracted eigentensor: %f\n", result.Frobenius_norm());

 
  // try testing against standard orthonormal directions.
  for (i = 0; i < 6; i++) {
    Tensor  x( 2, pdim, orthoEigen[i] );
    result = compliance("ijkl")*x("kl");
    //  printf("Frobenius norm of contracted orthoEigen %d: %f\n", i, result.Frobenius_norm());
  }
}
  

/**
 * A cheap and not perfect way to compare doubles. Best I can do */ 
int ReynoldsGlyph::sameTensor( Tensor t1, Tensor t2) {
  int same = 1;
  int i,j,k,m;
 
  double bigger;
  double smaller;

  if ( t1.rank() != 4 || t2.rank() != 4) {
    printf("ReynoldsGlyph::compareTensors, mistake in call,must be rank 4 tensors\n");
    return 0;
  }

  for ( i = 1; i <=4; i++ ) {
    if (t1.dim(i) != 3 || t2.dim(1) != 3) {
      printf("ReynoldsGlyph::compareTensors, mistake in call,must be dimension 3 tensors\n");
      return 0;
    }
  }
 
  fflush(stdout);
  for( i =1; i <=3; i++) {
    for( j = 1; j <=3; j++ ) {
      for(k = 1; k <= 3; k++ ) {
	for (m = 1; m <=3; m++ ) {

	  if ( t1.val4(i,j,k,m) > 0 ) { // positive value: establish squeeze
	    bigger =  t1.val4(i,j,k,m)*1.00001;
	    smaller = t1.val4(i,j,k,m)*0.99999;	      
	  }
	  else if (t1.val4(i,j,k,m) == 0.0 ) {
	    bigger = 0.000001;
	    smaller = -0.000001;	     
	  }
	  else {
	    bigger =  t1.val4(i,j,k,m)*0.99999;
	    smaller =  t1.val4(i,j,k,m)*1.00001;
	  }
	  
	  if ( t2.val4(i,j,k,m) > bigger || t2.val4(i,j,k,m) < smaller ){
	    //   printf("%f squeezes [ %f to %f ], against %f and fails\n",
	    //   t1.val4(i,j,k,m),smaller,bigger,  t2.val4(i,j,k,m));
	    // fflush(stdout);
	    return 0;
	  
	  }

	}
      }
    }
  }
  return same; // they are the same within  0.01
}

int ReynoldsGlyph::getChangedStiffness() {
  return changedStiffness;
}




/*function [T] = orthoprojection(V,W)
* function [T] = orthoprojection(V,W) V is column vector, 
* W is a column vector 
* T is the result of projecting V against W ;
* you would subtract this whole thing from V to get the actual projection
*/
void ReynoldsGlyph::orthoprojection(Vector V, Vector W, Vector& T ) 
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
 * Symmetric part of a matrix is 1/2(T + T^T)
 */
void ReynoldsGlyph::symmetrizeMatrix(Matrix source, Matrix & symm ) {
  int i,j,rows,cols;

  rows = source.noRows();
  cols = source.noCols();

  Matrix trans(rows,cols);
  transposeMatrix(source,trans);

  for (i = 0; i < rows; i++ )
    for (j = 0; j < cols; j++)
      symm(i,j) = 0.5*( source(i,j) + trans(i,j) );

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
int  ReynoldsGlyph::calcSymmEigens(  Matrix  theMatrix, int leadingDim, 
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


  delete []work;
  return info;
}

extern "C" int dgetrf_(int *M, int *N, double *A, int *LDA, 
		       int * IPIV, int * INFO);

/**
 * @returns 0 for success */
int ReynoldsGlyph::determinant( Matrix theMatrix, double & det ) {
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
    printf("ReynoldsGlyph::determinant, unable to allocate work arrays\n");
    exit(0);
  }
 
 // copy the matrix into 'a'
  k = 0;
  // copy the Matrix into a
  //  printf("calculating determinant for matrix:\n");
  for( i = 0; i < lda; i++) {
    for(j = 0; j < lda; j++) {
      a[k] =  theMatrix(i,j); // matrix also stored in column order
      k++;
  
    }
  
  }
 

  dgetrf_(&m, &n, a, &lda, ipiv, &info);

  k = 0;

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
 
  for (i = 0; i < n; i++) {
    if (ipiv[i] != (i+1)) { // fortran starts counting at 1, thus i+1
      det = -det;
    
    }
  }

 
  delete []ipiv;
  delete [] a;
  return info;
}

/**
 * Calclates polar Decomposition of D. 
 * tolerance: how many orders of magnitude between smallest
 * and largest eigen before forcing Gramm-Schmist 
 * R: rotation matrix
 * U: stretch matrix
 * T_m : eigenvectors of U
 * sqrtJ: eigenvalues of u
 * @returns 0 for success
 */ 
int ReynoldsGlyph::vanGelderPolar( Matrix D, Matrix &R, Matrix &U, double toler)
{
  double * T; // eigenvectors
  double * J; // jacobean like eigenvalues
  int i,j,k, success, zeroEigs;
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
 DTD = DTD*D;

 try {
   T = new double[n*n];
   J = new double[n];
 }
 catch( std::bad_alloc) {
    printf("GlyphViewer::vanGelderPolar, Unable to allocate T &J arrays\n");
    exit(0);
 }

 //[T,J] = eig(DTD);
 // DTD may not be symmetric
 // calculate asymm eigens and sort
 success = calcSymmEigens(  DTD, n, T, J ); // get T,J


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

  j = determinant( T_m, det );

  if (j != 0) {
    printf("GlyphViewer::vanGelderPolar, 1st determinant calculation failed\n");   delete [] T;
   delete [] J;
   return -1;

  }
  
  //if det(T) < 0
  //  T(:,1) = -T(:,1);
  if (j == 0 && det < 0 ) {
    for (i = 0; i < n; i++) {
      T[i] = -T[i];
      T_m(i,0) = -T_m(i,0); 
    }
  }

  // precondition: sorted eigenvalues, which lapack gives us
  zeroEigs = 0;
  for (i = 0; i < n; i++ ) {
    if (J[i] < 0) {
      J[i]= 0;
    }
    if (J[i] <= toler * eps * J[n-1]){
      zeroEigs = zeroEigs + 1;
      J[i]= 0;
    }
  }

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

  
    if (zeroEigs == 0 ) { // non-singular U, R = DU^(-1)
      //U0 = (T)*SQRTJ*(T');
   
      transposeMatrix( T_m, trans );
      U = T_m*sqrtJ*trans;

      // Ensure U0 is symmetric
      //U0 = (U0 + U0') / 2;    
      transposeMatrix(U,trans);
      U = 0.5 * (U + trans);

      //R = DU^(-1)
      j = U.Invert(DTD);
      if ( j != 0 ) {
	printf("GlyphViewer::vanGelderPolar, error inverting U\n");
	 success -1;
      }
      else {
	R = D*DTD; //  D*U^(-1)
	success = 0;
      }

      // calcSymmEigens(  U, n, T, J ); 

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

     // DTF = CF
      G = (D*T_m)*F;


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

	//C = [C1 G(:, 2:n)];
	//R = C*(T');
	//U2 = inv(R) * D;   % more accurate than assuming inv(R) = R'.
	//U = (U2 + U2') / 2;	  
	transposeMatrix( T_m, trans );
	  
	R = G*trans;
	R.Invert(trans);
	
	U = trans*D;
	transposeMatrix( U, trans );
	U = ( U + trans)*0.5;


    } // end of case 1 zero eigenvalue

    

    delete [] T;
    delete [] J;
    return success;
}

/**
 * check the symmetric part of the tensor for zero or
 * negative eigenvalues
 */
void ReynoldsGlyph::testBifurcation( Matrix D, int &  nonPositiveDefinite,
				     int &   indeterminate ) {

 int i,j,rows,cols, count;
 double * a;
 double * w;
 
 double eps =  2.2204e-16*10000; //matlab's epsilon plus tolerance

  rows = D.noRows();
  cols = D.noCols();



  try{ 
    a = new double[rows*cols];
    w = new double[rows];
  }
  catch( std::bad_alloc) {
    printf("ReynoldsGlyph::testBifurcation, Unable to allocate eigens\n");
    exit(0);
  }
  

  j = calcAsymmEigens(  D, rows, a, w);
  // NOTE: eigenvalues are sorted

  if (j != 0 )
    return;


  nonPositiveDefinite = 0;
  indeterminate = 0; 
  count = 0;

  for (i = 0; i < 6; i++ ) {
   
    if (fabs(w[i]) < eps ) // we consider this close to zero
      count++;

    if (w[i] < -0.5) 
      indeterminate = 1;

  }
  if (count > 0 )
    nonPositiveDefinite = 1;
  if (count > 1)
    indeterminate = 1;


  delete [] a;
  delete [] w;
}

double ReynoldsGlyph::getLowestEigen() {
  if (indeterminate)
    return 0;
  else 
    return(f_real_eigens[smallestEigen]);
}

/** go through f_real_eigens, count similar eigenvalues within 0.5
 * if 5 the same, it's spherical
 * @returns true for spherical */
bool ReynoldsGlyph::isSpherical( double real_eigens[6] ) {
  int i,count1;

  count1 = 1;
  for (i = 1; i < 6; i++) {
    if (sameDouble(real_eigens[0],real_eigens[i])){

	count1++;
    }

  }
  if (count1 == 5)
    return true;
  
  count1 = 1;
  for (i = 0; i < 5; i++) {
    if (sameDouble(real_eigens[5],real_eigens[i])){

	count1++;
    }
  }
  if (count1 == 5)
    return true;

  return false;
}

/**
 * A cheap way to compare doubles within a small threshold
* @author Alisa Neeman aneeman@cse.ucsc.edu 11/29/2007
 */ 
int ReynoldsGlyph::sameDouble( double t1, double t2) {
  
  double bigger;
  double smaller;

  if (t1 == t2 ) //doubtful, but it could happen
    return 1;

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
 * red: isotropic part ( color[0] )
 * green: shear part ( color[1] )
 * blue: axisymmetric part ( color[2] ) 
 * R + G + B = 1.0
 */
void ReynoldsGlyph::eigenMode( float color[4] ) {
  double * a;
  double * w;
  Matrix theMatrix(3,3);
  int i,j, count;
  double test, iso, aniso, red, green, blue;
  int min, max;

  double eps =  2.2204e-16*10000; // matlab's epsilon*10000

  double w_abs[3];

  // check f_real_eigens for 5 alike -- isotropy!!

  if ( isSpherical(f_real_eigens) ) {
    color[0] = 1;
    color[1] = color[2] = 0;
    return;
  } 
  


  // otherwise do eigendecomposition on first eigentensor
  try{ 
    a = new double[9];
    w = new double[3];
    
  }
  catch( std::bad_alloc) {
    printf("ReynoldsGlyph::eigenMode, Unable to allocate eigens\n");
    color[0] = color[1] = color[2] = 0;
    exit(-1);
  }
  
  //[xx xy xz yx yy yz zx zy zz ], symmetric, copy into matrix
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

    // check for isotropic part
    iso = (w[0] + w[1] + w[2])/3;
    if ( fabs(iso) > eps ) {
      
      // remove isotropy
      for (i = 0; i < 3; i++ )
	w[i] = w[i] - iso;
      
    }
    // now, let's get the anisotropic part
    // find abs(min eigenval) and abs(max eigenval), 
    // divide min by max (range will be 0 .. 0.5)
    
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
   
    // 0 .. 0.5 mapped to 0 .. 1
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

    // 1.0 - iso = aniso
    color[0] = red;
    color[1] = green;
    color[2] = blue;
  
  } // able to calculate eigens
  
  
  delete[] a;
  delete[] w;
  
}



/**
 * Color by lode angle, or grey for spherical */
void ReynoldsGlyph::eigenMode2 ( float color[4] )
{
  Matrix theMatrix(3,3);
  int i,j, count;
  double test, radius,theta;
  int min, max;
  double w_abs[3]; // absolute value of eigenvalues from eigentensor

  int dim[4] = {3,3,3,3};
  double * eigCopy;

  double eps =  2.2204e-16*10000; // matlab's epsilon*10000
  double piOverTwo = 1.57079633;
  double pi = 3.14159265;

  // check f_real_eigens for 5 alike -- isotropy!!
  test = radius = theta  = 0;

  color[0] = color[1] = color[2] = 0.8; // pale grey
  color[3] = 1;

  if ( isSpherical(f_real_eigens) ) {
    return;
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
   
    test = testOrientation(yieldNorm, elastic, t_eig);

    if ( test < 0 ) {
      for ( i = 0; i < 3; i++) {
	for( j = 0; j < 3; j++ ) {
	  theMatrix(i,j) = theMatrix(i,j)*(-1);
	}
      }
    }

  
    lodeAngle( theMatrix,  radius,  theta);

    if (radius < 0.000001 ) {
      color[0] = 0.7;
      color[1] = color[2] = 0.0; // dark red
      color[3] = 1.0;
      return;
    }
    theta = theta + 30; //degrees, rangle -30 to +30
    theta = theta/60; // -30 txc 0 shear +30 txe
    get_color_by_param ( color, theta );
  }
}
 

/**
 * Test orientation of eigenvector associated with low eigenvalue
 * and reverse direction if appropriate
 */
void ReynoldsGlyph::setOrientation() {
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



/**
 * Test an eigentensor's orientation (both t and -1*t are valid eigentensors)
 ** -- match with Lode angle
 * YieldSurfNormal_ij Elastic_ijkl Eigentensor_kl > 0
 */
double ReynoldsGlyph::testOrientation(BJtensor  yieldNorm,  BJtensor elastic, BJtensor eig) {
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

void  ReynoldsGlyph::setStress( int theNum ) {
  int i;
  double norm;
  Response * theResponse = NULL;
  Information eleInfo(1.0);
  char * argStress[] = { "stress" };
  
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
	setPetal(stress);
    } // 8 nodes
}

/**
 * get the stress tensor for this integration point */
void ReynoldsGlyph::getStress( Matrix & stress ) {
  int i,j;
  double norm;
  Response * theResponse = NULL;
  Information eleInfo(1.0);
  char * argStress[] = { "stress" };
  
  for ( i = 0; i < 3; i++) {
      for ( j = 0; j < 3; j++ ) {
	stress(i,j) = 0.0;
      }
  }
  if (intgrPoint == -1 || theElement == NULL ||  
      theElement->getNumExternalNodes() != 8 ) {
    return;
  }

  
  theResponse = 
    theElement->setResponse((const char **)argStress, 1, eleInfo); 
  if (  theResponse ) {
    theResponse->getResponse();
    
    Information &theInfo = theResponse->getInformation();
    const Vector &eleData = theInfo.getData();
    
    if ( eleData.Size() == 48 ) {
      // xx yy zz xy xz yz
      
      stress(0,0)    = eleData(intgrPoint*6);
      stress(1,1) = eleData(intgrPoint*6+1);
      stress(2,2) = eleData(intgrPoint*6+2);
      stress(0,1) = stress(1,0) = eleData(intgrPoint*6+3);
      stress(0,2) = stress(2,0) = eleData(intgrPoint*6+4);
      stress(1,2) = stress(2,1) = eleData(intgrPoint*6+5);
      
    } // size 48
    
    else if ( eleData.Size() == 49 ) {
      
      stress(0,0) = eleData(intgrPoint*6+1);
      stress(1,1) = eleData(intgrPoint*6+2);
      stress(2,2) = eleData(intgrPoint*6+3);
      stress(0,1)  = stress(1,0) = eleData(intgrPoint*6+4);
      stress(0,2) = stress(2,0) = eleData(intgrPoint*6+5);
      stress(1,2) = stress(2,1) = eleData(intgrPoint*6+6);
      
      
    } // size 49
    delete theResponse;
  } // non-null response

}

/**
 * get the eigentensor for this integration point. 
 */
void ReynoldsGlyph::getEigenTensor( Matrix & eigen ) {
  int i,j;

  for ( i = 0; i < 3; i++) {
      for ( j = 0; j < 3; j++ ) {
	 eigen(i,j) = 0.0;
      }
  }
  if (intgrPoint == -1 || theElement == NULL ||  
      theElement->getNumExternalNodes() != 8 ) {
    return;
  }
  
  // if indeterminate, zero stress tensor is fine
  if (indeterminate)
    return;

  // if five equal stiffness eigenvalues, identity
  if (isSpherical(f_real_eigens)) {
    for ( i = 0; i < 3; i++) {
      eigen(i,i) = 1.0;
    }
    return;
  }

  // otherwise get lowest eigentensor
  for ( i = 0; i < 3; i++) {
      for ( j = 0; j < 3; j++ ) {
	eigen(i,j) = f_matrices[0][i*3 +j];
      }
  }

}

int ReynoldsGlyph::getRotationTensor(Tensor & t_rotate) {
  return Matrix2TensorSysR4(*rotate, t_rotate); 
 
}

void  ReynoldsGlyph::getRotationMatrix(Matrix & theMatrix) {
  int i,j;
  if (theMatrix.noRows() != rotate->noRows()
      || theMatrix.noCols() != rotate->noCols() ) {
    printf("ReynoldsGlyph::getRotationMatrix, mismatched matrix dimensions\n");
    return;
  }
  for (i = 0 ; i < rotate->noRows(); i++ ) {
    for ( j = 0; j < rotate->noCols(); j++ ) {
      theMatrix(i,j) = ((*rotate)(i,j));
    }
  }

}

/**
 * Stiffness is arbitrary stiffness matrix, 6x6
 * eigenmode is 3x3 symmetric matrix. On error will be zero matrix
 * lowEigen is the smallest eigenvalue. Note that this will not
 * test whether the direction of the eigenvector is correct!!
 * (i.e. (-1)*eigenvector is also an eigenvector)
 * Need testOrientation(yieldNorm, elastic, t_eig) for that.
 */
void ReynoldsGlyph::eigenStress( Matrix stiffness, Matrix & eigenmode, double & lowEigen ) {
  int i,j;
  double * vr;
  Matrix rotate(6,6);
  Matrix stretch(6,6);
  double real_eigens[6]; //real eigenvalues associated with each eigentensor
  int indeterm, nonposdef;
  indeterm = nonposdef = 0;
  // run vanGelderpolar, check for isotropy case
  // send isotropy or lowest eigentensor

  lowEigen = 0.0;

   for ( i = 0; i < 3; i++ ) {
      for ( j = 0; j < 3; j++ ) {
	eigenmode(i,j) = 0;
      }
    } 
	 // TEST FOR NONPOSITIVE DEFINITE 


  testBifurcation( stiffness, nonposdef, indeterm );

  if (indeterm) {
	 // DONE TESTING, CONTINUE
    lowEigen = -1.0;
    return;
  }

	   // force Gramm Schmidt
  if (vanGelderPolar( stiffness, rotate, stretch, 1.0e12 ) != 0 ) {
    printf("ReynoldsGlyph::eigenStress, polar decomposition failed\n");

    return;
  }

  
  try {
    vr = new double[36];
  }
  catch( std::bad_alloc) {
    printf("ReynoldsGlyph::eigenStress, Unable to allocate eigenvectors\n");
    return;
  }

  if( calcSymmEigens(  stretch, 6, vr, real_eigens) != 0 ) {
    printf("ReynoldsGlyph::eigenStress, calculating eigenvalues failed\n");
   
    delete[] vr;
    return;
  }
	   
  if ( isSpherical(real_eigens) ) {
    // identity
    for ( i = 0; i < 3; i++ ) {
      for ( j = 0; j < 3; j++ ) {
	if (i == j )
	  eigenmode(i,j) = 1;
	else
	  eigenmode(i,j) = 0;
      }
    }
    delete[] vr;
    lowEigen = real_eigens[0];
   
   return;
  }
  //  printf("interpolated stiffness\n");
  // opserr << stiffness;

  //  printf("interpolated stretch\n");
  //opserr << stretch;

  // printf("interpolated stretch eigens\n");
  // for (i = 0; i < 6; i++)
  // printf("%f ", real_eigens[i]);
  // otherwise take first eigentensor (lowest eigenvalue)
 
  eigenmode(0,0) = vr[0];
  eigenmode(1,1) = vr[1];
  eigenmode(2,2) = vr[2];
  eigenmode(0,1) =  eigenmode(1,0) = vr[3];
  eigenmode(1,2) =  eigenmode(2,1) = vr[4];
  eigenmode(0,2) =  eigenmode(2,0) = vr[5];
  
  lowEigen = real_eigens[0];
  

}


/**
 * a way to classify deviatoric stress. r is like magnitude,
 * theta is between triaxial extension, shear, and triaxial 
 * compression. theta is between -30 and +30 if radius is non-zero
 * it will be -31 if radius is zero
 * Source: Notes on Quantifying modes of a second order tensor
 * By Rebecca Brannon
 */ 
void ReynoldsGlyph::lodeAngle( Matrix stress, double & radius, double & theta) 

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

    // change radians to degrees
    theta = radiansToDegrees(theta);
  }

}
