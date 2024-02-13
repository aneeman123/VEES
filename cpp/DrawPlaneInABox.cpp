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
** DrawPlaneInABox.cpp is a class to draw stress data at a finite element's
** integration points. The stress is drawn as a plane composed of two 
** eigenvectors, corresponding to the eigenvalues with the greatest 
** magnitudes.
**/
#include "DrawPlaneInABox.h"
#include <OPS_Stream.h>
#include <StandardStream.h>
#include <Recorder.h>
#include <ElementRecorder.h>
#include <Information.h>
#include <Response.h>
#include <f2c.h>

DrawPlaneInABox::DrawPlaneInABox( Element * el ) {
  int i;
  //if the hashtable is already there and its the right size, 
  //leave it. Othere waise allocate & fill it.
  if (MC_hash == NULL) {
	MC_hash = new HashTable(257);
	for( i = 0; i < 256; i++ ) {
	  // copy point into allocated, insert into hash 
	  MC_hash->insert( reverse_lookup[i][0], reverse_lookup[i][1]);
	}
  }
  if ( canDraw( el ) ) {
    boxes = NULL;
    triangles = NULL;
    totalTriangles = 0;
    
    theElement = el;
    
    allocateTriangles();
    num_triangles = 0;
    for ( i = 0; i < 8; i++ )
      e_num_tri[i] = 0;
    
    boxEightNodes(); //allocate and create boxes around gauss pts
    calculateEigens();
    for ( i = 0; i < 8; i++ ) { 
      diagonalize(i);
      
      sortEigens(i, true ); //sort by unsigned magnitude
      planeInABox(i); //test
      
    }
  }

}


DrawPlaneInABox::~DrawPlaneInABox() {
  int i;
  //free boxes.
  for ( i = 0; i < numPts; i++ )
    delete boxes[i];
  delete [] boxes;
  boxes = NULL;
}


/* draw all the triangles in the plane in
 * the appropriate color and opacity 
 * pre!: caller push and pop viewpoint matrix 
 *       caller does glBegin(GL_TRIANGLES), glEnd()
 * (more efficient of caller does it once.)
 */
 
void DrawPlaneInABox::drawPlane( int i, double * theColors, 
		       ColorRange range, int colorflag  ) {
  int j, pick;


 glGetIntegerv(GL_RENDER_MODE,&pick); 
  if( pick == GL_SELECT )
	glPushName( i ); //make gauss point pickable

 
  else if ( colorflag == BY_ELEMENT ) {
	set_color_by_param( theColors[theElement->getTag()], 1.0 );
  }
  
  else if ( colorflag == BY_GAUSS_PT ) {
	// punt!!
	set_color_by_param( theColors[theElement->getTag()], 1.0 );
  }

  glBegin(GL_TRIANGLES);
   for( j = 0; j < e_num_tri[i]; j++ ) {
     glNormal3dv(triangles[e_start_tri[i]+j].vertex[0].n);
     glVertex3dv(triangles[e_start_tri[i]+j].vertex[0].xyz);
	 
     /* tweak normal to get a little shading */
     glNormal3d(triangles[e_start_tri[i]+j].vertex[0].n[X] + 0.1,
				triangles[e_start_tri[i]+j].vertex[0].n[Y] + 0.1,
				triangles[e_start_tri[i]+j].vertex[0].n[Z] + 0.1);
     glVertex3dv(triangles[e_start_tri[i]+j].vertex[1].xyz);
	 
     glNormal3dv(triangles[e_start_tri[i]+j].vertex[0].n);
     glVertex3dv(triangles[e_start_tri[i]+j].vertex[2].xyz);
   }
   glEnd();

   glGetIntegerv(GL_RENDER_MODE,&pick); 
   if( pick == GL_SELECT )
	 glPopName(); 
}

/**
 * center is average of node coordinates in element */
void DrawPlaneInABox::setCenter(  ) {

  Node ** theNodes; // nodes for a single element
  Node * theNode;
  Vector v;
  int i,j,numNodes, numCrds;


  if (theElement == NULL ) {
	printf("Error, null element in call to  DrawPlaneInABox::setCenter\n");
	center[0] = center[1] = center[2] = 0;
	return;
  }
  
  theNodes = theElement->getNodePtrs();
  numNodes = theElement->getNumExternalNodes();

  //make the first node both min and max.
  //go through the rest of the nodes and get bigger or smaller x, y, z values
  // if no z values, just use zero for z. 
  if ( theNodes == NULL ) {
	printf("DrawPlaneInABox::setCenter,Error fetching element nodes\n");
	return;
  }
  theNode = theNodes[0];
  v = theNode->getCrds(); 
  numCrds = v.Size();

  center[0] = center[1] = center[2] = 0;

  for (j = 0; j < numNodes; j++ ){	  

	v = theNodes[j]->getCrds(); 

	for( i = 0; i < numCrds; i++ ) {

	  center[i] = center[i] + v[i];

	} 

  } // while more nodes

  //find average of all nodes (unweighted centroid :)
  for ( i = 0; i < numCrds; i++ ) {
	center[i] = center[i]/numNodes;
  }

}


/* given a marching cubes index and intersecting edge points,
 * and point index (for normals)
 * create triangles and place them in global triangle buffer 
 * Assumption: No ambiguous cases (no gradients needed) 
 */
void DrawPlaneInABox::marchingCube(int i, int mc_index,  double mc_edges[12][3] ) {
  int j,k;
  j = k = 0;

  if(triangles == NULL  ) {
    printf("DrawPlaneInABox::marchingCube,triangles haven't been allocated\n");
	return;
  }

  // loop through all active edges used to form triangles 
   while (tri_cases[mc_index][j] > -1 ) {
	 if (! (num_triangles < totalTriangles ) ) {
	   printf("DrawPlaneInABox::marchingCube,insufficient triangles \n");
	   return;
	 }
	 // grab three edges 
	 for( k = 0; k < 3; k++ ) {
	   
	   // intersecting point on edge 
	   triangles[num_triangles].vertex[k].xyz[X] = 
		 mc_edges[tri_cases[mc_index][j]][X];
	   triangles[num_triangles].vertex[k].xyz[Y] = 
		 mc_edges[tri_cases[mc_index][j]][Y];
	   triangles[num_triangles].vertex[k].xyz[Z] = 
		 mc_edges[tri_cases[mc_index][j]][Z];
	   
	   // write normals; should be same for all  
	   triangles[num_triangles].vertex[k].n[X] = eVectors[i][X];
	   triangles[num_triangles].vertex[k].n[Y] = eVectors[i][Y];
	   triangles[num_triangles].vertex[k].n[Z] = eVectors[i][Z];
	   
	   triangles[num_triangles].p_num = i;
	   
	   j++;
	 }// one triangle 
	 num_triangles++;
   }
}


/* Algorithm by Paul Bourke
 * http://astronomy.swin.edu.au/~pbourke/geometry/planeline/
 * for a given box around a point, find
 * where eigenvector plane intersects the box edges
 * Create index for edge intersection.
 * convert to marching cubes index.
 * Do marching cubes on intersection points.
 * @pre: MC_hash must be allocated and values inserted 
 * create triangles, add to list to be drawn.
 * NOTE: POINT MUST HAVE NON-NEGATIVE X,Y,Z!!!
 * @returns 1 for success, 0 for failure
 * @param i is the ith gauss point
 */
int DrawPlaneInABox::planeInABox( int i ) {
  double  A,B,C,D,t, i_t, t_denom;
  double diff; /* shift amount if t = 0 error */
  double delta[3]; /* length of box edge */
  double n0[3]; /* normal of point */
  double touch_plane[3]; /* where line touches plane */
  double box[8][3]; /* points at each corner of box */
  int j,count,index,k;
  double mc_edges[12][3]; /* passed to marching cubes */
  int * mc_index;

  double plus = 0.0; //ensure point & box are positive definite.
  
  if (boxes[i][X_LESS] < plus )
	plus = boxes[i][X_LESS];
  if (boxes[i][Y_LESS] < plus )
	plus = boxes[i][Y_LESS];
  if (boxes[i][Z_LESS] < plus )
	plus = boxes[i][Z_LESS];
  
  plus = fabs(plus);
  
  // put the correct corner points in the box 
  for( j = 0; j < 8; j++ ) {
    box[j][X] = boxes[i][box_points[j][X]] + plus;/* either X_LESS or X_MORE */
    box[j][Y] = boxes[i][box_points[j][Y]] + plus;/* either Y_LESS or Y_MORE */
    box[j][Z] = boxes[i][box_points[j][Z]] + plus;/* either Z_LESS or Z_MORE */
  } 
 
 diff = 0;
  
  for( index = 0; index <12; index++ )
    for( j = 0; j < 3; j++ )
      mc_edges[index][j] = 0;

  // Implicit form of plane normal dot (Po->P), P unknown.
  // Apx + Bpy + Cpz + D = 0 
  // A = n_x, D = -n_x*x0 -n_y*y0 -nz*z0 
  n0[X] = eVectors[i][X]; //ASSUMPTION: DESIRED NORMAL SORTED TO FIRST COLUMN
  n0[Y] = eVectors[i][Y];
  n0[Z] = eVectors[i][Z];  

  A = n0[X]; B = n0[Y]; C = n0[Z]; 

  D = (-1*n0[X]*(gPoints[i][X] + plus)) - 
  (n0[Y]*(gPoints[i][Y]+ plus)) - (n0[Z]*(gPoints[i][Z]+plus));
 
 if (D == 0) { // t cannot be calculated if D is 0 (line on plane) 
   
   diff = (box[1][X] - box[0][X])/1000; 
   
   // shift box, recalculate D with shifted points[i][X] 
   for( j = 0; j < 8; j++ ) {
	 box[j][X] = box[j][X] + diff;
   }
   
   D = (-1*n0[X]*(gPoints[i][X] +plus) + diff)
	 - (n0[Y]*(gPoints[i][Y] +plus)) - (n0[Z]*(gPoints[i][Z] + plus));
   
 }
 
 
 // for each edge, check for intersection and interpolate 
 index = 0;
 for( j = 0; j < 12; j++ ) {
   delta[X] = box[cube_edges[j][1]][X] - box[cube_edges[j][0]][X];
   delta[Y] = box[cube_edges[j][1]][Y] - box[cube_edges[j][0]][Y];
   delta[Z] = box[cube_edges[j][1]][Z] - box[cube_edges[j][0]][Z];
   
   // calculate numerator 
   // t = i_t = -1*(A*box[cube_edges[j][0]][X] + B*box[cube_edges[j][0]][Y] +
   //	C*box[cube_edges[j][0]][Z] + D);
   t = i_t = (A*box[cube_edges[j][0]][X] + B*box[cube_edges[j][0]][Y] + 
			  C*box[cube_edges[j][0]][Z] + D);
   
   t = (-1)*t;
   // denominator 
   t_denom = A*delta[X] + B*delta[Y] + C*delta[Z];
   //printf("%f ", t_denom); 
   t = t/t_denom;
   
   // intersection is x0 + dX*t 
   touch_plane[edge_dim[j]] = box[cube_edges[j][0]][edge_dim[j]] + 
	 delta[edge_dim[j]]*t;
   
   // Does plane intersect edge? 
   if (touch_plane[edge_dim[j]] > box[cube_edges[j][1]][edge_dim[j]]
	   || touch_plane[edge_dim[j]] < box[cube_edges[j][0]][edge_dim[j]])
	 ;
   else {  // edge and plane intersect 
	 if (t < 0.0 || t > 1.0 )  
	   ; 
	 // printf("error!!!\n");
	 else {
	   
	   // copy point- shift X back if necesssary 
	   // only one dimension is not the same for point 0 and 1 
	   mc_edges[j][X] =  box[cube_edges[j][1]][X];
	   mc_edges[j][Y] =  box[cube_edges[j][1]][Y];
	   mc_edges[j][Z] =  box[cube_edges[j][1]][Z];
	   
	   mc_edges[j][edge_dim[j]] =   touch_plane[edge_dim[j]];
	   
	   // check whether this point is already in the set 
	   
	   index = index + CASE_INDEX[j];
	 }
	 if( diff != 0 ) {
	   mc_edges[j][X] =  mc_edges[j][X] - diff; 
	 } 
	 // remove shift 
	 mc_edges[j][X] -= plus;
	 mc_edges[j][Y] -= plus;
	 mc_edges[j][Z] -= plus;
   }
 } // for each edge 
 //  printf("%d MC index: %d\n",i, index); 
 index = MC_hash->get_value( index ); // criss-cross MC reference 
 //allocate triangles as need; we will need to count theses as we make them
 // average case: 2 triangles per intgr. point. we will do 2.5 to be safe


 if(index == -1  ) {
   printf("Error trying to make marching cubes index! broken at point %d\n",
		  i);
   e_start_tri[i] = 0;
   e_num_tri[i] = 0;
   //     error_list[i] = 1;
   return 0;
 } 
 else {
   // set e_tri to triangle[num_triangles] 
   e_start_tri[i] = num_triangles; //current triangle
   j = num_triangles;
   marchingCube(i, index,  mc_edges ); // create triangles
   
   // check difference in num_triangles to find out how 
   // many triangles were made
   e_num_tri[i] = num_triangles - j;
   
 }
 printf("marched!!!\n");
 return 1; 
}
/* shift the box along the given dimension (X, Y, or Z)
 * before trying to form a marching cubes index.
 * returns: the marching cubes index or -1 if unsuccessful.
 */
int plane_in_a_box2(int i, int stretch_dim,   float box[8][3], float A, float B, float C, float D) {
   int index, j; /*
  float delta[3]; // edge length 
  float t, t_denom;
  float touch_plane[3]; // where line touches plane 
  float mc_edges[12][3]; // passed to marching cubes
  float diff;

  
  diff = (box[1][stretch_dim] - box[0][stretch_dim])/1000; // kludge
    for( j = 0; j < 8; j++ ) {
       box[j][stretch_dim] = box[j][stretch_dim] + diff;
    }

  

   for( j = 0; j < 12; j++ ) {
     delta[X] = box[cube_edges[j][1]][X] - box[cube_edges[j][0]][X];
     delta[Y] = box[cube_edges[j][1]][Y] - box[cube_edges[j][0]][Y];
     delta[Z] = box[cube_edges[j][1]][Z] - box[cube_edges[j][0]][Z]; 

     // calculate numerator 
 
     t = (A*box[cube_edges[j][0]][X] + B*box[cube_edges[j][0]][Y] + 
		C*box[cube_edges[j][0]][Z] + D);
    
     t = (-1)*t;
     // denominator 
     t_denom = A*delta[X] + B*delta[Y] + C*delta[Z];

     t = t/t_denom;
     
     // intersection is x0 + dX*t 
     touch_plane[edge_dim[j]] = box[cube_edges[j][0]][edge_dim[j]] + 
       delta[edge_dim[j]]*t;
     
     // Does plane intersect edge? 
     if (touch_plane[edge_dim[j]] > box[cube_edges[j][1]][edge_dim[j]]
	 || touch_plane[edge_dim[j]] < box[cube_edges[j][0]][edge_dim[j]])
       ;
     else {  // edge and plane intersect 
       if (t < 0.0 || t > 1.0 )  
	 ; 

       else {
       
		 // copy point- shift X back if necesssary 
		 // only one dimension is not the same for point 0 and 1 
	 mc_edges[j][X] =  box[cube_edges[j][1]][X];
	 mc_edges[j][Y] =  box[cube_edges[j][1]][Y];
	 mc_edges[j][Z] =  box[cube_edges[j][1]][Z];
	 
	 mc_edges[j][edge_dim[j]] =   touch_plane[edge_dim[j]];
	 
	 // check whether this point is already in the set 
	 
	 index = index + CASE_INDEX[j];
       }
	 if( diff != 0 ) {
	   mc_edges[j][X] =  mc_edges[j][X] - diff;
	   
	 } 
       }
   } // for each edge 
   //  printf("%d MC index: %d\n",i, index);  
   index = get_value( MC_hash, index ); //criss-cross MC reference

   // reset box in case this didn't fix it 
    for( j = 0; j < 8; j++ ) {
       box[j][stretch_dim] = box[j][stretch_dim] - diff;
    }
*/
   return index;
}


/**
 * Create boxes around each integration point
 * from closest node in brick and unweighted centroid 
 * of element
 */
void DrawPlaneInABox::boxEightNodes() {
  Response * theResponse = NULL;
  Node ** theNodes;
  double node[3];
  int i, j, closest;
  double curr,dist;
  Vector vNodes;
 
  theNodes = theElement->getNodePtrs();
  
  Vector v;//good
  Information info(v);//good
  
  int success = theElement->getResponse (  5, info );  // 5 gets gauss for UCDavis ele's  
  if (success == 0 ) { 
	const Vector & gaussData = info.getData();
	if ( &gaussData == NULL ) {
	  printf(" DrawPlaneInABox::boxEightNodes: No gauss data\n");
	  return;
	}
	numPts = (int)gaussData(0);
	
	setCenter();
	
	//allocate boxes
	if (boxes == NULL ) 
	  boxes = new double*[numPts]; //first entry is number of gauss pts

	if (boxes == NULL ) {
	  printf(" DrawPlaneInABox::boxEightNodes:");
	  printf("Out of memory allocating plane in a box boxes, element %d\n",
			 theElement->getTag());
	  return;
	}

	for( i = 0; i < numPts; i++ ) {
	  boxes[i] = new double[6];
	  if (boxes[i] == NULL){
	  printf(" DrawPlaneInABox::boxEightNodes:");
		printf("Out of memory allocating plane in a box boxes, element %d\n",
			   theElement->getTag());
		return;
	  }
	} 

	// set boundaries
	for ( i = 0; i < numPts; i++ ) {
	  
	  boxes[i][X_MORE] = boxes[i][X_LESS] = gaussData[i*3+1];
	  boxes[i][Y_MORE] = boxes[i][Y_LESS] = gaussData[i*3+2];
	  boxes[i][Z_MORE] = boxes[i][Z_LESS] = gaussData[i*3+3];
	  
	  gPoints[i][X] = gaussData[i*3+1];
	  gPoints[i][Y] = gaussData[i*3+2];
	  gPoints[i][Z] = gaussData[i*3+3];
	  
	  //find closest Node. bounding box comes from 
	  //nearest Node and centroid
	  curr = 100000000;
	  closest = -1;
	  for (j = 0; j < numPts; j++ ) {
		vNodes = theNodes[j]->getCrds();  
		
		dist = distance( gPoints[i], vNodes );
		if (dist < curr) {
		  curr = dist;
		  closest = j;
		}
	  } // find closest node

	  vNodes = theNodes[closest]->getCrds();  
	  
	  // now, the box. go through closest x,y,z Node
	  
	  if ( vNodes[X] > boxes[i][X_MORE] ) 
		boxes[i][X_MORE] = vNodes[X];
	
	  else if ( vNodes[X] < boxes[i][X_LESS] ) 
		boxes[i][X_LESS] = vNodes[X];
	  
	  if ( vNodes[Y] > boxes[i][Y_MORE] ) 
		boxes[i][Y_MORE] = vNodes[Y];

	  
	  else if ( vNodes[Y] < boxes[i][Y_LESS] ) 
		boxes[i][Y_LESS] = vNodes[Y];

	  if ( vNodes[Z] > boxes[i][Z_MORE] ) 
		boxes[i][Z_MORE] = vNodes[Z];

	  else if ( vNodes[Z] < boxes[i][Z_LESS] )
		boxes[i][Z_LESS] = vNodes[Z];

	  // and go through centroid x,y,z
	  if ( center[X] > boxes[i][X_MORE] )
		boxes[i][X_MORE] = center[X];
	
	  else if ( center[X] < boxes[i][X_LESS] )
		boxes[i][X_LESS] = center[X];
	  
	  if ( center[Y] > boxes[i][Y_MORE] )
		boxes[i][Y_MORE] = center[Y];
	
	  else if ( center[Y] < boxes[i][Y_LESS] ) 
		boxes[i][Y_LESS] = center[Y];

	  
	  if ( center[Z] > boxes[i][Z_MORE] ) 
		boxes[i][Z_MORE] = center[Z];
	
	  else if ( center[Z] < boxes[i][Z_LESS] ) 
		boxes[i][Z_LESS] = center[Z];
	 
	} // for each integration pt.
  } //if rec'd Response

  else printf("No response, element %d request %d\n", theElement->getTag(), 5);
  //check how recorder frees data
}

//lapack/blas function
#ifdef _WIN32

extern "C" int  DSYEVD(char *jobz, char *uplo, int *n, 
			       double *a, int *lda, double *w, double * work, int *lwork,
					  int * iwork, int * liwork, int *info);

#else

extern "C"  int  dsyevd_(char *jobz, char *uplo, int *n, 
			       double *a, int *lda, double *w, double * work, int *lwork,
					  int * iwork, int * liwork, int *info);

#endif

/**
 * fetch stress tensor, calculate eigenvectors and eigenvalues
 */
void DrawPlaneInABox::calculateEigens( ) {
  Response * theResponse = NULL;
  int i,j,success, numPts;

  
  if ( theElement == NULL ) {
	printf("DrawPlaneInABox::calculateEigens( ):");
	printf("Need to set element to get gauss planes\n");
	return;
  }
  Vector v;//good
  Information info(v);//good
  success = theElement->getResponse ( 4, info );  // 4 gets stress for UCDavis ele's 
  if (success != 0 ) {
	printf("DrawPlaneInABox::calculateEigens( ):");
	printf(" No stress response 4 from element %d\n", theElement->getTag());
	return;
  }
  const Vector & InfoS = info.getData();
  if ( &InfoS == NULL ) {
	printf("DrawPlaneInABox::calculateEigens( ):");
	printf( "No gauss stress data\n");
	return;
  }
  numPts = (int)InfoS(0);
  for (i = 0; i < 8; i++ ) {
	for (j = 0; j < 9; j++ ) {
	  eVectors[i][j] = 0.0;
	}
  }

  //single array of size i*j (rows * cols)
  //where to put sigma_ij: i+(j*cols) 
  //get eigenvectors from element
  for (i = 0; i < numPts; i ++ ) {
	// upper triangle of matrix which lives in array
	eVectors[i][0] = InfoS(i*6+1); //sigma_xx
	eVectors[i][4] = InfoS(i*6+2); //sigma_yy
	eVectors[i][8] = InfoS(i*6+3); //sigma_zz
	eVectors[i][1] = eVectors[i][3] = InfoS(i*6+4); //sigma_xy
	eVectors[i][2] = eVectors[i][6] = InfoS(i*6+5); // sigma_xz
	eVectors[i][5] = eVectors[i][7] = InfoS(i*6+6); // sigma_yz	
		

  }

}


/**
 *  calculate eigenvectors and eigenvalues for
 *  ith gauss pt's stress or strain 
 * @pre eVector matrix has been filled with stress or
 * strain values from ith gauss point.
*/ 
int DrawPlaneInABox::diagonalize(int i) {
  
  char *jobz = "V"; // Compute eigenvalues and eigenvectors
  char *uplo = "U";  // Upper triagle of matrix is stored
  int n = MATRIX_ORDER; //order of the matrix

  int lda = MATRIX_ORDER; //length of row of a

  // work is scratch space
  double * work = new 
	double[1 + (6*MATRIX_ORDER) + (2*MATRIX_ORDER*MATRIX_ORDER)];

  int lwork = 1 + (6*MATRIX_ORDER) + (2*MATRIX_ORDER*MATRIX_ORDER);
  int * iwork = new int[3 + (5*MATRIX_ORDER)];
  int liwork = 3 + (5*MATRIX_ORDER);

  int info = 0; // returns success or failure

#ifdef _WIN32
  unsigned int sizeC = 1;
  DSYEVD(jobz, uplo, &n, eVectors[i], &lda, eValues[i], work, &lwork,
					  iwork,  &liwork, &info);
#else
 dsyevd_(jobz, uplo, &n, eVectors[i], &lda, eValues[i], work, &lwork,
					  iwork,  &liwork, &info);
#endif

 delete []work;
 delete []iwork;
 // do something with eVectors and eValues;
 return info;
}

/**
 * Sort gauss point's eigenvalues and eigenvectors
 * by ascending value. Note: lapack writes
 * vectors to the array column by column
 * @pre: eigenvectors & values have been set 
 * @param i  gauss point ID
 * @param unsignedVal  sort by signed or unsigned value
*/
void DrawPlaneInABox::sortEigens(int i, bool unsignedVal ) {
  int j,k;
  int order[3]; //MIN, MAX, MED
  double vals[9];
  int place;
  double temp;

 
  for (j = 0; j < 3; j++ ) {
   	if ( unsignedVal )
	  vals[j] = fabs(eValues[i][j]);
	else 
	  vals[j] = eValues[i][j];
  }

 //determine sorted order of values
  //simple scheme: if other eigen is greater, subtract 1
  // if lesser, add one. sums to +2, 0, or -2
  for (j = 0; j < 3; j++) {
	place = 0;
	for (k = 0; k < 3; k++) {
	  if (j != k ) {
		if (vals[j] > vals[k]) 
		  place +=1;

		else if (vals[j] < vals[k]) 
		  place -=1;
	  }
	}
   
	order[j] = place/2 + 1; // shift to 0,1,2
  }
  
  // swap eigenvalues
  for( j = 0; j < 3; j++ ) {
	vals[j] = eValues[i][order[j]];
  }
  for( j = 0; j < 3; j++ ) {
	eValues[i][j] = vals[j];
  }

  //swap eigenvectors.
  for( j = 0; j < 3; j++ ) {
	for (k = 0; k < 3; k++) {
	  vals[j*3 + k] = eVectors[i][order[j]*3 + k];
	}
  }
  for( j = 0; j < 9; j++ ) {
	eVectors[i][j] = vals[j];
  }
}

/**
 * we expect 2.5 triangles per box. One box per gauss point
 * set counter to 1st triangle
 */
void DrawPlaneInABox::allocateTriangles() {
  int num = (int)(2.5*theElement->getNumExternalNodes());

  if ( num > totalTriangles ) { //allocate more triangles

	if ( triangles != NULL )
	  free( triangles );
	
	triangles = (Triangle *)malloc(	sizeof(Triangle) * num );
	
	if ( triangles == NULL ) {
	  printf("DrawPlaneInABox::allocateTriangles, unable to allocate.\n"); 
	  exit (-1);
	}

	totalTriangles = num;
  }
}

/**
 * @pre! all planes have been set, triangles allocated and calculated, etc.
 */
void DrawPlaneInABox::drawElement(Element * el, double * colors, 
							  ColorRange range, int colorStyleFlag )
{
  int i, num;
  if (theElement ->getTag() != el->getTag() ) {
	printf(" DrawPlaneInABox::drawElement, wrong tag number\n");
	return;
  }
  num = theElement->getNumExternalNodes();
  for (i = 0; i < num; i++ )
	drawPlane( i, colors,  range, colorStyleFlag );
}

/** return if this method is capable of drawing 
 * this element type */
bool  DrawPlaneInABox::canDraw( Element * el )
{ 
  int tag = el->getClassTag();
  return (tag ==   ELE_TAG_EightNodeBrick ||
	  tag == ELE_TAG_EightNodeBrick_u_p_U);
}

/** distance = sqrt( dx^2 + dy^2 + dz^2) */
double  DrawPlaneInABox::distance( double p1[], Vector p2 ) {
  double delta[3];
  int i;
  double sum = 0;
  
  for(i = 0; i < p2.Size(); i++ ) {

	delta[i] = p1[i] - p2[i];
	delta[i] = delta[i]*delta[i];
	sum += delta[i];
  }


  return (sqrt( sum )); 
}

/**
 * Find out what DrawMethod this is */
int DrawPlaneInABox::getClassTag(){
  return DRAW_PLANE_IN_A_BOX;
}

