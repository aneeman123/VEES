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
** DESIGNER:          Boris Jeremic, Zhaohui Yang and Xiaoyan Wu
** PROGRAMMER:        Boris Jeremic, Zhaohui Yang  and Xiaoyan Wu
** Alisa Neeman (aneeman@cse.ucsc.edu) 
**
** GaussCoord3D provides the set of gauss integration point locations 
** for 8 and 20 node isoparametric brick.
** Appropriate for all OpenSees eight node brick types as of November, 2006.
**/

#include "GaussCoord3D.h"


#define FixedOrder 2

Vector Gauss8(FixedOrder*FixedOrder*FixedOrder*3+1); //Gauss point coordinates

GaussCoord3D:: GaussCoord3D() {

  r_integration_order = s_integration_order = t_integration_order = FixedOrder;
}

/**
 * H is the displacement interpolation tensor,
 * and it depends on the desired point in local coordinates
 * r1, r2, r3 are in -1 to +1 coord system
 */
tensor GaussCoord3D::H_3D_8Node(double r1, double r2, double r3) {
  int dimension[] = {24,3};
  
  tensor H(2, dimension, 0.0);

  // influence of the node number 8
  H.val(22,1) = (1.0 + r1)*(1.0 - r2)*(1.0 - r3)*0.125;
  H.val(23,2)= H.val(22,1);
  H.val(24,3)= H.val(22,1);
  
  // influence of the node number 7
  H.val(19,1) = (1.0 - r1)*(1.0 - r2)*(1.0 - r3)*0.125;
  H.val(20,2) = H.val(19,1);
  H.val(21,3) = H.val(19,1);
  
  // influence of the node number 6
  H.val(16,1) = (1.0 - r1)*(1.0 + r2)*(1.0 - r3)*0.125;
  H.val(17,2) = H.val(16,1); 
  H.val(18,3) = H.val(16,1);
  
  // influence of the node number 5
  H.val(13,1) = (1.0 + r1)*(1.0 + r2)*(1.0 - r3)*0.125;
  H.val(14,2) = H.val(13,1);
  H.val(15,3) = H.val(13,1);
  
  // influence of the node number 4
  H.val(10,1) = (1.0 + r1)*(1.0 - r2)*(1.0 + r3)*0.125;
  H.val(11,2) = H.val(10,1);
  H.val(12,3) = H.val(10,1);
  
  // influence of the node number 3            
  H.val(7,1) = (1.0 - r1)*(1.0 - r2)*(1.0 + r3)*0.125;
  H.val(8,2) = H.val(7,1); 
  H.val(9,3) = H.val(7,1);
  
  // influence of the node number 2
  H.val(4,1) = (1.0 - r1)*(1.0 + r2)*(1.0 + r3)*0.125;
  H.val(5,2) = H.val(4,1); 
  H.val(6,3) = H.val(4,1); 
  
  // influence of the node number 1
  H.val(1,1) = (1.0 + r1)*(1.0 + r2)*(1.0 + r3)*0.125;
  H.val(2,2) = H.val(1,1);
  H.val(3,3) = H.val(1,1);
  
  return H;

}

double GaussCoord3D::get_Gauss_p_c(short order, short point_numb){
//  Abscissae coefficient of the Gaussian quadrature formula
// starting from 1 not from 0
    static double Gauss_coordinates[7][7];

    Gauss_coordinates[1][1] = 0.0 ;

    Gauss_coordinates[2][1] = -0.577350269189626;
    Gauss_coordinates[2][2] = -Gauss_coordinates[2][1];

    Gauss_coordinates[3][1] = -0.774596669241483;
    Gauss_coordinates[3][2] = 0.0;
    Gauss_coordinates[3][3] = -Gauss_coordinates[3][1];

    Gauss_coordinates[4][1] = -0.861136311594053;
    Gauss_coordinates[4][2] = -0.339981043584856;
    Gauss_coordinates[4][3] = -Gauss_coordinates[4][2];
    Gauss_coordinates[4][4] = -Gauss_coordinates[4][1];

    Gauss_coordinates[5][1] = -0.906179845938664;
    Gauss_coordinates[5][2] = -0.538469310105683;
    Gauss_coordinates[5][3] = 0.0;
    Gauss_coordinates[5][4] = -Gauss_coordinates[5][2];
    Gauss_coordinates[5][5] = -Gauss_coordinates[5][1];

    Gauss_coordinates[6][1] = -0.932469514203152;
    Gauss_coordinates[6][2] = -0.661209386466265;
    Gauss_coordinates[6][3] = -0.238619186083197;
    Gauss_coordinates[6][4] = -Gauss_coordinates[6][3];
    Gauss_coordinates[6][5] = -Gauss_coordinates[6][2];
    Gauss_coordinates[6][6] = -Gauss_coordinates[6][1];

    return Gauss_coordinates[order][point_numb];
}

/**
 * This works the same as getResponse from Elements
 * @returns 0 for success
 */
int 
GaussCoord3D::computeGaussPoint8Node( Node ** theNodes, Information &eleInfo ) 
{   
    int count;
    count = FixedOrder*FixedOrder*FixedOrder;
    //Vector Gsc(count*3+1); //+1: number of Gauss point in element
    Gauss8(0) = count;

    double r  = 0.0;
    double s  = 0.0;
    double t  = 0.0;

    short where = 0;

    // special case for 8 nodes only
    static const int dim[] = {3, 8}; // static-> see ARM pp289-290
    static const int dimM[] = {3,  FixedOrder*FixedOrder*FixedOrder}; 
    tensor NodalCoord(2, dim, 0.0);
    tensor matpointCoord(2, dimM, 0.0);
    int h_dim[] = {24,3};   // Xiaoyan changed from {60,3} to {24,3} for 8 nodes
    tensor H(2, h_dim, 0.0);
    
    //Zhaohui using node pointers, which come from the Domain
    const Vector &nd1Crds = theNodes[0]->getCrds();
    const Vector &nd2Crds = theNodes[1]->getCrds();
    const Vector &nd3Crds = theNodes[2]->getCrds();
    const Vector &nd4Crds = theNodes[3]->getCrds();
    const Vector &nd5Crds = theNodes[4]->getCrds();
    const Vector &nd6Crds = theNodes[5]->getCrds();
    const Vector &nd7Crds = theNodes[6]->getCrds();
    const Vector &nd8Crds = theNodes[7]->getCrds();
    
    NodalCoord.val(1,1) = nd1Crds(0); 
    NodalCoord.val(2,1) = nd1Crds(1); 
    NodalCoord.val(3,1) = nd1Crds(2);

    NodalCoord.val(1,2) = nd2Crds(0); 
    NodalCoord.val(2,2) = nd2Crds(1); 
    NodalCoord.val(3,2) = nd2Crds(2);

    NodalCoord.val(1,3) = nd3Crds(0); 
    NodalCoord.val(2,3) = nd3Crds(1); 
    NodalCoord.val(3,3) = nd3Crds(2);

    NodalCoord.val(1,4) = nd4Crds(0); 
    NodalCoord.val(2,4) = nd4Crds(1); 
    NodalCoord.val(3,4) = nd4Crds(2);

    NodalCoord.val(1,5) = nd5Crds(0); 
    NodalCoord.val(2,5) = nd5Crds(1); 
    NodalCoord.val(3,5) = nd5Crds(2);
 
    NodalCoord.val(1,6) = nd6Crds(0); 
    NodalCoord.val(2,6) = nd6Crds(1); 
    NodalCoord.val(3,6) = nd6Crds(2);

    NodalCoord.val(1,7) = nd7Crds(0); 
    NodalCoord.val(2,7) = nd7Crds(1); 
    NodalCoord.val(3,7) = nd7Crds(2);

    NodalCoord.val(1,8) = nd8Crds(0); 
    NodalCoord.val(2,8) = nd8Crds(1); 
    NodalCoord.val(3,8) = nd8Crds(2);

    for( short GP_c_r = 1 ; GP_c_r <= r_integration_order ; GP_c_r++ )
      {
        r = get_Gauss_p_c( r_integration_order, GP_c_r );
        for( short GP_c_s = 1 ; GP_c_s <= s_integration_order ; GP_c_s++ )
          {
            s = get_Gauss_p_c( s_integration_order, GP_c_s );
            for( short GP_c_t = 1 ; GP_c_t <= t_integration_order ; GP_c_t++ )
              {
                t = get_Gauss_p_c( t_integration_order, GP_c_t );
                // this short routine is supposed to calculate position of
                // Gauss point from 3D array of short's
                where =
		         ((GP_c_r-1)*s_integration_order+GP_c_s-1)*t_integration_order+GP_c_t-1;
                // derivatives of local coordinates with respect to local coordinates
				
		H = H_3D_8Node(r,s,t);
		
		for (int encount=1 ; encount <= 8 ; encount++ )  { // move the coord by the node
		                                                    // weighted by the interpolation function
		  matpointCoord.val(1,where+1) +=
		    NodalCoord.val(1,encount) * H.val(encount*3-2,1); //X
		  
		  matpointCoord.val(2,where+1) +=
		    NodalCoord.val(2,encount) * H.val(encount*3-1,2); //Y
		  
		  matpointCoord.val(3,where+1) +=
		    NodalCoord.val(3,encount) * H.val(encount*3-0,3); //Z
		  
		}
		
		Gauss8(where*3+1) = matpointCoord.val(1,where+1);
		Gauss8(where*3+2) = matpointCoord.val(2,where+1);
		Gauss8(where*3+3) = matpointCoord.val(3,where+1);
				
				//matpoint[where].reportTensor("");
              }
          }
      }

	return  eleInfo.setVector( Gauss8 );
}

