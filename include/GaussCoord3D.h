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
** Appropriate for all OpenSees eight node bricks as of November, 2006.
**/

#ifndef _GAUSS_COORD_3D_
#define _GAUSS_COORD_3D_
// includes for the domain classes
#include <Domain.h>
#include <Node.h>
#include <Element.h>

#include <Matrix.h>
#include <Vector.h>
#include <basics.h>
#include <BJtensor.h>
#include <nDarray.h>
#include <stresst.h>
#include <straint.h>

#include <Information.h>

/* helper classes for draw */
#include <SingleDomEleIter.h>
#include <SingleDomNodIter.h>
#include <NodeIter.h>

#include <math.h>
#include <stdio.h>
#include <stdio.h>

#ifdef WIN32
#include <windows.h>
#endif

class GaussCoord3D{
 public:
  GaussCoord3D();

  tensor H_3D_8Node( double r1, double r2, double r3 );
  double get_Gauss_p_c( short order, short point_numb );

  int computeGaussPoint8Node( Node ** theNodes, Information &eleInfo );


 private:
    int r_integration_order; // Gauss-Legendre integration order in r direction
    int s_integration_order; // Gauss-Legendre integration order in s direction
    int t_integration_order; // Gauss-Legendre integration order in t direction


};


#endif
