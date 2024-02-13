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
** ScalarVolume is a class to turn gauss point data from irregular
** finite elements into a regular grid.
**/
#include <stdlib.h>
#include <OPS_Globals.h>
#include <StandardStream.h>
#include <ArrayOfTaggedObjects.h>

// includes for the domain classes
#include <Domain.h>
#include <Node.h>


/* helper classes for draw */
#include <SingleDomEleIter.h>
#include <SingleDomNodIter.h>
#include <NodeIter.h>

#include <math.h>
#include <stdio.h>

class ScalarVolume {
 public:
  ScalarVolume(Domain * d);
  ~ScalarVolume();
  void saveVolumeToFile();
  void setValue( double value, int x, int y, int z );
 private:
  double *** grid;
  Vector bounds;
  double min; // min distance between nodes
  double spaceX, spaceY, spaceZ; // grid spacing

};
