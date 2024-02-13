
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
** ColorRange.cpp: 
** class that handles color map parameters and scaling
** checks whether a color is within the selected range 
** to give a boolean decision whether to draw the object.
**/

#include "ColorRange.h"
#include <stdio.h>

ColorRange::ColorRange(){
  min = 0.0;
  max = 1.0;
  baseVal = 1.0;   // true range of values, either +/-baseVal or 0 to baseVal 
  logBaseVal = 0.0;
  positiveDefinite = true; // actual values >= 0
  logScaled = false; 
  invert = false;
  data = NULL;
  colorSchemeName = -1;
}

/**
 * Set attributes of color map: the range of values, whether 
 * they are positive definite, and whether to calculate 
 * colors on linear or log scale
 * @param baseValue must be positive
 */
ColorRange::ColorRange( double baseValue, bool posDef, bool logScaling )
{
  min = 0.0;
  max = 1.0;
  invert = false;

  baseVal = baseValue;
  logBaseVal = logf(baseValue);

  positiveDefinite = posDef;
  logScaled = logScaling;

}



/**
 * Set cutoff range parameters. Caller must check that min <= max
 */
void 
ColorRange::setRange( double minRange, double maxRange, bool colInvert)
{
  if( minRange > maxRange )
	printf("Error, minimum grater than maximum, range unchanged.\n");
  else {
	min = minRange;
	max = maxRange;
	invert = colInvert;
  }

}

/**
 * Check if a given value is within the range
 */
bool ColorRange::inRange(double color) {
  if ( invert == false &&
	   color >= min && color <= max )
	return true;

  if( invert == true &&
	  ( color < min || color > max ))
	return true;

  return false;

}

/**
 * Deterministically decide whether to draw an element or not.
 * The logic: if there exists a node that is out of range, don't draw.
 * The inverse: if all nodes are out of (inverted) range, don't draw.
 */
bool ColorRange::nodesInRange( double * theColors, Element * theElement ) {
  int i;
  Node ** theNodes; // nodes for a single element
  int   numNodes = theElement->getNumExternalNodes(); 
  int count = 0;
  //check for within range
 
  theNodes = theElement->getNodePtrs();
  for ( i = 0; i < numNodes; i++ ) {
	if ( invert ) {
	  if ( ! inRange(theColors[theNodes[i]->getTag()]) )
		count ++;
	}
	
	else if ( !inRange(theColors[theNodes[i]->getTag()]) )
	  return false; // one in range; bail out.
  }
  if (count == numNodes ) 
	return false; //all in inverted range, don't draw
  
  return true;
}


/**
 * Is cutoff range inverted?
 */
bool ColorRange::isInverted() {
  return invert;
}


/**
 * Set attributes of color map: the range of values, whether 
 * they are positive definite, and whether to calculate 
 * colors on linear or log scale
 * @param baseValue must be positive
 */
void ColorRange::setBasis( double baseValue, bool posDef, bool logScaling )
{
  baseVal = baseValue;
  logBaseVal = logf(baseValue);

  positiveDefinite = posDef;
  logScaled = logScaling;

}

void ColorRange::setDataName(int name) {
  colorSchemeName = name;
}

void ColorRange::setData( double * theData, int len ) {
  if ( data != NULL && data != theData ) {
    printf("ColorRange %d, deleting old data set\n", colorSchemeName);
    delete[] data;

  }
  data = theData;
  dataLength = len;
}

int ColorRange::getDataName() {
  return colorSchemeName;
}

bool ColorRange::isPositiveDefinite() const{
  return positiveDefinite;
}

bool ColorRange::isLogScaled() const{
  return logScaled;
}

double ColorRange::getBaseValue() const {
  return baseVal;
}

/**
 * Returns a value between 0 and 1 for use as
 * color mapping. If the value is outside the 
 * given range, it will get the clamped value (-1, 0 or +1).
 * @param value assumed to be raw value, not log of value.
 * @pre logScaled, positiveDefinate, baseVal and logBaseVal are set
 */
double ColorRange::getColor( double value ) {
  double color;
 
  if (!logScaled) {
    if (positiveDefinite)
      color = value/baseVal;
    else 
      color = (value + baseVal)/(2*baseVal);
  }
  else { // log scaled
   
    if( value > -1 && value < 1) { // proper fraction makes negative logs
       color = 0.0;
    }

    else { // improper fraction makes correct log
      
      color = fabs(value); // no log of neg. numbers
      
      if ( positiveDefinite) {
	color = logf(color)/logBaseVal;
      }
      
      else { // -1 to + 1 shifted to 0 to 1
	color =  logf(color);

	if ( value < 0.0 )   
	  color = -color; 

	  // add twice
	  color = (logBaseVal + color)/(2*logBaseVal);
      }
    }// improper fraction
  }// is log scaled
  return color;
}

/**
 * Use existing range to scale data values. Clamp into 0 to 1
 * for positive definite data.
 */
void ColorRange::scaleDataSet() {
  int i;
  bool clamped = false;
  bool negVals = false;

  double min = 5000000;
  double max = -500000;
  if (data == NULL) {
    printf("Error in ColorRange::scaleDataSet(), data set is NULL\n");
    return;
  }

 
 for( i = 0; i < dataLength; i++ ) {    
    if (data[i] < min )
      min = data[i];
    if ( data[i] > max)
      max = data[i];
  }
   printf("  ColorRange::scaleDataSet(), data min %f & data max %f\n",
	  min, max ); 
  /*if (logScaled) {

    printf("ColorRange::scaleDataSet: Error, log scale not yet implemented\n");
    return;
    } */

 if ( data != NULL ) {
  
    for( i = 0; i < dataLength; i++ ) {      
      data[i] = getColor( data[i] ); // handles log scale as well
      if (data[i] > 1.0 ) {
	data[i] = 1.0;
	clamped = true;
      }
      if( data[i] < 0.0 ) {
	data[i] = 0.0;
	negVals = true;
      } 
    }// for each value
    if (clamped)
      printf("Warning: please expand data color range; had to clamp\n");
    if( negVals )
      printf("Warning: some scaled data values  are not positive\n");
  } //positiveDefinite data set

  min = 5000000;
  max = -500000;
  for( i = 0; i < dataLength; i++ ) {    
    if (data[i] < min )
      min = data[i];
    if ( data[i] > max)
      max = data[i];
  }
  printf("ending min %f max %f basVal %f logBaseVal %f\n", 
	 min,max, baseVal, logBaseVal);

}


/**
 * ignoring log scale for the moment */
void ColorRange::scaleAndSetRange( double * scalar, int length, bool posDef, bool logScaling) 
{

  int i;
  double max, min, range;

  bool pos;

  max = -100000000;
  min =  100000000;

  pos = true;
 
  // establish range 
  for ( i = 0; i < length; i++ ) {
	if ( scalar[i] > max )
	  max = scalar[i];
	if( scalar[i] < min )
	  min = scalar[i];
	if (scalar[i] < 0.0 )
	  pos = false;
  } 

  //  printf("before scaling, min %f max %f baseval %f\n", min,max,baseVal);

  if (posDef && !pos ) {

    printf("ColorRange::scaleAndSetRange:\n");
    printf("Data is not positive definite as specified\n");
    printf("No scaling done\n");
    return;
  
  }

  positiveDefinite = posDef;

  /* // special case: all data is same, nothing to scale
  if ( (max - min) < 0.00001  &&
       (max - min) > -0.00001 ) {
    printf("M");
    baseVal  = fabs(min);
    if ( min == 0.0 ) {
      baseVal = 1.0; //default value
      logBaseVal = 0.0;
      return;
    }
    if (baseVal > 1.0) // log of proper fraction is negative
     logBaseVal = logf( baseVal);
    else 
      logBaseVal = 0.0;

    return;
    } */

  
  if(fabs(min) >= fabs(max))
	range = fabs( min );
  else range = fabs(max);
 
  if (range < 1.0)
    baseVal = 1.0;
  else 
    baseVal = range;

  if (baseVal > 1.0) // log of proper fraction is negative
    logBaseVal = logf( range );
  else 
    logBaseVal = 1.0;

  for( i = 0; i < length; i++) {

    scalar[i] = getColor( scalar[i] );
  }

 
  data = scalar;
  dataLength = length;

  max = -100000000;
  min =  100000000;
 for ( i = 0; i < length; i++ ) {
	if ( scalar[i] > max )
	  max = scalar[i];
	if( scalar[i] < min )
	  min = scalar[i];
	if (scalar[i] < 0.0 )
	  pos = false;
  } 
 //  printf("after scaling, min %f max %f scale %f\n", min,max,baseVal);

}

double * ColorRange::getData() {
  return data;
}
