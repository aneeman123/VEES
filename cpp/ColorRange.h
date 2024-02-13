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
** ColorRange.h: Header for
** class that handles color map parameters and scaling
** checks whether a color is within the selected range 
** to give a boolean decision whether to draw the object.
**  
**/


#ifndef COLOR_RANGE
#define COLOR_RANGE
#include <Element.h>
#include <Node.h>
class ColorRange {
 public:
  ColorRange();
  ColorRange( double baseValue, bool posDef, bool logScaling );

  // range cutoff values
  void setRange( double minRange, double maxRange, bool colInvert);
  bool inRange(double color);
  bool nodesInRange( double * theColors, Element * theElement );
  bool isInverted();


  // color mapping methods
  void setBasis( double baseValue, bool posDef, bool logScaling );
  void setDataName(int name);
  void setData( double * data, int length );
  void scaleDataSet();
  void scaleAndSetRange( double * scalar, int length, bool posDef, bool logScaling);
 
  //accessor functions
  int getDataName();
  bool isPositiveDefinite() const;
  bool isLogScaled() const;
  double getBaseValue() const; //determines max color value 
  double getColor( double value );
  double * getData(); // make this const later


 private:
  double min; // 0 to 1 value for draw/No draw
  double max; // 0 to 1 value for draw/No draw
  bool invert;  // invert selection for draw/no draw

  //colormap
  double baseVal;   // true range of values, either +/-baseVal or 0 to baseVal
  double logBaseVal;
  bool positiveDefinite; // actual values must be >= 0
  bool logScaled;          // color by log-scaled


  double * data;
  int dataLength;
  int colorSchemeName;

};




#endif
