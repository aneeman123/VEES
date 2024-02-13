#include "StressGlyph.h"
#include <Information.h>
#include <Response.h>
StressGlyph::StressGlyph() {
  for (i =0; i < 9; i++)
    f_stress[i] = 0;
}

void  StressGlyph::setIntegrPoint( int theNum ) {

  Response * theResponse = NULL;
  Information eleInfo(1.0);
  char * argStress[] = { "stress" };

  if (theNum != -1 ) { // actual hit
    intgrPoint  = theNum;
    
    for ( i = 0; i < 9; i++ ) {
      f_stress[i]  = 0.0;
      
    }
    
    if (theElement != NULL ) {  
      if ( theElement->getNumExternalNodes() == 8 ) {
	//try get stress
	theResponse = theElement->setResponse((const char **)argStress, 1, eleInfo); 
	if (  theResponse ) {
	  theResponse->getResponse();
	  
	  Information &theInfo = theResponse->getInformation();
	  const Vector &eleData = theInfo.getData();
	  
	  if ( eleData.Size() == 48 ) {
	    // xx yy zz xy xz yz
	    
	    f_stress[0] = eleData(intgrPoint*6);
	    f_stress[4] = eleData(intgrPoint*6+1);
	    f_stress[8] = eleData(intgrPoint*6+2);
	    f_stress[1] = f_stress[3] = eleData(intgrPoint*6+3);
	    f_stress[2] = f_stress[6] = eleData(intgrPoint*6+4);
	    f_stress[5] = f_stress[7] = eleData(intgrPoint*6+5);
	  }
	} // size 48
	
	else if ( eleData.Size() == 49 ) {
	  
	  f_stress[0] = eleData(intgrPoint*6+1);
	  f_stress[4] = eleData(intgrPoint*6+2);
	  f_stress[8] = eleData(intgrPoint*6+3);
	  f_stress[1] = f_stress[3] = eleData(intgrPoint*6+4);
	  f_stress[2] = f_stress[6] = eleData(intgrPoint*6+5);
	  f_stress[5] = f_stress[7] = eleData(intgrPoint*6+6);
	  
	  
	} // size 49
      } // 8 nodes
    } //non-null response
    
    setPetal(f_stress);
  } // hit


	 
       

}
