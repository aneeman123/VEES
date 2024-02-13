/* ****************************************************************** **
**    OpenSees - Open System for Earthquake Engineering Simulation    **
**          Pacific Earthquake Engineering Research Center            **
**                                                                    **
**                                                                    **
** (C) Copyright 1999, The Regents of the University of California    **
** All Rights Reserved.                                               **
**                                                                    **
** Commercial use of this program without express permission of the   **
** University of California, Berkeley, is strictly prohibited.  See   **
** file 'COPYRIGHT'  in main directory for information on usage and   **
** redistribution,  and for a DISCLAIMER OF ALL WARRANTIES.           **
**                                                                    **
** Developed by:                                                      **
**       Alisa Neeman aneeman@cse.ucsc.edu                                     **
**                                                                    **
**                                                                    **
** ****************************************************************** */


#include <Domain.h>
#include <TclModelBuilder.h>


#include <tcl.h>
#include <tk.h>

#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include "GlutWin.h"
extern ModelBuilder *theBuilder;
extern Domain *theDomain;


int 
TclModelBuilderVeesCommand(ClientData clientData, 
			      Tcl_Interp *interp, int argc,    
			      TCL_Char **argv, 
			      Domain *theDomain, TclModelBuilder *theTclBuilder)
{

   opserr << "VEES!!!\n";
 	GlutWin * g = & GlutWin::Instance();
  	g->createGUI( theDomain);
	return TCL_OK;
}