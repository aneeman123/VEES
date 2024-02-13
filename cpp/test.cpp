/* a small driver to test xml and opensees together */

#include <Domain.h>
#include <OPS_Globals.h>
#include <OPS_Stream.h>
#include <stdio.h>
#include <StandardStream.h>

#include <ArrayOfTaggedObjects.h>

// includes for the domain classes
#include <Domain.h>
#include <Node.h>
#include <Truss.h>
#include <ElasticMaterial.h>
#include <SP_Constraint.h>
#include <LoadPattern.h>
#include <LinearSeries.h>
#include <NodalLoad.h>

// includes for the analysis classes
#include <StaticAnalysis.h>
#include <AnalysisModel.h>
#include <Linear.h>
#include <PenaltyConstraintHandler.h>
#include <DOF_Numberer.h>
#include <RCM.h>
#include <LoadControl.h>
#include <BandSPDLinSOE.h>
#include <BandSPDLinLapackSolver.h>

#include <xercesc/dom/DOM.hpp>
#include <xercesc/parsers/XercesDOMParser.hpp>
#include "DomainBuilder.h"

XERCES_CPP_NAMESPACE_USE

#if defined(XERCES_NEW_IOSTREAMS)
#include <iostream>
#else
#include <iostream.h>
#endif

double        ops_Dt = 0;
Domain       *ops_TheActiveDomain = 0;
Element      *ops_TheActiveElement = 0;

#ifndef WIN32
StandardStream sserr;
OPS_Stream opserr =&sserr;
#endif

#ifdef WIN32
StandardStream sserr;
OPS_Stream opserr;
#endif


int init(Domain * theDomain)
{


    //	now create a domain and a modelbuilder
    //  and build the model
    //

 
    
    // create the nodes using constructor: 
    //		Node(tag, ndof, crd1, crd2)
    // and then add them to the domain
    
    Node *node1 = new Node(1, 2,   0.0,  0.0);
    Node *node2 = new Node(2, 2, 144.0,  0.0);
    Node *node3 = new Node(3, 2, 168.0,  0.0);    
    Node *node4 = new Node(4, 2,  72.0, 96.0);        
    theDomain->addNode(node1);
    theDomain->addNode(node2);
    theDomain->addNode(node3);
    theDomain->addNode(node4);
    printf("added nodes\n");
	/*
	node1->printXML(sserr);
	node2->printXML(sserr);
	node3->printXML(sserr);
	node4->printXML(sserr);
	*/
    // create an elastic material using constriuctor:  
    //		ElasticMaterialModel(tag, E)

    UniaxialMaterial *theMaterial = new ElasticMaterial(1, 3000);
    
    // create the truss elements using constructor:
    //		Truss(tag, dim, nd1, nd2, Material &,A)
    // and then add them to the domain
    
    Truss *truss1 = new Truss(1, 2, 1, 4, *theMaterial, 10.0);
    Truss *truss2 = new Truss(2, 2, 2, 4, *theMaterial,  5.0);    
    Truss *truss3 = new Truss(3, 2, 3, 4, *theMaterial,  5.0);        
    theDomain->addElement(truss1);
    theDomain->addElement(truss2);
    theDomain->addElement(truss3);    
    printf("added elements\n");
    // create the single-point constraint objects using constructor:
    //		SP_Constraint(tag, nodeTag, dofID, value)
    // and then add them to the domain
    
    SP_Constraint *sp1 = new SP_Constraint(1, 1, 0, 0.0);
    SP_Constraint *sp2 = new SP_Constraint(2, 1, 1, 0.0);    
    SP_Constraint *sp3 = new SP_Constraint(3, 2, 0, 0.0);
    SP_Constraint *sp4 = new SP_Constraint(4, 2, 1, 0.0);    
    SP_Constraint *sp5 = new SP_Constraint(5, 3, 0, 0.0);
    SP_Constraint *sp6 = new SP_Constraint(6, 3, 1, 0.0);        
    theDomain->addSP_Constraint(sp1);
    theDomain->addSP_Constraint(sp2);
    theDomain->addSP_Constraint(sp3);
    theDomain->addSP_Constraint(sp4);    
    theDomain->addSP_Constraint(sp5);    
    theDomain->addSP_Constraint(sp6);    

    // construct a linear time series object using constructor:
    //		LinearSeries()
    printf("added constraints\n");
    TimeSeries *theSeries = new LinearSeries();
    
    // construct a load pattren using constructor:
    //		LoadPattern(tag)
    // and then set it's TimeSeries and add it to the domain
    
    LoadPattern *theLoadPattern = new LoadPattern(1);
    theLoadPattern->setTimeSeries(theSeries);
    theDomain->addLoadPattern(theLoadPattern);
    
    // construct a nodal load using constructor:
    //		NodalLoad(tag, nodeID, Vector &)
    // first construct a Vector of size 2 and set the values NOTE C 
	//INDEXING
    // then construct the load and add it to the domain
    
    Vector theLoadValues(2);
    theLoadValues(0) = 100.0;
    theLoadValues(1) = -50.0;
    NodalLoad *theLoad = new NodalLoad(1, 4, theLoadValues);
    theDomain->addNodalLoad(theLoad, 1);

    // create an Analysis object to perform a static analysis of the model
    //  - constructs:
    //    AnalysisModel of type AnalysisModel,
    //	  EquiSolnAlgo of type Linear
    //	  StaticIntegrator of type LoadControl
    //	  ConstraintHandler of type Penalty
    //    DOF_Numberer which uses RCM
    //    LinearSOE of type Band SPD
    // and then printf('the StaticAnalysis object
    printf("abaout to analyze");
    AnalysisModel     *theModel = new AnalysisModel();
    EquiSolnAlgo      *theSolnAlgo = new Linear();
    StaticIntegrator  *theIntegrator = new LoadControl(1.0, 1, 1.0, 1.0);
    ConstraintHandler *theHandler = new 
PenaltyConstraintHandler(1.0e8,1.0e8);
    RCM               *theRCM = new RCM();
    DOF_Numberer      *theNumberer = new DOF_Numberer(*theRCM);    
    BandSPDLinSolver  *theSolver = new BandSPDLinLapackSolver();       
    LinearSOE         *theSOE = new BandSPDLinSOE(*theSolver);        

    StaticAnalysis    theAnalysis( *theDomain,
				  *theHandler,
				  *theNumberer,
				  *theModel,
				  *theSolnAlgo,
				  *theSOE,
				  *theIntegrator);

    // perform the analysis & print out the results for the domain
    int numSteps = 1;
	theAnalysis.analyze(numSteps);

	printf("analysis complete\n");
    //sserr << *theDomain;

	return 0;
}



// boilerplate DOM loading example.
int main(int argc, char *argv[])
{
  bool bFailed = false;
  printf("test..\n");
  if (argc < 2) {
	printf("usage: ./test <xml file>\n");
   return -1;
  }
  // initialize the XML library.
  //double check that file name ends in .xml
	Domain * d;
  DomainBuilder db;
  d = db.parseDomain( argv[1] );
  //  d->Print( opserr , 0);
  printf("done!\n");

	d = new Domain();
	init( d );

 
  return 0;
}
