/**
 * These are the coloring methods for node, element and integration
 * point visulizations. Used by DomainViewer, GUI, and some element viewers.
 * Order by what is hopefully most useful to the user 
 * This separate file is to make the names more global, since they are
 * accessed through ColorRange objects
 */
#ifndef COLOR_SCHEME_NAMES
#define  COLOR_SCHEME_NAMES

static char * colorSchemeNames[] = { "Node Number", "Element Number", 
		    	   "Element Type", "Displacement", 
				    	 "Mean Stress", "Deviatoric Stress",
				     "Min Stiffness", "Stiffness Eigenmode",
				     "DOF 0"};

enum scalarColorSchemes{ BY_NODE_NUM, BY_ELEMENT_NUM, BY_ELEMENT_TYPE, 
			 BY_NODE_DISPLACEMENT, BY_MEAN_STRESS, BY_DEVIATORIC,
			 BY_STIFFNESS, BY_STIFF_EIGENMODE, BY_DOF0 };

#define NUM_COLOR_SCHEMES 9




#endif
