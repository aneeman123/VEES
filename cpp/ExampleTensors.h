#ifndef EXAMPLE_TENSORS
#define EXAMPLE_TENSORS

//cross anisotropic stiffness matrix
static double cross[36] =
  {   14.8226,    4.5123,    0.2378,         0,         0,         0,
	  4.5123,   14.8226,    0.2378,         0,         0,         0,
	  0.2378,    0.2378,    0.5559,         0,         0,         0,
	  0,         0,         0,    1.6600,         0,         0,
	  0,         0,         0,         0,    1.6600,         0,
	  0,         0,         0,         0,         0,    6.6225 };


static double elasticIso[36] = //guanzhou elastic isotropic material
  {  //no normalizing constants
	27925925.925926, 15037037.037037, 15037037.037037, 0.000000, 0.000000, 0.000000,
	15037037.037037, 27925925.925926, 15037037.037037, 0.000000, 0.000000, 0.000000,
	15037037.037037, 15037037.037037, 27925925.925926, 0.000000, 0.000000, 0.000000,
	0.000000, 0.000000, 0.000000, 51555554, 0.000000, 0.000000,
	0.000000, 0.000000, 0.000000, 0.000000, 51555554, 0.000000,
	0.000000, 0.000000, 0.000000, 0.000000, 0.000000,51555554
  };



//ideal isotropy from Basser paper
static double isoSet[6][9] =
  {
    {-0.8, 0.0, 0.0, 0.0, 0.4, 0.0, 0.0, 0.0, 0.4},
    {0.577, 0.000, 0.000, 0.000, 0.577, 0.000, 0.000, 0.000, 0.577},
    {0.0, 0.0, 0.0, 0.0, -0.7, 0.0, 0.0, 0.0, 0.7},
    {0.0, 1.0, 0.0, 1.0, 0.0, 0.0, 0.0, 0.0, 0.0 },
    {0.0, 0.0, 0.0, 0.0, 0.0, 1.0, 0.0, 1.0, 0.0 },
    {0.0, 0.0, 1.0, 0.0, 0.0, 0.0, 1.0, 0.0, 0.0 }
    
  };

static double paperIsoSet[6][9] = //normalized
  {
	{0.5773,0,0,0,0.5773,0,0,0,0.5773},
	{.4082,0,0,0,-.8164,0,0,0,.4082},
	{.7071,0,0,0,0,0,0,0,-.7071},
	{0,.7071,0,.7071,0,0,0,0,0},
	{0,0,0,0,0,.7071,0,.7071,0},
	{0,0,.7071,0,0,0,.7071,0,0}
  };
static double identityMatrix[36]={
  1,0,0,0,0,0,
  0,1,0,0,0,0,
  0,0,1,0,0,0,
  0,0,0,1,0,0,
  0,0,0,0,1,0,
  0,0,0,0,0,1
};


//cross anisotropic stiffness tensor "natural unrolling" - glyph beautifully wrong
// isotropic part looks ok, but other parts appear off by 90 degrees on 1 axis
/*
static double crossSet[9][9] =  // natural unrolling of 3^4 tensor: do not use
  {{0,0,0, 0,0,.7071, 0,-0.7071,0},
   {0, .6568,.2618,.6568,0,0,-0.2618,0,0},
   {0, -.2618,0.6568, 0.2618,0,0,-0.6568,0,0},
   {-0.0127,0,0,0,-0.0127,0,0,0,0.9998},
   {0,0,0,0,0,-0.7071, 0,-0.7071, 0},
   {0,0.7071,0,0.7071,0,0,0,0,0},
   {0.7071,0,0,0,-0.7071,0,0,0,0},
   {0,0,-0.7071, 0,0,0,-0.7071,0,0},
   {-0.7070,0,0,0,-0.7071,0,0,0,-0.0179}};
*/

static double tCrossSet[6][9] = //tarantalo eigentensors
  {{-0.0127,0,0,0,-0.0127,0,0,0,0.9998},
   {0,1,0,1,0,0,0,0,0},
   {-0.7071,0,0,0,0.7071,0,0,0,0},
   {0,0,1,0,0,0,1,0,0},
   {-0.7071,0,0,0,-0.7071,0,0,0,-0.0179},
   {0,0,0,0,0,1,0,1,0}
  };

//cubic symmetry (copper) from Sutcliffe paper
static double copper[36] =
  {17,  12.3, 12.3, 0, 0, 0,
   12.3,17,   12.3, 0, 0, 0,
   12.3,12.3, 17,   0, 0, 0,
   0,   0,    0,    7.5, 0,0,
   0,   0,    0,   0,  7.5, 0,
   0,   0,    0,   0,  0,  7.5
  };
static double copperSet[6][9] = //cubic
  {{0.5774, 0,0,0,0.5774,0,0,0,0.5774}, 
   {0.0017, 0,0, 0,0.7063, 0, 0,0,-0.7079},
   {-0.8165, 0,0, 0, 0.4097, 0, 0, 0, 0.4068},
   {0,0.7071,0,0.7071,0,0,0,0,0},
   {0,0,0,0,0,0.7071,0,0.7071,0},
   {0,0,0.7071,0,0,0,0.7071,0,0} 
  };

//tetragonal (Sutcliffe)
static double tin[36] = {
  8.391, 4.870, 2.810, 0, 0, 0,
  4.870, 8.381, 2.810, 0, 0, 0,
  2.810, 2.810, 9.665, 0, 0, 0,
  0,     0,     0, 1.754, 0, 0,
  0,     0,     0, 0, 1.754, 0,
  0,     0,     0, 0, 0, 0.7407
};

static double tinSet[6][9] =
  {{0.5492, 0 ,0, 0,0.5492, 0, 0, 0,  0.5492},
   {-0.3834, 0, 0, 0,-0.3834, 0, 0,0, 0.8403 },
   {-0.7071, 0, 0, 0, 0.7071, 0, 0,0,0},

   {0,0,0.7071,0,0,0,0.7071,0,0},
   {0,0,0,0,0,0.7071,0,0.7071,0},
   {0,0.7071,0,0.7071,0,0,0,0,0}

   };


static double stretchCrossSet[6][9] = //tarantalo eigentensors*eigenvalues
  {{-0.00698373,0,0,0,-0.00698373,0,0,0,0.54979}, //0.5499 eigenvalue
   {0,3.32,0,3.32,0,0,0,0,0},                  //3.32
   {-7.290413,0,0,0,7.290413,0,0,0,0},       //10.3103
   {0,0,13.245,0,0,0,13.245,0,0},                  //13.245
   {-13.67595,0,0,0,-13.67595,0,0,0,-0.346202},//19.3409
   {0,0,0,0,0,3.32 ,0,3.32,0}                   //3.32
  };


//shear stress
static double stress[9] = 
  {0, 1, 1,
   1, 0, 1,
   1, 1, 0  };

//diagonalized stress
static double stress2[9] = 
  {1,0,0,
   0,2,0,
  0,0,3};

//isotropic stress
static double iso_stress[9] = 
  {1,0,0,
   0,1,0,
  0,0,1};



#endif
