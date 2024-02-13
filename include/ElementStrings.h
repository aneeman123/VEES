#include <classTags.h>
/** All the tags so that we can get the names and put them in a drop down 
 * list so user can show or hide them 
 * Note: there are some doubles of classtag numbers in here - I am ignoring the
 * issue for now, just following OpenSees specs.
*/


#ifndef _ELEMENT_STRINGS_
#define _ELEMENT_STRINGS_

#define numEleTypes  37 //array length of eleStrings

static char * eleStrings[] = {
  "BbarBrick",
  "beam2d02",
  "beam2d03",
  "beam2d04",
  "beam3d01",
  "beam3d02",
  "BeamWithHinges2d",
  "BeamWithHinges3d",
  "Brick",
  "ConstantPressureVolumeQuad",
  "CorotTruss",
  "CorotTrussSection",
  "DispBeamColumn2d",
  "DispBeamColumn3d",
  "EightNodeBrick",
  "EightNodeBrick_u_p_U",
  "ElasticBeam2d",
  "ElasticBeam3d",
  "EnhancedQuad",
  "fElmt02",
  "fElmt05",
  "ForceBeamColumn2d",
  "ForceBeamColumn3d",
  "FourNodeQuad",
  "FourNodeQuadUP",
  "Joint2D",
  "NineNodeMixedQuad",
  "NLBeamColumn2d",
  "NLBeamColumn3d",
  "ShellMITC4",
  "Subdomain",
  "Truss",
  "TrussSection",
  "TwentyNodeBrick",
  "TwentyNodeBrick_u_p_U",
  "ZeroLength",
    "ZeroLengthSection" 
}; 

static int eleTags[] = {
  ELE_TAG_BbarBrick,
  ELE_TAG_beam2d02,
  ELE_TAG_beam2d03,
  ELE_TAG_beam2d04,
  ELE_TAG_beam3d01,
  ELE_TAG_beam3d02,
  ELE_TAG_BeamWithHinges2d,
  ELE_TAG_BeamWithHinges3d,
  ELE_TAG_Brick,
  ELE_TAG_ConstantPressureVolumeQuad,
  ELE_TAG_CorotTruss,
  ELE_TAG_CorotTrussSection,
  ELE_TAG_DispBeamColumn2d,
  ELE_TAG_DispBeamColumn3d,
  ELE_TAG_EightNodeBrick,
  ELE_TAG_EightNodeBrick_u_p_U,
  ELE_TAG_ElasticBeam2d,
  ELE_TAG_ElasticBeam3d,
  ELE_TAG_EnhancedQuad,
  ELE_TAG_fElmt02,
  ELE_TAG_fElmt05,
  ELE_TAG_ForceBeamColumn2d,
  ELE_TAG_ForceBeamColumn3d,
  ELE_TAG_FourNodeQuad,
  ELE_TAG_FourNodeQuadUP,
  ELE_TAG_Joint2D,
  ELE_TAG_NineNodeMixedQuad,
  ELE_TAG_NLBeamColumn2d,
  ELE_TAG_NLBeamColumn3d,
  ELE_TAG_ShellMITC4,
  ELE_TAG_Subdomain,
  ELE_TAG_Truss,
  ELE_TAG_TrussSection,
  ELE_TAG_TwentyNodeBrick,
  ELE_TAG_TwentyNodeBrick_u_p_U,
  ELE_TAG_ZeroLength,
  ELE_TAG_ZeroLengthSection 
};


#endif
