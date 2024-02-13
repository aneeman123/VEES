#ifndef _STRESS_GLYPH__
#define _STRESS_GLYPH_
/**
 * A Reynolds Glyph representing stress, rather than stiffness per se
 * Keeps stress tensor and applies Reynolds glyph methods
 */

#include "ReynoldsGlyph.h"

class StressGlyph:ReynoldsGlyph {
  StressGlyph();
  void  setIntegrPoint( int theNum );
 private:
  double f_stress[9];
}
#endif
