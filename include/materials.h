/*****************************************************************
**
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
** utilities for setting colorscale colors
**/

#ifdef WIN32
#include <windows.h>
#endif

#ifndef _MATERIALS_
#define _MATERIALS_
/* Expecting array with 4 (RGBA) ambient, 4 diffuse, 
 * 4 specular, and 1 shininess qualities
 * set the crayon color, and draw! 
 */

enum scale_type{ LINEAR_COLOR_SCALE, LOG_COLOR_SCALE, POW_COLOR_SCALE };

void set_color(float color[13]);
void set_color_by_param( float heat, float opacity );
void set_color_minus_plus( float color, float opacity );
void set_color_by_param_shiny( float heat, float opacity );
void set_saturation_by_param(float heat, float opacity);
void set_saturation_by_param_shiny(float heat, float opacity);
void set_color_by_param_2D( float heat );
void set_color_by_param_minus_plus( float heat, float opacity );
float get_log_scale_color( float heat, float log_range  );
float get_normal_scale_color( float heat, float linear_range );
void get_color_by_param ( float color[4], float heat );

#endif
