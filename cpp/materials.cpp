/*******************************************************************
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

#if defined(__APPLE__)&& defined(__MACH__)
#include <GLUT/glut.h>
#include <OpenGL/gl.h>
#include <OpenGL/glu.h>
#else
#include <GL/gl.h>			/* OpenGL header file */
#include <GL/glu.h>			/* OpenGL utilities header file */
#include <GL/glut.h>	
#endif


#include <stdio.h>
#include <math.h>
/* Expecting array with 4 (RGBA) ambient, 4 diffuse, 
 * 4 specular, and 1 shininess qualities
 * set the crayon color, and draw! 
 */
void set_color(float color[13]) {
  float mat_ambient[4];
  float mat_diffuse[4];
  float mat_specular[4];
  int i;

  for( i = 0; i < 4; i++ ) {
    mat_ambient[i] = color[i];
    mat_diffuse[i] = color [i+4];
    mat_specular[i] = color[i+8];
  }

  glMaterialfv(GL_FRONT_AND_BACK, GL_AMBIENT, mat_ambient);
  glMaterialfv(GL_FRONT_AND_BACK, GL_DIFFUSE, mat_diffuse);
  glMaterialfv(GL_FRONT_AND_BACK, GL_SPECULAR, mat_specular);
  glMaterialf(GL_FRONT_AND_BACK, GL_SHININESS, color[12]);

  
} 

/**
 * allows caller to set hue and opacity afterward
 * @param color gets set to an RGB color mapped to zero to 1
 * rainbow color map 
*/
void get_color_by_param ( float color[4], float heat ) {
  int r,g,b;
   
  color[3] = 1;

  r = 0;
  g = 1;
  b  = 2; 


  if (heat <= 0.25 ) {  // blue to cyan; increase green
    color[b] = 1.0;
    color[g] = heat*4.0;
    color[r] = 0.0;
  }

  else if (heat > 0.25 && heat <= 0.5 ){ //cyan to green
	color[g] = 1.0;
	color[r] = 0.0;
    color[b] = 1.0 -( (heat - 0.25)*4);    //reduce blue as we near 0.5
  }

  else if ( heat > 0.5 && heat <= 0.75 ) { //green to yellow
	color[g] = 1.0;
	color[b] = 0.0;
	color[r] = (heat - 0.5)*4; // increase red
	
  }
  else if ( heat > 0.75 && heat <= 1.11 ) { //yellow to red
    if (heat > 1.0 )
      heat = 1.0; //clamp
	color[r] = 1.0;
	color[b] = 0.0;
	color[g] = 1.0 -( (heat - 0.75)*4); //decrease green
  }
  else {
    printf("Error setting color in get_color_by_param; heat %f\n", heat);
    return;
  }

}


/* 
 * create a custom rgb color where pure red is hot and pure
 * blue is cold. Function will set ambient, diffuse and specular for you.
 * | pure blue | blue + green | pure green | red + green | pure red |
 * !pre: both float and heat are in the range of 0-1
 */
void set_color_by_param( float heat, float opacity ) {
  float color[4]; /* r,g,b,a */
  int r,g,b;
 
 float spec_color[4] = {0.4,0.4,0.4,1.0};//{1.0,1.0,1.0,1.0};
   
  color[3] = opacity;

  r = 0;
  g = 1;
  b  = 2; 
  

  if (heat <= 0.25 ) {  // blue to cyan; increase green
    color[b] = 1.0;
    color[g] = heat*4.0;
    color[r] = 0.0;
  }

  else if (heat > 0.25 && heat <= 0.5 ){ //cyan to green
	color[g] = 1.0;
	color[r] = 0.0;
    color[b] = 1.0 -( (heat - 0.25)*4);    //reduce blue as we near 0.5
  }

  else if ( heat > 0.5 && heat <= 0.75 ) { //green to yellow
	color[g] = 1.0;
	color[b] = 0.0;
	color[r] = (heat - 0.5)*4; // increase red
	
  }
  else if ( heat > 0.75 && heat <= 1.11 ) { //yellow to red
    if (heat > 1.0 )
      heat = 1.0; //clamp
	color[r] = 1.0;
	color[b] = 0.0;
	color[g] = 1.0 -( (heat - 0.75)*4); //decrease green
  }
  else {
    printf("Error setting color in set_color_by_param; heat %f\n", heat);
    return;
  }


  glMaterialfv( GL_FRONT_AND_BACK, GL_AMBIENT, color );
  glMaterialfv( GL_FRONT_AND_BACK, GL_DIFFUSE, color );

}



/**
 * create a custom rgb color where pure red is hot and pure
 * blue is cold. Function will set ambient, diffuse and specular for you.
 * | pure blue | blue + green | pure green | red + green | pure red |
 * !pre:color is in the range of -1 to +1, opacity 0-1
 */
void set_color_minus_plus( float heat, float opacity ) {
  float color[4]; /* r,g,b,a */
  int r,g,b;
 
 float spec_color[4] = {0.4,0.4,0.4,1.0};//{1.0,1.0,1.0,1.0};
   
  color[3] = opacity;

  r = 0;
  g = 1;
  b  = 2; 
  

  if ( heat >= -1.0 && heat <= -0.5 ) {  // blue to cyan; increase green
    color[b] = 1.0;
    color[g] = (1.0 + heat)*2; //-1 to -0.5
    color[r] = 0.0;
  }

  else if (heat > -0.5 && heat <= 0.0 ){ //cyan to green
	color[g] = 1.0;
	color[r] = 0.0;
    color[b] = 1.0 -( heat*2 );    //reduce blue as we near 0.0
  }

  else if ( heat > 0.0 && heat <= 0.5 ) { //green to yellow
	color[g] = 1.0;
	color[b] = 0.0;
	color[r] = (heat *2); // increase red
	
  }
  else if ( heat > 0.5 && heat <= 1.0 ) { //yellow to red
	color[r] = 1.0;
	color[b] = 0.0;
	color[g] = (1.0 - heat)*2; //decrease green
  }
  else {
    // printf("Error setting color in set_color_by_param; heat %f\n", heat);
    // clamp
    if (heat > 1.0) {
      	color[r] = 1.0;
	color[b] = 0.0;
	color[g] = 0.0;
    }
    if (heat < -1.0) {
      color[b] = 1.0;
      color[g] = 0.0;	
      color[r] = 0.0;    
    }
  }

  glMaterialfv( GL_FRONT_AND_BACK, GL_AMBIENT, color );
  glMaterialfv( GL_FRONT_AND_BACK, GL_DIFFUSE, color );
}


/* 
 * same as set_color_by_param but with specular and shininess 
 * components
 *
 */
void set_color_by_param_shiny( float heat, float opacity ) {
  float color[4]; /* r,g,b,a */
  int r,g,b;
  float shiny = 120; //very shiny
  float spec_color[4] = {0.4,0.4,0.4,1.0};//{1.0,1.0,1.0,1.0};
   
  color[3] = opacity;

  r = 0;
  g = 1;
  b  = 2; 
  

  if (heat <= 0.25 ) {  // blue to cyan; increase green
    color[b] = 1.0;
    color[g] = heat*4.0;
    color[r] = 0.0;
  }

  else if (heat > 0.25 && heat <= 0.5 ){ //cyan to green
	color[g] = 1.0;
	color[r] = 0.0;
    color[b] = 1.0 -( (heat - 0.25)*4);    //reduce blue as we near 0.5
  }

  else if ( heat > 0.5 && heat <= 0.75 ) { //green to yellow
	color[g] = 1.0;
	color[b] = 0.0;
	color[r] = (heat - 0.5)*4; // increase red
	
  }
  else if ( heat > 0.75 && heat <= 1.0 ) { //yellow to red
	color[r] = 1.0;
	color[b] = 0.0;
	color[g] = 1.0 -( (heat - 0.75)*4); //decrease green
  }
  else {
    printf("Error setting color in set_color_by_param; heat %f\n", heat);
    return;
  }


  glMaterialfv( GL_FRONT_AND_BACK, GL_AMBIENT, color );
  glMaterialfv( GL_FRONT_AND_BACK, GL_DIFFUSE, color );
  glMaterialfv( GL_FRONT_AND_BACK, GL_SPECULAR, spec_color );
  glMaterialfv( GL_FRONT_AND_BACK, GL_SHININESS, &shiny );

}

/**
 * 0 is black, 1 is white */
void set_saturation_by_param_shiny(float heat, float opacity) {

  float shiny = 30; //very shiny
  float spec_color[4] = {0.1,0.1,0.1,1.0};//{1.0,1.0,1.0,1.0};

  float color[4]; /* r,g,b,a */
  color[0] = color[1] = color[2] = heat;
  color[3] = opacity;

  glMaterialfv( GL_FRONT_AND_BACK, GL_AMBIENT, color );
  glMaterialfv( GL_FRONT_AND_BACK, GL_DIFFUSE, color );
  glMaterialfv( GL_FRONT_AND_BACK, GL_SPECULAR, spec_color );
  glMaterialfv( GL_FRONT_AND_BACK, GL_SHININESS, &shiny );
}

/**
 * 0 is black, 1 is white */
void set_saturation_by_param(float heat, float opacity) {
  float color[4]; /* r,g,b,a */
  color[0] = color[1] = color[2] = heat;
  color[3] = opacity;

  glMaterialfv( GL_FRONT_AND_BACK, GL_AMBIENT, color );
  glMaterialfv( GL_FRONT_AND_BACK, GL_DIFFUSE, color );
}



/* 
 * get a color value that is the log of the current log, 
 * shifted over 0-1 and normalized
 */
float get_log_scale_color( float heat, float log_range  )
{    
     float log_heat;
        

     /* get in +n to -n range */
     if(heat > -1 && heat < 1) { 
       log_heat = 0.0;
       return log_heat;
     }
     else {
       if ( heat < 0.0 )
	 log_heat = -1*(log(fabs( heat )));
       else 
	 log_heat = log( heat ); 
     }

     
     log_heat = log_heat + log_range; /* shift over 0-1 */
     log_heat = log_heat/(log_range*2); /* normalize */
     return log_heat;
}

/* 
 * linear range of values shifted into 0-1 where 0.5 maps to 0
 */
float get_normal_scale_color( float heat, float linear_range ) 
{
  float normal_heat;
  
  normal_heat = heat + linear_range; /* shift over 0-1 */
  normal_heat = normal_heat/(linear_range*2); /* normalize */
  return normal_heat;
} 


/* 
 * create a custom rgb color where pure red is hot and pure
 * blue is cold. Function will set ambient, diffuse and specular for you.
 * | pure blue | blue + green | pure green | red + green | pure red |
 * !pre: both float and heat are in the range of -1 to +1
 */
void set_color_by_param_minus_plus( float heat, float opacity ) {
  float color[4]; /* r,g,b,a */
  float scale;
  int r,g,b;
  float spec_color[4] = {1.0,1.0,1.0,1.0};

  color[3] = opacity;
  r = 0;
  g = 1;
  b  = 2; 
  if (heat <= 0.0 ) {      
    scale = fabs(heat);
    color[b] = scale;
    color[g] = 1.0 - scale;
    color[r] = 0.0;
  }
  else if (heat > 0.0 && heat <= 1.0 ){
    scale = heat;
    color[r] = scale;
    color[g] = 1.0 - scale;
    color[b] = 0.0;
  }
  else {
    printf("Error setting color in set_color_by_param_minus_plus; heat %f\n", 
	   heat);
    return;
  }
  glMaterialfv(GL_FRONT_AND_BACK, GL_AMBIENT, color );
  glMaterialfv(GL_FRONT_AND_BACK, GL_DIFFUSE, color);
  glMaterialfv(GL_FRONT_AND_BACK, GL_SPECULAR, spec_color);

}




/* In glOrtho2D colors are flat, not "materials"
 * with ambient and diffuse qualities. This function sets
 * unshaded colors from 0 to 1 (cold to hot)
 */
void set_color_by_param_2D( float heat ) {
  float color[3]; /* r,g,b */
  float scale;
  int r,g,b;
 
  r = 0;
  g = 1;
  b  = 2; 

  if (heat <= 0.5 ) {      
    scale = heat*2.0;
    color[b] = 1.0 - scale;
    color[g] = scale;
    color[r] = 0.0;
  }
  else if (heat > 0.5 && heat <= 1.0 ){
    scale = (heat - 0.5)*2;
    color[r] = scale;
    color[g] = 1.0 - scale;
    color[b] = 0.0;
  }
  else {
    printf("Error setting color in set_color_by_param; heat %f\n", heat);
    return;
  }
  glColor3f(color[0], color[1],color[2]);

}
