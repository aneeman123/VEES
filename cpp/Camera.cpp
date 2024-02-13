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
**
** Data needed for basic camera, and associated methods.
** camera does zoom, but object must rotate itself
**/

#include "Camera.h"

Camera::Camera() {
  NO_ZOOM = 200; //default position of camera
  eye[X] = 0;
  eye[Y] = 0;
  eye[Z] = NO_ZOOM;
  at[X] = at[Y] = at[Z] = 0; // look at origin
  up[X] = 0;
  up[Y] = 1; //+Y is up
  up[Z] = 0;

}

/* VisWin should z length based on 
  size of volume to be viewed */
void Camera::initCamera( int z_length ) {
  NO_ZOOM  = z_length;
  
}

void Camera::lookAt() {
 gluLookAt( eye[X],eye[Y],eye[Z],
            at[X],at[Y],at[Z],
			up[X],up[Y], up[Z]);
}

void Camera::calculateCameraPosition() {
  eye[Z] = NO_ZOOM - curZoom;
}

void Camera::setZoom( double zoom ) {
  curZoom = zoom;
}

/** so orthographic projection can also have zoom */  
double Camera::getZoom() {
  return curZoom;
}

void Camera::zoomIn() {
  curZoom +=3;
}

void Camera::zoomOut() {
  curZoom -=3;
}
