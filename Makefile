#Makefile for $(HOME)/OpenSees/SRC/VEESlinux

include ../../Makefile.def  #include OpenSees libs

#Linux -s et your path to glut
#INC_DIRS=-I/usr/include -I/glut-3.6/include -I./ -I/usr/local/include/GLUI 
#LIB_DIRS=-L/usr/lib -L/glut-3.6/lib/glut -L/usr/X11R6/lib 


#LDLIBS= -lGL -lGLU -lglut -lXmu -lXt -lSM -lICE -lXext -lX11 \
#	-lXi -lXext -lX11 -lm  -lglui 


OBJS =  Viewer.o GlutSubWindow.o GlutWin.o Camera.o  materials.o auxiliary.o \
	ColorRange.o VisWin.o DrawElement.o DrawSolid.o GaussCoord3D.o \
	DrawSolidDisplaced.o DrawWireFrame.o HashLink.o HashTable.o\
	HideElement.o 	DrawDisplacedWireFrame.o DrawPlaneInABox.o \
	DrawGaussPoints.o DomainViewer.o ElementViewer.o GlyphViewer.o \
	ColorBarViewer.o  DrawStiffness.o ReynoldsGlyph.o

#	 HashLink.o HashTable.o DrawPlaneInABox.o GlyphViwer.o

all: $(OBJS)

tidy:	
	@$(RM) $(RMFLAGS) Makefile.bak *~ #*# core

clean: tidy
	@$(RM) $(RMFLAGS) $(OBJS) *.o

spotless: clean
	@$(RM) $(RMFLAGS)

wipe: spotless
