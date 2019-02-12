###############################################################
#
#
# Makefile for buildRateTree 
###############################################################

#everybody here needs root
ROOTFLAGS=`root-config --cflags`

#we also need cdmsbats
CDMSBATSDIR=/home/phys/mast/cdms/processing/cdmsbats
BUILDDIR=$(CDMSBATSDIR)/BUILD
INCDIR=$(BUILDDIR)/include
LIBDIR=$(BUILDDIR)/lib
LINALGINCDIR=$(BUILDDIR)/linalgebra
CDMSBATSINC= -I$(INCDIR) -I$(LINALGINCDIR)/linalgebra -I$(LINALGINCDIR)/boost-numeric/include
CDMSBATSLIB= -L$(LIBDIR)

#CPPFLAGS+= $(ROOTFLAGS)
#CPPFLAGS+= $(CDMSBATSINC)

buildRateTree: buildRateTree.cpp
	g++ -o buildRateTree $(CDMSBATSINC) buildRateTree.cpp `root-config --cflags --glibs` $(CDMSBATSLIB) -lCdmsbats
	
clean:
	rm -f buildRateTree
	rm -f *.o
	rm -f *.so

