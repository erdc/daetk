# ACL:license
#  ----------------------------------------------------------------------
#  This software and ancillary information (herein called "SOFTWARE")
#  called PETE (Portable Expression Template Engine) is
#  made available under the terms described here.  The SOFTWARE has been
#  approved for release with associated LA-CC Number LA-CC-99-5.
#  
#  Unless otherwise indicated, this SOFTWARE has been authored by an
#  employee or employees of the University of California, operator of the
#  Los Alamos National Laboratory under Contract No.  W-7405-ENG-36 with
#  the U.S. Department of Energy.  The U.S. Government has rights to use,
#  reproduce, and distribute this SOFTWARE. The public may copy, distribute,
#  prepare derivative works and publicly display this SOFTWARE without 
#  charge, provided that this Notice and any statement of authorship are 
#  reproduced on all copies.  Neither the Government nor the University 
#  makes any warranty, express or implied, or assumes any liability or 
#  responsibility for the use of this SOFTWARE.
#  
#  If SOFTWARE is modified to produce derivative works, such modified
#  SOFTWARE should be clearly marked, so as not to confuse it with the
#  version available from LANL.
#  
#  For more information about PETE, send e-mail to pete@acl.lanl.gov,
#  or visit the PETE web page at http://www.acl.lanl.gov/pete/.
#  ----------------------------------------------------------------------
# ACL:license

ifndef PETEARCH
PETEARCH = SGI64KCC
endif

PROJECTROOT = $(shell cd ../..; pwd)
include $(PROJECTROOT)/config/$(PETEARCH)

default: RGB

operators: $(MAKEEXPR)
	$(MAKEEXPR) --classes RGBDefs.in \
	--guard RGBOPS_H > RGBOperators.h

RGB: RGB.cpp RGB.h RGBOperators.h
	$(CXX) -o RGB $(CXXINC) $(CXXOPTARG) RGB.cpp 

$(MAKEEXPR):
	make -C$(PROJECTROOT) tools

clean:
	rm -f RGB
	rm -rf ti_files
	rm -rf ii_files

run: RGB
	RGB

# ACL:rcsinfo
#  ----------------------------------------------------------------------
#  $RCSfile: makefile,v $   $Author: ckees $
#  $Revision: 1.1 $   $Date: 2001/08/10 15:22:47 $
#  ----------------------------------------------------------------------
# ACL:rcsinfo
