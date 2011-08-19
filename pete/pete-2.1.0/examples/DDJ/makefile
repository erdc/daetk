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

default: Listing1 Listing2 Listing3 Listing4

operators: $(MAKEEXPR)
	$(MAKEEXPR) --classes Listing1Defs.in --assign-ops \
	> Listing1Operators.h
	$(MAKEEXPR) --classes Listing2Defs.in --assign-ops \
	> Listing2Operators.h

Listing1: Listing1.cpp Listing1Operators.h
	$(CXX) -o Listing1 $(CXXINC) $(CXXOPTARG) Listing1.cpp -lm

Listing2: Listing2.cpp Listing2Operators.h
	$(CXX) -o Listing2 $(CXXINC) $(CXXOPTARG) Listing2.cpp -lm

Listing3: Listing3.cpp Listing1Operators.h
	$(CXX) -o Listing3 $(CXXINC) $(CXXOPTARG) Listing3.cpp -lm

Listing4: Listing4.cpp Listing2Operators.h
	$(CXX) -o Listing4 $(CXXINC) $(CXXOPTARG) Listing4.cpp -lm

$(MAKEEXPR):
	make -C$(PROJECTROOT) tools

clean:
	rm -f Listing1 Listing2 Listing3 Listing4
	rm -rf ti_files
	rm -rf ii_files

run: Listing1 Listing2 Listing3 Listing4
	Listing1
	Listing2
	Listing3
	Listing4

# ACL:rcsinfo
#  ----------------------------------------------------------------------
#  $RCSfile: makefile,v $   $Author: ckees $
#  $Revision: 1.1 $   $Date: 2001/08/10 15:22:47 $
#  ----------------------------------------------------------------------
# ACL:rcsinfo
