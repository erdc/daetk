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
#  reproduce, and distribute this SOFTWARE, and to allow others to do so.
#  The public may copy and use this SOFTWARE, FOR NONCOMMERCIAL USE ONLY,
#  without charge, provided that this Notice and any statement of
#  authorship are reproduced on all copies.  Neither the Government nor
#  the University makes any warranty, express or implied, or assumes any
#  liability or responsibility for the use of this SOFTWARE.
#  
#  If SOFTWARE is modified to produce derivative works, such modified
#  SOFTWARE should be clearly marked, so as not to confuse it with the
#  version available from LANL.
#  
#  For more information about PETE, send e-mail to pete@acl.lanl.gov,
#  or visit the PETE web page at http://www.acl.lanl.gov/pete/.
#  ----------------------------------------------------------------------
# ACL:license

PROJECTROOT = $(shell cd ../..; pwd)
include $(PROJECTROOT)/config/$(PETEARCH)

default: CMRVecTest

operators: $(MAKEEXPR)
	$(MAKEEXPR) --classes CMRVecDefs.in > CMRVecOperators.h

CMRVecTest: CMRVecTest.cpp Chronograph.cpp CMRVec.h CMRVecOperators.h
	$(CXX) -I/home/cekees/codes/CC/daetk -DDEBUG_TRACE -o CMRVecTest $(CXXINC) $(CXXOPTARG) \
	CMRVecTest.cpp Chronograph.cpp CMRVecBlasd.C \
	-lblas

$(MAKEEXPR):
	$(MAKE) -C$(PROJECTROOT) tools

clean:
	rm -f CMRVecTest
	rm -rf ti_files
	rm -rf ii_files

run: CMRVecTest
	CMRVecTest

# ACL:rcsinfo
#  ----------------------------------------------------------------------
#  $RCSfile: makefile,v $   $Author: cekees $
#  $Revision: 1.2 $   $Date: 2002/09/16 16:13:31 $
#  ----------------------------------------------------------------------
# ACL:rcsinfo
