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

ifndef PETEARCH
PETEARCH = SGI64KCC
endif

PROJECTROOT = $(shell cd ..; pwd)
include $(PROJECTROOT)/config/$(PETEARCH)

default:
	make -CRGB
	make -CTArray
	make -CVec3
	make -CVector
	make -CMixed
	make -CDDJ

operators:
	make -CRGB operators
	make -CTArray operators
	make -CVec3 operators
	make -CVector operators
	make -CMixed operators
	make -CDDJ operators

run:
	make -CRGB run
	make -CTArray run
	make -CVec3 run
	make -CVector run
	make -CMixed run
	make -CDDJ run

clean:
	make -CRGB clean
	make -CTArray clean
	make -CVec3 clean
	make -CVector clean
	make -CMixed clean
	make -CDDJ clean

# ACL:rcsinfo
#  ----------------------------------------------------------------------
#  $RCSfile: makefile,v $   $Author: ckees $
#  $Revision: 1.1 $   $Date: 2001/08/10 15:22:47 $
#  ----------------------------------------------------------------------
# ACL:rcsinfo
