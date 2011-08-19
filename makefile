include $(DAETK_DIR)/config/$(DAETK_ARCH)
include $(DAETK_DIR)/config/$(DAETK_ARCH).sources
include sources.h

LIBNAME		= libdaetk
ARCHIVE		= $(LIBNAME)$(ARCHIVE_SUFFIX)

#################
# Target Rules
#################

all: $(ARCHIVE)

include ${DAETK_DIR}/config/$(DAETK_ARCH).archive

install: all
	cp -f $(ARCHIVE) ${PREFIX}/lib
	cp -f *.h ${PREFIX}/include
	cp -rf pete/pete-2.1.0/src/PETE	${PREFIX}/include

clean:
	$(RM) $(OBJS) $(FILESTOCLEAN)
	$(RM) -rf ti_files

clobber: clean
	$(RM) $(ARCHIVE)

depend:
	$(MAKEDEPEND) $(MAKEDEPENDFLAGS) -- $(CXXFLAGS) $(CCFLAGS) $(FCFLAGS) $(INCLUDES) -- $(SRCS) -f dep.txt 

.SUFFIXES: .cpp .c .f .o

.cpp.o:
	$(CXX) $(CXXFLAGS) $(INCLUDES) $(DEFS) -c $<

.c.o:
	$(CC) $(CFLAGS) $(INCLUDES) $(DEFS) -c $<

.f.o:
	$(FC) $(FCFLAGS) -c $<

include dep.txt
