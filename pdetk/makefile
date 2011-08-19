include $(DAETK_DIR)/config/$(DAETK_ARCH)
include $(DAETK_DIR)/config/$(DAETK_ARCH).sources
include sources.h

LIBNAME		= libpdetk
ARCHIVE		= $(LIBNAME)$(ARCHIVE_SUFFIX)

#################
# Target Rules
#################

all: $(ARCHIVE)

include ${DAETK_DIR}/config/$(DAETK_ARCH).archive

clean:
	$(RM) $(OBJS) $(FILESTOCLEAN)
	$(RM) -rf ti_files

clobber: clean
	$(RM) $(ARCHIVE)

depend:
	$(MAKEDEPEND) $(MAKEDEPENDFLAGS) -- $(DEFS) $(INCLUDES) -- $(SRCS) -f dep.txt 

.SUFFIXES: .cpp .c .f .o

.cpp.o:
	$(CXX) $(CXXFLAGS) $(INCLUDES) $(DEFS) -c $<

.c.o:
	$(CC) $(CFLAGS) $(INCLUDES) $(DEFS) -c $<

.f.o:
	$(FC) $(FCFLAGS) -c $<

include dep.txt
