# Begin CVS Header 
#   $Source: /Volumes/Home/Users/shoops/cvs/copasi_dev/copasi/libs/libs.pro,v $ 
#   $Revision: 1.8 $ 
#   $Name:  $ 
#   $Author: shoops $ 
#   $Date: 2012/06/15 15:31:50 $ 
# End CVS Header 

# Copyright (C) 2012 - 2010 by Pedro Mendes, Virginia Tech Intellectual 
# Properties, Inc., University of Heidelberg, and The University 
# of Manchester. 
# All rights reserved. 

# Copyright (C) 2008 by Pedro Mendes, Virginia Tech Intellectual 
# Properties, Inc., EML Research, gGmbH, University of Heidelberg, 
# and The University of Manchester. 
# All rights reserved. 

######################################################################
# Automatically generated by qmake (1.07a) Fri Oct 27 13:31:24 2006
######################################################################

TEMPLATE = subdirs

include(../common.pri)

# Directories
SUBDIRS += COPASISE
!contains(BUILD_GUI, no) {
  SUBDIRS += COPASIUI
}

DISTDIRS = $${SUBDIRS}

DISTFILES += \
        libs.pro \
        lib.pri

include(../srcDistribution.pri)

src_distribution.commands = \
  rm -rf ../../copasi_src/copasi/libs; \
  $(CHK_DIR_EXISTS) ../../copasi_src || $(MKDIR) ../../copasi_src; \
  $(CHK_DIR_EXISTS) ../../copasi_src/copasi || $(MKDIR) ../../copasi_src/copasi; \
  $(CHK_DIR_EXISTS) ../../copasi_src/copasi/libs || $(MKDIR) ../../copasi_src/copasi/libs; \
  cp $${DISTFILES} ../../copasi_src/copasi/libs; \
  $$join(DISTDIRS, "; $(MAKE) -f $(MAKEFILE) $@; cd ..; cd ", "cd ", "; $(MAKE) -f $(MAKEFILE) $@; cd ..;")

QMAKE_EXTRA_TARGETS += src_distribution