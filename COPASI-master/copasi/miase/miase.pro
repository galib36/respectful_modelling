# Begin CVS Header 
#   $Source: /Volumes/Home/Users/shoops/cvs/copasi_dev/copasi/miase/miase.pro,v $ 
#   $Revision: 1.6 $ 
#   $Name:  $ 
#   $Author: shoops $ 
#   $Date: 2011/03/07 19:30:22 $ 
# End CVS Header 

# Copyright (C) 2011 - 2010 by Pedro Mendes, Virginia Tech Intellectual 
# Properties, Inc., University of Heidelberg, and The University 
# of Manchester. 
# All rights reserved. 

# Copyright (C) 2008 by Pedro Mendes, Virginia Tech Intellectual 
# Properties, Inc., EML Research, gGmbH, University of Heidelberg, 
# and The University of Manchester. 
# All rights reserved. 

LIB = miase
SRC_TARGET = miase
DISTFILES = miase.pro

# Input
HEADERS +=  CMiaseParser.h \
			CMiaseBase.h \
			CMiaseSimulation.h \
			CMiaseModel.h \
			CMiaseTask.h \
			CMiaseOutput.h \
			CMiaseML.h
           
           
SOURCES +=  CMiaseParser.cpp \
			CMiaseBase.cpp \
			CMiaseSimulation.cpp \
			CMiaseModel.cpp \
			CMiaseTask.cpp \
			CMiaseOutput.cpp \
			CMiaseML.cpp

include(../lib.pri)
include(../common.pri)
include(../srcDistribution.pri)
