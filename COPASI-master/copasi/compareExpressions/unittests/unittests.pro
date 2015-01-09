# Begin CVS Header 
#   $Source: /Volumes/Home/Users/shoops/cvs/copasi_dev/copasi/compareExpressions/unittests/unittests.pro,v $ 
#   $Revision: 1.19 $ 
#   $Name:  $ 
#   $Author: shoops $ 
#   $Date: 2011/05/25 15:12:19 $ 
# End CVS Header 

# Copyright (C) 2011 - 2010 by Pedro Mendes, Virginia Tech Intellectual 
# Properties, Inc., University of Heidelberg, and The University 
# of Manchester. 
# All rights reserved. 

# Copyright (C) 2008 by Pedro Mendes, Virginia Tech Intellectual 
# Properties, Inc., EML Research, gGmbH, University of Heidelberg, 
# and The University of Manchester. 
# All rights reserved. 

# Copyright (C) 2001 - 2007 by Pedro Mendes, Virginia Tech Intellectual 
# Properties, Inc. and EML Research, gGmbH. 
# All rights reserved. 

TEMPLATE = app

CONFIG -= qt

include(../../common.pri)

LIBS += -L../../lib/ 


include(../../app.pri)

INCLUDEPATH += ../../..

COPASI_LIBS += $${COPASI_LIBS_SE}

contains(BUILD_OS, WIN32) {
#  LIBS += $$join(COPASI_LIBS, ".lib  ../lib/", ../lib/, .lib)

#  POST_TARGETDEPS += $$join(COPASI_LIBS, ".lib  ../lib/", ../lib/, .lib)
}

contains(BUILD_OS, Linux) {
!isEmpty(CPPUNIT_PATH) {
  LIBS += -L$${CPPUNIT_PATH}/lib -lcppunit ../../sbml/unittests/utilities.o
  INCLUDEPATH += $${CPPUNIT_PATH}/include

  LIBS = $$join(COPASI_LIBS, ".a  ../../lib/lib", ../../lib/lib, .a) \
         $${LIBS}
   
  POST_TARGETDEPS += $$join(COPASI_LIBS, ".a  ../../lib/lib", ../../lib/lib, .a)

  POST_TARGETDEPS += ../../sbml/unittests/utilities.o
}

}

contains(BUILD_OS, SunOS) {
  QMAKE_LFLAGS += -z rescan

!isEmpty(CPPUNIT_PATH) {
  LIBS += -L$${CPPUNIT_PATH}/lib -lcppunit ../../sbml/unittests/utilities.o
  INCLUDEPATH += $${CPPUNIT_PATH}/include

  LIBS = $$join(COPASI_LIBS, ".a  ../../lib/lib", ../../lib/lib, .a) \
         $${LIBS}
   
  POST_TARGETDEPS += $$join(COPASI_LIBS, ".a  ../../lib/lib", ../../lib/lib, .a)
  POST_TARGETDEPS += ../../sbml/unittests/utilities.o
   }
}  

contains(BUILD_OS, Darwin){
  QMAKE_LFLAGS += -Wl,-search_paths_first
  
!isEmpty(CPPUNIT_PATH) {
  LIBS += -L$${CPPUNIT_PATH}/lib -lcppunit ../../sbml/unittests/utilities.o
  INCLUDEPATH += $${CPPUNIT_PATH}/include

  LIBS = $$join(COPASI_LIBS, ".a  ../../lib/lib", ../../lib/lib, .a) \
         $${LIBS}
   
  POST_TARGETDEPS += $$join(COPASI_LIBS, ".a  ../../lib/lib", ../../lib/lib, .a)
  POST_TARGETDEPS += ../../sbml/unittests/utilities.o

}

}


# Input
HEADERS +=  test_depth_first_iterator.h \
            test_expression_comparison.hpp \
            test_node_conversion.hpp \
            test_cnormallogical.hpp \
            test_normalform.hpp \
            test_simplify.h \
            test_compare_utilities.h \
            test_biomodels217_expression.h

SOURCES += test.cpp \
           test_depth_first_iterator.cpp \ 
           test_expression_comparison.cpp \
           test_node_conversion.cpp \
           test_cnormallogical.cpp \
           test_normalform.cpp \
           test_simplify.cpp \
           test_compare_utilities.cpp \
           test_biomodels217_expression.cpp


