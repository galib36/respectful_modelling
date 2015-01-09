# Begin CVS Header 
#   $Source: /Volumes/Home/Users/shoops/cvs/copasi_dev/copasi/bindings/bindings.pro,v $ 
#   $Revision: 1.12 $ 
#   $Name:  $ 
#   $Author: shoops $ 
#   $Date: 2011/12/19 16:20:20 $ 
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

TEMPLATE = subdirs
DISTFILES = bindings.pro

include(../common.pri)

SUBDIRS += cpp_examples

contains(ENABLE_JAVA,yes){
  SUBDIRS += java
  !isEmpty(JUNIT_PATH){
    SUBDIRS += java/unittests
  }
}

contains(ENABLE_PYTHON,yes){
  SUBDIRS += python
}

contains(ENABLE_OCTAVE,yes){
 SUBDIRS += octave
}

contains(ENABLE_R,yes){
 SUBDIRS += R
}

contains(ENABLE_PERL,yes){
 SUBDIRS += perl
}

contains(ENABLE_CSHARP,yes){
 SUBDIRS += csharp
}


src_distribution.commands = @echo "Bindings are not distributed."


QMAKE_EXTRA_TARGETS += src_distribution
