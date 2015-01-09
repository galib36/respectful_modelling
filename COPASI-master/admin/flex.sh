#!/bin/bash
# Copyright (C) 2010 - 2014 by Pedro Mendes, Virginia Tech Intellectual 
# Properties, Inc., University of Heidelberg, and The University 
# of Manchester. 
# All rights reserved. 

# Copyright (C) 2008 - 2009 by Pedro Mendes, Virginia Tech Intellectual 
# Properties, Inc., EML Research, gGmbH, University of Heidelberg, 
# and The University of Manchester. 
# All rights reserved. 

# Copyright (C) 2005 - 2007 by Pedro Mendes, Virginia Tech Intellectual 
# Properties, Inc. and EML Research, gGmbH. 
# All rights reserved. 


PATH=$PATH:/bin:/usr/bin:/usr/local/bin
AWK=${COPASI_AWK:-gawk}
SED=${COPASI_SED:-sed}
FLEX="$1"

for arg in $@; do
  SOURCE_FILE=$arg
done

PREFIX=`${AWK} -- '$0 ~ "%option prefix=" {split($2, out, "\""); print out[2]}' $SOURCE_FILE`

TARGET_FILE=${SOURCE_FILE/%.*/_lex.cpp}

echo compiling $SOURCE_FILE '==>' $TARGET_FILE

${FLEX}  -t $SOURCE_FILE | \
  ${SED} -e 's?<FlexLexer.h>?"copasi/FlexLexer.h"?' \
         -e 's/include <fstream>/include <iostream>/' \
         -e 's/using std::istream;/using namespace std;/' \
         -e '/using std::ostream;/d' \
         -e '/#include <unistd.h>/d' \
         > $TARGET_FILE

if [ x`uname -a | grep -ic cygwin` = x"1" ]; then
  unix2dos $TARGET_FILE
fi
