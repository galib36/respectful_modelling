// Begin CVS Header 
//   $Source: /fs/turing/cvs/copasi_dev/copasi/bindings/swig/CFitProblem.i,v $ 
//   $Revision: 1.3 $ 
//   $Name:  $ 
//   $Author: shoops $ 
//   $Date: 2010/07/16 18:56:27 $ 
// End CVS Header 

// Copyright (C) 2010 by Pedro Mendes, Virginia Tech Intellectual 
// Properties, Inc., University of Heidelberg, and The University 
// of Manchester. 
// All rights reserved. 

// Copyright (C) 2008 by Pedro Mendes, Virginia Tech Intellectual 
// Properties, Inc., EML Research, gGmbH, University of Heidelberg, 
// and The University of Manchester. 
// All rights reserved. 

%{

#include "parameterFitting/CFitProblem.h"

%}

%rename (output) CFitProblem::print;
%ignore CFitProblem::setCallBack;

%include "parameterFitting/CFitProblem.h"



