// Begin CVS Header 
//   $Source: /fs/turing/cvs/copasi_dev/copasi/bindings/swig/CMCAMethod.i,v $ 
//   $Revision: 1.2 $ 
//   $Name:  $ 
//   $Author: shoops $ 
//   $Date: 2011/12/19 16:20:17 $ 
// End CVS Header 

// Copyright (C) 2011 by Pedro Mendes, Virginia Tech Intellectual 
// Properties, Inc., University of Heidelberg, and The University 
// of Manchester. 
// All rights reserved. 

%{

#include "steadystate/CMCAMethod.h"

%}

%newobject CMCAMethod::createMethod(CCopasiMethod::SubType subType = CCopasiMethod::mcaMethodReder);


%include "steadystate/CMCAMethod.h"

