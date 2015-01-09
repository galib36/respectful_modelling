// Begin CVS Header
//   $Source: /Volumes/Home/Users/shoops/cvs/copasi_dev/copasi/compareExpressions/unittests/test.cpp,v $
//   $Revision: 1.18 $
//   $Name:  $
//   $Author: shoops $
//   $Date: 2011/03/07 19:26:44 $
// End CVS Header

// Copyright (C) 2011 - 2010 by Pedro Mendes, Virginia Tech Intellectual
// Properties, Inc., University of Heidelberg, and The University
// of Manchester.
// All rights reserved.

// Copyright (C) 2008 by Pedro Mendes, Virginia Tech Intellectual
// Properties, Inc., EML Research, gGmbH, University of Heidelberg,
// and The University of Manchester.
// All rights reserved.

// Copyright (C) 2001 - 2007 by Pedro Mendes, Virginia Tech Intellectual
// Properties, Inc. and EML Research, gGmbH.
// All rights reserved.
#define COPASI_MAIN

#include <cppunit/ui/text/TestRunner.h>

#define COPASI_MAIN

#include "copasi.h"

#include "test_depth_first_iterator.h"
#include "test_node_conversion.hpp"
#include "test_expression_comparison.hpp"
#include "test_normalform.hpp"
#include "test_cnormallogical.hpp"
#include "test_simplify.h"
#include "test_compare_utilities.h"
#include "test_biomodels217_expression.h"

int main(int /*argc*/, char ** /*argv*/)
{
  CppUnit::TextUi::TestRunner runner;
  runner.addTest(test_depth_first_iterator::suite());
  runner.addTest(test_cnormallogical::suite());
  runner.addTest(test_normalform::suite());
  runner.addTest(test_node_conversion::suite());
  runner.addTest(test_simplify::suite());
  runner.addTest(test_expression_comparison::suite());
  runner.addTest(test_compare_utilities::suite());
  //runner.addTest(test_biomodels217_expression::suite());
  runner.run();
  return 0;
}
