// Begin CVS Header
//   $Source: /Volumes/Home/Users/shoops/cvs/copasi_dev/copasi/sbml/unittests/test000068.h,v $
//   $Revision: 1.2 $
//   $Name:  $
//   $Author: gauges $
//   $Date: 2009/02/18 20:41:02 $
// End CVS Header

// Copyright (C) 2008 by Pedro Mendes, Virginia Tech Intellectual
// Properties, Inc., EML Research, gGmbH, University of Heidelberg,
// and The University of Manchester.
// All rights reserved.

#ifndef TEST_000068_H__
#define TEST_000068_H__

#include <cppunit/TestFixture.h>
#include <cppunit/TestSuite.h>
#include <cppunit/TestResult.h>
#include <cppunit/extensions/HelperMacros.h>

class CCopasiDataModel;

// Bug 1068

class test000068 : public CppUnit::TestFixture
  {
    CPPUNIT_TEST_SUITE(test000068);
    CPPUNIT_TEST(test_bug1068);
    CPPUNIT_TEST_SUITE_END();

  protected:
    static const char* MODEL_STRING1;
    static CCopasiDataModel* pCOPASIDATAMODEL;

  public:
    void setUp();

    void tearDown();

    void test_bug1068();
  };

#endif /* TEST000068_H__ */
