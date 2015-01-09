// Begin CVS Header
//   $Source: /Volumes/Home/Users/shoops/cvs/copasi_dev/copasi/sbml/unittests/test000056.cpp,v $
//   $Revision: 1.6 $
//   $Name:  $
//   $Author: gauges $
//   $Date: 2010/03/11 11:52:00 $
// End CVS Header

// Copyright (C) 2010 by Pedro Mendes, Virginia Tech Intellectual
// Properties, Inc., University of Heidelberg, and The University
// of Manchester.
// All rights reserved.

// Copyright (C) 2008 by Pedro Mendes, Virginia Tech Intellectual
// Properties, Inc., EML Research, gGmbH, University of Heidelberg,
// and The University of Manchester.
// All rights reserved.

#include "test000056.h"

#include <sstream>
#include "utilities.hpp"
#include "copasi/CopasiDataModel/CCopasiDataModel.h"

#include "sbml/SBMLDocument.h"
#include "sbml/Model.h"
#include "sbml/Rule.h"
#include "sbml/Species.h"
#include "sbml/Parameter.h"
#include "sbml/math/ASTNode.h"

#include "copasi/report/CCopasiRootContainer.h"

CCopasiDataModel* test000056::pCOPASIDATAMODEL = NULL;

void test000056::setUp()
{
  // Create the root container.
  CCopasiRootContainer::init(0, NULL, false);
  // Create the global data model.
  pCOPASIDATAMODEL = CCopasiRootContainer::addDatamodel();
}

void test000056::tearDown()
{
  CCopasiRootContainer::destroy();
}

void test000056::test_bug1005()
{
  // make sure models without compartments don't lead to a crash
  CCopasiDataModel* pDataModel = pCOPASIDATAMODEL;

  try
    {
      CPPUNIT_ASSERT(pDataModel->importSBMLFromString(test000056::MODEL_STRING));
    }
  catch (...)
    {
      CPPUNIT_ASSERT(false);
    }

  const SBMLDocument* pDocument = pDataModel->getCurrentSBMLDocument();

  const Model* pSBMLModel = pDocument->getModel();

  CPPUNIT_ASSERT(pSBMLModel != NULL);

  CPPUNIT_ASSERT(pSBMLModel->getNumCompartments() == 0);

  CPPUNIT_ASSERT(pSBMLModel->getNumSpecies() == 0);

  CPPUNIT_ASSERT(pSBMLModel->getNumReactions() == 0);

  CPPUNIT_ASSERT(pSBMLModel->getNumInitialAssignments() == 0);

  CPPUNIT_ASSERT(pSBMLModel->getNumParameters() == 1);

  const Parameter* pParameter1 = pSBMLModel->getParameter(0);

  CPPUNIT_ASSERT(pParameter1 != NULL);

  CPPUNIT_ASSERT(pParameter1->getConstant() == true);

  CPPUNIT_ASSERT(pParameter1->getId() == "parameter_1");

  CPPUNIT_ASSERT(pParameter1->getName() == "A");

  CPPUNIT_ASSERT(pParameter1->isSetValue() == true);

  CPPUNIT_ASSERT(fabs((pParameter1->getValue() - 12.0) / 12.0) < 1e-15);
}

const char* test000056::MODEL_STRING =
  "<?xml version=\"1.0\" encoding=\"UTF-8\"?>\n"
  "<sbml xmlns=\"http://www.sbml.org/sbml/level2\" level=\"2\" version=\"1\">\n"
  "  <model id=\"Model_1\" name=\"New Model\">\n"
  "    <listOfUnitDefinitions>\n"
  "      <unitDefinition id=\"volume\">\n"
  "        <listOfUnits>\n"
  "          <unit kind=\"litre\" scale=\"-3\"/>\n"
  "        </listOfUnits>\n"
  "      </unitDefinition>\n"
  "      <unitDefinition id=\"substance\">\n"
  "        <listOfUnits>\n"
  "          <unit kind=\"mole\" scale=\"-3\"/>\n"
  "        </listOfUnits>\n"
  "      </unitDefinition>\n"
  "    </listOfUnitDefinitions>\n"
  "    <listOfParameters>\n"
  "      <parameter id=\"parameter_1\" name=\"A\" value=\"12\" constant=\"true\"/>\n"
  "    </listOfParameters>\n"
  "  </model>\n"
  "</sbml>\n"
  ;
