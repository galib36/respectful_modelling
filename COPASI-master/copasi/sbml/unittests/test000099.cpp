// Begin CVS Header
//   $Source: /Volumes/Home/Users/shoops/cvs/copasi_dev/copasi/sbml/unittests/test000099.cpp,v $
//   $Revision: 1.2 $
//   $Name:  $
//   $Author: shoops $
//   $Date: 2011/08/02 20:44:07 $
// End CVS Header

// Copyright (C) 2011 by Pedro Mendes, Virginia Tech Intellectual
// Properties, Inc., University of Heidelberg, and The University
// of Manchester.
// All rights reserved.

#include "test000099.h"

#include <vector>
#include "utilities.hpp"

#include "copasi/report/CCopasiContainer.h"
#include "copasi/report/CCopasiObjectName.h"
#include "copasi/report/CCopasiRootContainer.h"
#include "copasi/CopasiDataModel/CCopasiDataModel.h"
#include "copasi/model/CModel.h"

void test000099::setUp()
{
  // Create the root container.
  CCopasiRootContainer::init(0, NULL, false);
  pDataModel = CCopasiRootContainer::addDatamodel();
}

void test000099::tearDown()
{
  CCopasiRootContainer::destroy();
}


// tests whether we are importing notes on all elemnts
void test000099::test_bug1675()
{
  CPPUNIT_ASSERT(pDataModel != NULL);
  CModel* pModel = NULL;
  std::vector< CCopasiContainer * > listOfContainer;

  // try to import the file which should lead to an exception
  try
    {
      pDataModel->importSBMLFromString(SBML_MODEL_BAD);
      // importing should lead to an exception, so we should never
      // reach this point
      assert(false);
    }
  catch (...)
    {
      CPPUNIT_ASSERT(pModel == NULL);
      const CModel * pModel2 =
        dynamic_cast< const CModel * >(pDataModel->ObjectFromName(listOfContainer, CCopasiObjectName("CN=Root,Model=my test model")));
      CPPUNIT_ASSERT(pModel2 == NULL);
    }

  // import the good model
  try
    {
      CPPUNIT_ASSERT(pDataModel->importSBMLFromString(SBML_MODEL_GOOD) == true);
    }
  catch (...)
    {
      // this should not lead to an exception
      assert(false);
    }

  // now check if we find the object for the model name
  pModel = pDataModel->getModel();
  CPPUNIT_ASSERT(pModel != NULL);
  const CModel * pModel2 =
    dynamic_cast< const CModel * >(pDataModel->ObjectFromName(listOfContainer, pDataModel->getModel()->getCN()));
  CPPUNIT_ASSERT(pModel == pModel2);
}

const char* test000099::SBML_MODEL_BAD =
  "<?xml version=\"1.0\" encoding=\"UTF-8\"?>\n"
  "<sbml xmlns=\"http://www.sbml.org/sbml/level2\" level=\"2\" version=\"1\">\n"
  "  <model name=\"my test model\">\n"
  "    <listOfCompartments>\n"
  "      <compartment id=\"compartment_1\" size=\"1.0\"/>\n"
  "    </listOfCompartments>\n"
  "    <listOfSpecies>\n"
  "      <species id=\"species_1\" compartment=\"compartment_1\" initialConcentration=\"1.0\"/>\n"
  "      <species id=\"species_2\" compartment=\"compartment_1\" initialConcentration=\"1.0\"/>\n"
  "    </listOfSpecies>\n"
  "    <listOfReactions>\n"
  "      <reaction id=\"reaction_1\">\n"
  "        <listOfReactants>\n"
  "          <speciesReference species=\"species_1\"/>\n"
  "        </listOfReactants>\n"
  "        <kineticLaw>\n"
  "          <math xmlns=\"http://www.w3.org/1998/Math/MathML\">\n"
  "            <apply>\n"
  "              <times/>\n"
  "              <ci> species_1 </ci>\n"
  "              <ci> species_2 </ci>\n"
  "            </apply>\n"
  "          </math>\n"
  "        </kineticLaw>\n"
  "      </reaction>\n"
  "    </listOfReactions>\n"
  "  </model>\n"
  "</sbml>\n"
  ;


const char* test000099::SBML_MODEL_GOOD =
  "<?xml version=\"1.0\" encoding=\"UTF-8\"?>\n"
  "<sbml xmlns=\"http://www.sbml.org/sbml/level2\" level=\"2\" version=\"1\">\n"
  "  <model name=\"my test model\">\n"
  "    <listOfCompartments>\n"
  "      <compartment id=\"compartment_1\" size=\"1.0\"/>\n"
  "    </listOfCompartments>\n"
  "    <listOfSpecies>\n"
  "      <species id=\"species_1\" compartment=\"compartment_1\" initialConcentration=\"1.0\"/>\n"
  "      <species id=\"species_2\" compartment=\"compartment_1\" initialConcentration=\"1.0\"/>\n"
  "    </listOfSpecies>\n"
  "    <listOfReactions>\n"
  "      <reaction id=\"reaction_1\">\n"
  "        <listOfReactants>\n"
  "          <speciesReference species=\"species_1\"/>\n"
  "        </listOfReactants>\n"
  "        <listOfModifiers>\n"
  "          <modifierSpeciesReference species=\"species_2\"/>\n"
  "        </listOfModifiers>\n"
  "        <kineticLaw>\n"
  "          <math xmlns=\"http://www.w3.org/1998/Math/MathML\">\n"
  "            <apply>\n"
  "              <times/>\n"
  "              <ci> species_1 </ci>\n"
  "              <ci> species_2 </ci>\n"
  "            </apply>\n"
  "          </math>\n"
  "        </kineticLaw>\n"
  "      </reaction>\n"
  "    </listOfReactions>\n"
  "  </model>\n"
  "</sbml>\n"
  ;



