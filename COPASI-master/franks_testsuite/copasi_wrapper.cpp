// Begin CVS Header
//   $Source: /Volumes/Home/Users/shoops/cvs/copasi_dev/franks_testsuite/copasi_wrapper.cpp,v $
//   $Revision: 1.11 $
//   $Name:  $
//   $Author: bergmann $
//   $Date: 2012/04/10 09:51:05 $
// End CVS Header

// Copyright (C) 2012 - 2010 by Pedro Mendes, Virginia Tech Intellectual
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

#include <iostream>
#include <stdlib.h>

#include "copasi/copasi.h"
#include "copasi/CopasiDataModel/CCopasiDataModel.h"
#include "copasi/report/CCopasiRootContainer.h"
#include "copasi/model/CMetab.h"
#include "copasi/report/CCopasiObjectName.h"
#include "copasi/utilities/CCopasiVector.h"
#include "copasi/model/CModel.h"
#include "copasi/utilities/CCopasiException.h"
#include "copasi/commandline/COptionParser.h"
#include "copasi/commandline/COptions.h"

#include "copasi/trajectory/CTrajectoryTask.h"
#include "copasi/trajectory/CTrajectoryMethod.h"
#include "copasi/trajectory/CTrajectoryProblem.h"
#include "copasi/report/CReportDefinitionVector.h"
#include "copasi/report/CReportDefinition.h"

int main(int argc, char *argv[])
{
  // Parse the commandline options
  // first argument is the SBML filename
  // second argument is the endtime
  // third argument is the step number
  // fourth argument is the filename where the results are to be written
  // fifth argument is the tmp directory (this is not needed)
  // the rest of the arguments are species names for the result
  try
    {
      // Create the root container.
      CCopasiRootContainer::init(0, NULL, false);
    }

  catch (copasi::autoexcept &e)
    {}

  catch (copasi::option_error &e)
    {}

  if (argc < 5)
    {
      std::cout << "Usage: batch_wrapper SBMLFILENAME STARTTIME ENDTIME STEPNUMBER OUTFILENAME" << std::endl;
      exit(1);
    }

  char* pSBMLFilename = argv[1];
  const char * pStartTime = argv[2];
  const char * pEndTime = argv[3];
  const char * pStepNumber = argv[4];
  char* pOutputFilename = argv[5];
  CTrajectoryTask* pTrajectoryTask = NULL;

  std::string CWD = COptions::getPWD();
  double startTime = strToDouble(pStartTime, NULL);
  double endTime = strToDouble(pEndTime, NULL);
  double stepNumber = strToDouble(pStepNumber, NULL);

  if (startTime < 0.0)
    {
      std::cerr << "Invalid endtime " << pEndTime << std::endl;
      exit(1);
    }

  if (endTime <= 0.0)
    {
      std::cerr << "Invalid endtime " << pEndTime << std::endl;
      exit(1);
    }

  if (stepNumber <= 0.0)
    {
      std::cerr << "Invalid step number " << pStepNumber << std::endl;
      exit(1);
    }

  try
    {
      // Create the global data model.
      CCopasiDataModel* pDataModel = CCopasiRootContainer::addDatamodel();

      // Import the SBML File
      pDataModel->importSBML(pSBMLFilename);

      // create a report with the correct filename and all the species against
      // time.
      CReportDefinitionVector* pReports = pDataModel->getReportDefinitionList();
      CReportDefinition* pReport = pReports->createReportDefinition("Report", "Output for batch run");
      pReport->setTaskType(CCopasiTask::timeCourse);
      pReport->setIsTable(false);
      pReport->setSeparator(CCopasiReportSeparator(", "));

      std::vector<CRegisteredObjectName>* pHeader = pReport->getHeaderAddr();
      std::vector<CRegisteredObjectName>* pBody = pReport->getBodyAddr();
      pBody->push_back(CCopasiObjectName(pDataModel->getModel()->getCN() + ",Reference=Time"));
      pBody->push_back(CRegisteredObjectName(pReport->getSeparator().getCN()));
      pHeader->push_back(CCopasiStaticString("time").getCN());
      pHeader->push_back(pReport->getSeparator().getCN());
      const CCopasiVectorNS<CCompartment>& compartments = pDataModel->getModel()->getCompartments();
      unsigned int j, jMax = compartments.size();
      /*
      for (j = 0; j < jMax;++j)
      {
        if(compartments[j]->getStatus()!=CModelEntity::FIXED)
        {
          pBody->push_back(compartments[j]->getObject(CCopasiObjectName("Reference=Volume"))->getCN());
          pBody->push_back(pReport->getSeparator().getCN());
          pHeader->push_back(CCopasiStaticString(compartments[j]->getSBMLId()).getCN());
          pHeader->push_back(pReport->getSeparator().getCN());
        }
      }
      */
      const CCopasiVector<CMetab>& metabolites = pDataModel->getModel()->getMetabolites();
      jMax = metabolites.size();

      for (j = 0; j < jMax; ++j)
        {
          if (metabolites[j]->getStatus() != CModelEntity::FIXED)
            {
              pBody->push_back(metabolites[j]->getObject(CCopasiObjectName("Reference=Concentration"))->getCN());
              pBody->push_back(pReport->getSeparator().getCN());
              pHeader->push_back(CCopasiStaticString(metabolites[j]->getSBMLId()).getCN());
              pHeader->push_back(pReport->getSeparator().getCN());
            }
        }

      /*
      const CCopasiVectorN<CModelValue>& parameters = pDataModel->getModel()->getModelValues();
      jMax = parameters.size();
      for (j = 0; j < jMax;++j)
      {
        if(parameters[j]->getStatus()!=CModelEntity::FIXED)
        {
          pBody->push_back(parameters[j]->getObject(CCopasiObjectName("Reference=Value"))->getCN());
          pBody->push_back(pReport->getSeparator().getCN());
          pHeader->push_back(CCopasiStaticString(parameters[j]->getSBMLId()).getCN());
          pHeader->push_back(pReport->getSeparator().getCN());
        }
      }
      const CCopasiVectorNS<CReaction>& reactions = pDataModel->getModel()->getReactions();
      jMax = reactions.size();
      for (j = 0; j < jMax;++j)
      {
        pBody->push_back(reactions[j]->getObject(CCopasiObjectName("Reference=Flux"))->getCN());
        pBody->push_back(pReport->getSeparator().getCN());
        pHeader->push_back(CCopasiStaticString(reactions[j]->getSBMLId()).getCN());
        pHeader->push_back(pReport->getSeparator().getCN());
      }
      */

      // delete the last separator
      if ((*pBody->rbegin()) == pReport->getSeparator().getCN())
        {
          pBody->erase(--pBody->end());
        }

      if ((*pHeader->rbegin()) == pReport->getSeparator().getCN())
        {
          pHeader->erase(--pHeader->end());
        }

      // create a trajectory task
      pTrajectoryTask = new CTrajectoryTask();
      // use LSODAR from now on since we will have events pretty soon
      pTrajectoryTask->setMethodType(CCopasiMethod::LSODAR);
      pTrajectoryTask->getProblem()->setModel(pDataModel->getModel());

      pTrajectoryTask->setScheduled(true);

      pTrajectoryTask->getReport().setReportDefinition(pReport);
      pTrajectoryTask->getReport().setTarget(CWD + "/" + pOutputFilename);
      pTrajectoryTask->getReport().setAppend(false);

      CTrajectoryProblem* pProblem = dynamic_cast<CTrajectoryProblem*>(pTrajectoryTask->getProblem());

      pProblem->setStepNumber((const unsigned C_INT32)stepNumber);
      pDataModel->getModel()->setInitialTime((const C_FLOAT64)startTime);
      pProblem->setDuration((const C_FLOAT64)endTime - startTime);
      pProblem->setTimeSeriesRequested(true);

      CTrajectoryMethod* pMethod = dynamic_cast<CTrajectoryMethod*>(pTrajectoryTask->getMethod());

      pMethod->getParameter("Absolute Tolerance")->setValue(1.0e-12);

      CCopasiVectorN< CCopasiTask > & TaskList = * pDataModel->getTaskList();

      TaskList.remove("Time-Course");
      TaskList.add(pTrajectoryTask, true);

      // save the file for control purposes
      //std::string saveFilename = pSBMLFilename;
      //saveFilename = saveFilename.substr(0, saveFilename.length() - 4) + ".cps";
      //pDataModel->saveModel(saveFilename, NULL, true);

      // Run the trajectory task

      //pTrajectoryTask->initialize(CCopasiTask::OUTPUT_UI, NULL,NULL);
      //pTrajectoryTask->process(true);
      //pTrajectoryTask->restore();

      // create another report that will write to the directory where the input file came from
      // this can be used for debugging
      // create a trajectory task
      pTrajectoryTask->getReport().setTarget(pOutputFilename);

      pTrajectoryTask->initialize(CCopasiTask::OUTPUT_UI, pDataModel, NULL);
      pTrajectoryTask->process(true);
      pTrajectoryTask->restore();
    }
  catch (CCopasiException Exception)
    {
      std::cerr << Exception.getMessage().getText() << std::endl;
    }

  std::string Text = "";

  while (CCopasiMessage::size() > 0)
    {
      const CCopasiMessage& message = CCopasiMessage::getLastMessage();

      if (message.getType() < CCopasiMessage::RAW_FILTERED)
        {
          if (Text != "") Text += "\n";

          Text += message.getText();
        }
    }

  if (Text != "") std::cerr << Text << std::endl;

  CCopasiRootContainer::destroy();

  return 0;
}
