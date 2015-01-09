// Copyright (C) 2010 - 2014 by Pedro Mendes, Virginia Tech Intellectual
// Properties, Inc., University of Heidelberg, and The University
// of Manchester.
// All rights reserved.

// Copyright (C) 2008 - 2009 by Pedro Mendes, Virginia Tech Intellectual
// Properties, Inc., EML Research, gGmbH, University of Heidelberg,
// and The University of Manchester.
// All rights reserved.

#define COPASI_MAIN

#include <iostream>
#include <fstream>
#include <string>
#include <list>
#include <sstream>
#include <stdlib.h>

#include "copasi/copasi.h"
#include "copasi/report/CCopasiRootContainer.h"
#include "copasi/CopasiDataModel/CCopasiDataModel.h"
#include "copasi/report/CCopasiContainer.h"
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
#include "copasi/report/CCopasiStaticString.h"

void parse_items(const std::string& text, std::list<std::string>& itemList)
{
  // parse the items
  //std::cout << "parsing: " << text << std::endl;
  std::istringstream iss(text);
  std::string token;
  std::string delimiter = " \t\r\n";
  std::size_t begin, end;

  while (getline(iss, token, ','))
    {
      // strip whitespaces
      begin = token.find_first_not_of(delimiter);
      end = token.find_last_not_of(delimiter);

      if (begin != std::string::npos)
        {
          //std::cout << token.substr(begin,end-begin+1) << std::endl;
          itemList.push_back(token.substr(begin, end - begin + 1));
        }
    }
}

int main(int argc, char** argv)
{
  std::string in_dir;
  std::string out_dir;
  std::string test_name;
  unsigned int level = 0;
  unsigned int version = 0;
  double start_time = 0.0;
  double duration = 0.0;
  unsigned int steps = 0;
  std::list<std::string> variables;
  std::set<std::string> amounts;
  std::set<std::string> concentrations;
  double absolute_error = 0.0;
  double relative_error = 0.0;

  if (argc != 6)
    {
      std::cerr << "Usage: sbml-testsuite INPUT_DIRECTORY TESTNAME OUTPUTDIRECTORY SBMLLEVEL SBMLVERSION" << std::endl;
      exit(1);
    }

  try
    {
      // Create the root container.
      CCopasiRootContainer::init(argc, argv, false);
    }
  catch (copasi::autoexcept &e)
    {}

  catch (copasi::option_error &e)
    {}

  in_dir = argv[1];
  test_name = argv[2];
  out_dir = argv[3];
  level = strtol(argv[4], NULL, 10);
  version = strtol(argv[5], NULL, 10);

  if (level < 1 || level > 3)
    {
      std::cerr << "SBML Level " << level << " Version " << version << " not supported." << std::endl;
      exit(1);
    }

  if (level == 1 && (version < 1 || version > 2))
    {
      std::cerr << "SBML Level " << level << " Version " << version << " not supported." << std::endl;
      exit(1);
    }

  std::ostringstream os;
  os << in_dir << "/" << test_name << "/" << test_name << "-settings.txt";
  // open the settings file and parse it
  std::cout << "parsing file: " << os.str() << std::endl;
  std::ifstream f;
  f.open(os.str().c_str());

  if (f.is_open())
    {
      std::string line;

      while (! f.eof())
        {
          getline(f, line);
          size_t pos = line.find(':');

          if (pos != std::string::npos)
            {
              std::string tag = line.substr(0, pos);
              std::string value = line.substr(pos + 1);

              if (tag == "start")
                {
                  start_time = strToDouble(value.c_str(), NULL);
                }
              else if (tag == "duration")
                {
                  duration = strToDouble(value.c_str(), NULL);
                }
              else if (tag == "steps")
                {
                  steps = strtol(value.c_str(), NULL, 10);
                }
              else if (tag == "variables")
                {
                  // parse the variables
                  parse_items(value, variables);
                }
              else if (tag == "amount")
                {
                  // parse the variables
                  std::list<std::string> tmp;
                  parse_items(value, tmp);
                  amounts.insert(tmp.begin(), tmp.end());
                }
              else if (tag == "concentration")
                {
                  // parse the variables
                  std::list<std::string> tmp;
                  parse_items(value, tmp);
                  concentrations.insert(tmp.begin(), tmp.end());
                }
              else if (tag == "absolute")
                {
                  absolute_error = strToDouble(value.c_str(), NULL);
                }
              else if (tag == "relative")
                {
                  relative_error = strToDouble(value.c_str(), NULL);
                }
              else
                {
                  std::cerr << "Warning. Unknown tag \"" << tag << "\"." << std::endl;
                }
            }
        }

      f.close();
    }
  else
    {
      std::cerr << "Error. Unable to open file \"" << os.str() << "\"." << std::endl;
      exit(1);
    }

  if (duration <= 0.0)
    {
      std::cerr << "Error. Duration < 0.0." << std::endl;
      exit(1);
    }

  if (steps == 0)
    {
      std::cerr << "Error. Number of steps set to 0." << std::endl;
      exit(1);
    }

  if (variables.empty())
    {
      std::cerr << "Error. No variables specified." << std::endl;
      exit(1);
    }

  if (start_time < 0.0)
    {
      std::cerr << "Error. Start time < 0.0." << std::endl;
      exit(1);
    }

  if (absolute_error < 0.0)
    {
      std::cerr << "Error. Absolute error < 0.0." << std::endl;
      exit(1);
    }

  if (relative_error < 0.0)
    {
      std::cerr << "Error. Relative error < 0.0." << std::endl;
      exit(1);
    }

  /*
  std::cout << "Simulating " << in_dir << "/" << test_name << "/" << test_name << "-sbml-l" << level << "v" << version << ".xml" << std::endl;
  std::cout << "Start time: " << start_time << std::endl;
  std::cout << "Duration: " << duration << std::endl;
  std::cout << "Absolute error: " << absolute_error << std::endl;
  std::cout << "Relative error: " << relative_error << std::endl;
  std::list<std::string>::const_iterator it=variables.begin(),endit=variables.end();
  std::cout << "Outputting variables: ";
  while(it!=endit)
  {
      std::cout << *it << ", ";
      ++it;
  }
  std::cout << std::endl;
  */
  os.str("");
  os << in_dir << "/" << test_name << "/" << test_name << "-sbml-l" << level << "v" << version << ".xml";
  std::string sbml_filename = os.str();
  double end_time = start_time + duration;
  os.str("");

  os << out_dir << "/" << test_name << ".csv";
  std::string output_filename = os.str();
  CTrajectoryTask* pTrajectoryTask = NULL;

  std::string CWD = COptions::getPWD();

  if (end_time == 0.0)
    {
      std::cerr << "Error. Invalid endtime " << end_time << std::endl;
      exit(1);
    }

  try
    {
      // Create the global data model.
      CCopasiDataModel* pDataModel = CCopasiRootContainer::addDatamodel();

      // Import the SBML File
      pDataModel->importSBML(sbml_filename.c_str());
      // create a report with the correct filename and all the species against
      // time.
      CReportDefinitionVector* pReports = pDataModel->getReportDefinitionList();
      CReportDefinition* pReport = pReports->createReportDefinition("Report", "Output for SBML testsuite run");
      pReport->setSeparator(CCopasiReportSeparator(","));
      pReport->setTaskType(CCopasiTask::timeCourse);
      pReport->setIsTable(false);
      std::vector<CRegisteredObjectName>* pHeaderAddr = pReport->getHeaderAddr();
      std::vector<CRegisteredObjectName>* pBodyAddr = pReport->getBodyAddr();

      //std::vector<CRegisteredObjectName>* pTable = pReport->getTableAddr();
      pHeaderAddr->push_back(CCopasiStaticString("time").getCN());
      pBodyAddr->push_back(CCopasiObjectName(pDataModel->getModel()->getValueReference()->getCN()));
      pHeaderAddr->push_back(pReport->getSeparator().getCN());
      pBodyAddr->push_back(pReport->getSeparator().getCN());
      // create a map of all possible variables
      std::map<std::string, const CModelEntity*> variableMap;
      const CCopasiVector<CMetab>& metabolites = pDataModel->getModel()->getMetabolites();
      unsigned int i, iMax = metabolites.size();

      for (i = 0; i < iMax; ++i)
        {
          variableMap.insert(std::pair<std::string, const CModelEntity*>(metabolites[i]->getSBMLId(), metabolites[i]));
        }

      const CCopasiVector<CCompartment>& compartments = pDataModel->getModel()->getCompartments();

      iMax = compartments.size();

      for (i = 0; i < iMax; ++i)
        {
          variableMap.insert(std::pair<std::string, const CModelEntity*>(compartments[i]->getSBMLId(), compartments[i]));
        }

      const CCopasiVector<CModelValue>& modelValues = pDataModel->getModel()->getModelValues();

      iMax = modelValues.size();

      for (i = 0; i < iMax; ++i)
        {
          variableMap.insert(std::pair<std::string, const CModelEntity*>(modelValues[i]->getSBMLId(), modelValues[i]));
        }

      std::list<std::string>::const_iterator it = variables.begin(), endit = variables.end();
      std::map<std::string, const CModelEntity*>::const_iterator pos;
      unsigned int dummyCount = 1;
      std::ostringstream os;

      while (it != endit)
        {
          pos = variableMap.find(*it);

          if (pos == variableMap.end())
            {
              std::cerr << "Could not find a model entity for the SBML id " << *it << std::endl;
              exit(1);
            }

          os.str("");
          os << pos->second->getSBMLId();
          pHeaderAddr->push_back(CCopasiStaticString(os.str()).getCN());

          if (dynamic_cast<const CMetab*>(pos->second) != NULL)
            {
              if (amounts.find(*it) != amounts.end())
                {
                  // create a new global value with an assignment
                  std::stringstream ss;
                  ss << "dummy_modelvalue_" << dummyCount;
                  CModelValue* pTmpMV = pDataModel->getModel()->createModelValue(ss.str());
                  ++dummyCount;
                  ss.str("");
                  assert(pTmpMV);
                  // create an assignment that takes the concentration of the
                  // metabolite and multiplies it by the compartment
                  ss << "<" << dynamic_cast<const CMetab*>(pos->second)->getConcentrationReference()->getCN() << "> * <" << dynamic_cast<const CMetab*>(pos->second)->getCompartment()->getValueReference()->getCN() << ">";
                  pTmpMV->setStatus(CModelEntity::ASSIGNMENT);
                  bool tmpRes = pTmpMV->setExpression(ss.str());
                  assert(tmpRes == true);
                  pBodyAddr->push_back(pTmpMV->getValueReference()->getCN());
                }
              else if (concentrations.find(*it) != concentrations.end())
                {
                  pBodyAddr->push_back(dynamic_cast<const CMetab*>(pos->second)->getConcentrationReference()->getCN());
                }
              else
                {
                  std::cerr << "Species \"" << *it << "\" neither given in the amounts or the conentration section." << std::endl;
                  exit(1);
                }
            }
          else if (dynamic_cast<const CCompartment*>(pos->second) != NULL)
            {
              pBodyAddr->push_back(pos->second->getValueReference()->getCN());
            }
          else if (dynamic_cast<const CModelValue*>(pos->second) != NULL)
            {
              pBodyAddr->push_back(pos->second->getValueReference()->getCN());
            }

          ++it;

          if (it != endit)
            {
              pHeaderAddr->push_back(pReport->getSeparator().getCN());
              pBodyAddr->push_back(pReport->getSeparator().getCN());
            }
        }

      // create a trajectory task
      pTrajectoryTask =  new CTrajectoryTask();
      pTrajectoryTask->getProblem()->setModel(pDataModel->getModel());

      pTrajectoryTask->setScheduled(true);

      pTrajectoryTask->getReport().setReportDefinition(pReport);
      pTrajectoryTask->getReport().setTarget(output_filename);
      pTrajectoryTask->getReport().setAppend(false);

      CTrajectoryProblem* pProblem = dynamic_cast<CTrajectoryProblem*>(pTrajectoryTask->getProblem());

      pProblem->setStepNumber((const unsigned C_INT32)steps);
      pProblem->setDuration((const C_FLOAT64)end_time);
      pProblem->setTimeSeriesRequested(true);
      pProblem->setContinueSimultaneousEvents(true);

      CTrajectoryMethod* pMethod = dynamic_cast<CTrajectoryMethod*>(pTrajectoryTask->getMethod());

      pMethod->getParameter("Absolute Tolerance")->setValue(1.0e-12);

      CCopasiVectorN< CCopasiTask > & TaskList = * pDataModel->getTaskList();

      TaskList.remove("Time-Course");
      TaskList.add(pTrajectoryTask, true);

      // save the file for control purposes
      std::string saveFilename = sbml_filename.substr(0, sbml_filename.length() - 4) + ".cps";
      pDataModel->saveModel(saveFilename, NULL, true, false);

      // Run the trajectory task

      if (!pTrajectoryTask->initialize(CCopasiTask::OUTPUT_SE, pDataModel, NULL))
        {
          std::cerr << CCopasiMessage::getAllMessageText();
        }
      else if (!pTrajectoryTask->process(true))
        {
          std::cerr << CCopasiMessage::getAllMessageText();
        }
      else if (!pTrajectoryTask->restore())
        {
          std::cerr << CCopasiMessage::getAllMessageText();
        }

      // create another report that will write to the directory where the input file came from
      // this can be used for debugging
      // create a trajectory task
      //pTrajectoryTask->getReport().setTarget(pOutputFilename);

      //pTrajectoryTask->initialize(CCopasiTask::OUTPUT_SE, pDataModel, NULL);
      //pTrajectoryTask->process(true);
      //pTrajectoryTask->restore();
    }

  catch (CCopasiException & Exception)
    {
      std::cout << Exception.getMessage().getText() << std::endl;
    }

  CCopasiRootContainer::destroy();
  return 0;
}
