#include <iostream>
#include <fstream>
#include <string>
#include <cmath>
#include <sstream>
#include <iomanip>
#include <stdexcept>
#include <fstream>
#include <random>
#include <chrono>


#define COPASI_MAIN

#include <copasi/copasi.h>
#include <copasi/report/CCopasiRootContainer.h>
#include <copasi/CopasiDataModel/CCopasiDataModel.h>
#include <copasi/model/CModel.h>
#include <copasi/model/CModelValue.h>
#include <copasi/utilities/CAnnotatedMatrix.h>
#include <copasi/steadystate/CSteadyStateTask.h>
#include <copasi/steadystate/CMCATask.h>
#include <copasi/steadystate/CMCAMethod.h>

#include "copasi/report/CCopasiRootContainer.h"

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
#include "copasi/utilities/CCopasiParameter.h"


#include <sbml/SBMLTypes.h>





using namespace std;

int main(int argc, char *argv[])
{
   if(argc < 1)
   {
       cout << "No given SBML file." <<endl;
       exit(1);
   }
   else if(argc < 2)
   {
       cout << " Your command line should have the following format : " <<endl;
       cout << " ./ResearchP model_name parameters_file_name" << endl;
       exit(0);
   }
   string sourceSMBLFileName = argv[1];

    CCopasiRootContainer::init(0, NULL);
    CCopasiDataModel* pDataModel = CCopasiRootContainer::addDatamodel();
    assert(CCopasiRootContainer::getDatamodelList()->size() == 1);
    bool successImport = pDataModel->importSBML(sourceSMBLFileName, NULL);
    if(!successImport)
    {
        cout << "Cannot import the SBML file " << sourceSMBLFileName << endl;
        exit(1);
    }
    CModel* cModel = pDataModel->getModel();
    int nbrReactions = cModel->getReactions().size();
    int nbrSpecies = cModel->getMetabolites().size();
    int nbrParam;
    CReaction* currentReaction;
    string currentParamName;
    CCopasiParameter::Value currentParamValue;
    string outFilename = argv[2];
    ofstream outputParam(outFilename);

    for(int j=0; j<nbrReactions; j++)
        {
            currentReaction = cModel->getReactions()[j];
            nbrParam = currentReaction->getParameters().size();
            for(int k=0; k<nbrParam; k++)
            {
                currentParamName = currentReaction->getSBMLId()+"_"+currentReaction->getParameters().getName(k);
                currentParamValue = currentReaction->getParameters().getValue(currentReaction->getParameters().getName(k));
                outputParam<<currentParamName<<"\t"<<*currentParamValue.pDOUBLE<<endl;

            }
        }

}
