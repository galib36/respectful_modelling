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

#include "copasi/steadystate/CSteadyStateTask.h"
#include "copasi/steadystate/CSteadyStateMethod.h"
#include "copasi/steadystate/CSteadyStateProblem.h"

#include "copasi/report/CReportDefinitionVector.h"
#include "copasi/report/CReportDefinition.h"

#include <sbml/SBMLTypes.h>

#include "exprtk.hpp"

#define _USE_MATH_DEFINES

using namespace std;

/**  STATISTICS FUNCTIONS  **/

// Code from http://www.johndcook.com/normal_cdf_inverse.html
double RationalApproximation(double t)
{
    // Abramowitz and Stegun formula 26.2.23.
    // The absolute value of the error should be less than 4.5 e-4.
    double c[] = {2.515517, 0.802853, 0.010328};
    double d[] = {1.432788, 0.189269, 0.001308};
    return t - ((c[2]*t + c[1])*t + c[0]) /
               (((d[2]*t + d[1])*t + d[0])*t + 1.0);
}

double NormalCDFInverse(double p)
{
    if (p <= 0.0 || p >= 1.0)
    {
        std::stringstream os;
        os << "Invalid input argument (" << p
           << "); must be larger than 0 but less than 1.";
        throw std::invalid_argument( os.str() );
    }

    // See article above for explanation of this section.
    if (p < 0.5)
    {
        // F^-1(p) = - G^-1(p)
        return -RationalApproximation( sqrt(-2.0*log(p)) );
    }
    else
    {
        // F^-1(p) = G^-1(1-p)
        return RationalApproximation( sqrt(-2.0*log(1-p)) );
    }
}
// End of code from http://www.johndcook.com/normal_cdf_inverse.html

//Code from Fiona originally written in Python
pair<double,double> getNormalParamFromLogNormParam_median(double median, double s)
{
	double mu = log(median);
	double x = exp(2*mu);
	double sigma_2 = log( (x + sqrt(x*(4*pow(s,2)+x))) / (2*x) );
	double sigma = sqrt(sigma_2);
	return make_pair(mu, sigma);
}

pair<double,double> getNormalParamFromLogNormParam_maxci_median( double median, double maxci)
{
    double z = NormalCDFInverse(1 - 0.05/2);
	double mu = log(median);
	double sigma=(log(maxci)-mu)/z;
	return make_pair(mu, sigma);
}

pair<double,double> getNormalParamFromLogNormParam_minci_median( double median, double minci)
{
    double z = - NormalCDFInverse(1 - 0.05/2);
	double mu = log(median);
	double sigma=(log(minci)-mu)/z;
	return make_pair(mu, sigma);
}
//End of code from Fiona

/**  READING PARAMETERS FILE FUNCTIONS  **/
struct parameterDistribution
{
    string distribution_type;
    double param1; // mu for log-normal and normal distributions, min for uniform distributions
    double param2; // sigma for log-normal and normal distributions, max for uniform distributions
};

/**Analysis of the parameters distributions file
    Layout
        delimitation between information = \t = tabulation
        For each kind of distribution type (possibles are lognormal defined by mean and standard deviation, lognormal defined by mean and maxci, normal and uniform) :
        1rst row  = name of the distribution (log-normal distribution median, log-normal distribution maxci_median, normal distribution or uniform distribution)
        2dn row = parameters' name : reactionName_parameterName with reationName and parameterName the same as in the SBML file
        3rd row = first parameter value (for lognormal and normal = mean, for uniform = minimum)
        4th row = second parameter value (for lognormal and normal = standard deviation of maxci, for uniform = maximum)
**/

vector<pair<string, struct parameterDistribution>> readParameterFile(string inputFileName)
{
	ifstream file(inputFileName, ios::in);
	string temp;
	unsigned posTab;
	getline(file,temp);
	string paramName;
	double param1, param2;
	string distributionType;
	double minCI, maxCI;
	string equation = "";
	struct parameterDistribution paramDistrib;
	vector<pair<string, struct parameterDistribution>> parametersMap;
	pair<double, double> muSigmaTemp;

	while(getline(file,temp))
	{
	    posTab = temp.find('\t');
	    paramName = temp.substr(0,posTab);
	    temp = temp.substr(posTab+1);

	    posTab = temp.find('\t');
	    param1 = atof(temp.substr(0,posTab).c_str());
	    temp = temp.substr(posTab+1);

	    posTab = temp.find('\t');
	    param2 = atof(temp.substr(0,posTab).c_str());
	    temp = temp.substr(posTab+1);

	    posTab = temp.find('\t');
	    distributionType = temp.substr(0,posTab);
	    temp = temp.substr(posTab+1);

	    posTab = temp.find('\t');
	    minCI = atof(temp.substr(0,posTab).c_str());
	    temp = temp.substr(posTab+1);

	    posTab = temp.find('\t');
	    maxCI = atof(temp.substr(0,posTab).c_str());
	    temp = temp.substr(posTab+1);

	    posTab = temp.find('\t');
	    equation = temp.substr(0,posTab);

		if(minCI != 0)
		{
			//Calculation of mu and sigma
            muSigmaTemp = getNormalParamFromLogNormParam_minci_median(param1, minCI);
            //Saving in struct
            paramDistrib.distribution_type = "lognormal";
            paramDistrib.param1 = muSigmaTemp.first;
            paramDistrib.param2 = muSigmaTemp.second;

		}
		else if (maxCI !=0)
		{
			//Calculation of mu and sigma
            muSigmaTemp = getNormalParamFromLogNormParam_maxci_median(param1, maxCI);
            //Saving in struct
            paramDistrib.distribution_type = "lognormal";
            paramDistrib.param1 = muSigmaTemp.first;
            paramDistrib.param2 = muSigmaTemp.second;
		}
		else if(distributionType == "calculated")
		{
			paramDistrib.distribution_type = equation;
			paramDistrib.param1 = 0;
            paramDistrib.param2 = 0;
		}
		else
		{
			if(distributionType == "lognormal")
			{
                //Calculation of mu and sigma
                muSigmaTemp = getNormalParamFromLogNormParam_median(param1, param2);
                //Saving in struct
                paramDistrib.distribution_type = "lognormal";
                paramDistrib.param1 = muSigmaTemp.first;
                paramDistrib.param2 = muSigmaTemp.second;
			}
			else
			{
                paramDistrib.distribution_type = distributionType;
                paramDistrib.param1 = param1;
                paramDistrib.param2 = param2;
			}
		}
		parametersMap.push_back(make_pair(paramName,paramDistrib));
	}
    return parametersMap;
}

double calcExpr(map <string, double> paramMap, string expr)
{
    map <string, double>::iterator it;

    typedef exprtk::symbol_table<double>      symbol_table_t;
    typedef exprtk::expression<double>          expression_t;
    typedef exprtk::parser<double>                  parser_t;
    typedef exprtk::parser_error::type               error_t;

    symbol_table_t symbol_table;
    expression_t expression;
    expression.register_symbol_table(symbol_table);

    for(it=paramMap.begin(); it!=paramMap.end(); it++)
    {
        symbol_table.add_variable(it->first, it->second);
    }
    symbol_table.add_constants();

    parser_t parser;
    parser.compile(expr,expression);

    return expression.value();
}

int main(int argc, char *argv[])
{

/*To run the software
./ResearchP model_name parameters_file_name major_reaction_for_MCA type_of_elasticities_for_MCA number_of_iterations TimeSeries_Duration TimeSeries_StepSize TimeSeries_Number "*/
   if(argc == 1)
   {
       cout << "No given SBML file." <<endl;
       exit(1);
   }
   else if(argc == 2)
   {
       cout << " Your command line should have the following format : " <<endl;
       cout << " If you want the sotware to ample the parameters : " << endl;
       cout << " ./ResearchP model_name parameters_file_name major_reaction_for_MCA type_of_elasticities_for_MCA number_of_iterations" << endl;
       cout << " If your parameters file contains sampled parameters : " << endl;
       cout << " ./ResearchP model_name parameters_file_name major_reaction_for_MCA type_of_elasticities_for_MCA" << endl;
       exit(0);
   }
   else if(argc < 3)
   {
       cout << "Wrong parameters" <<endl;
       cout << " Your command line should have the following format : " <<endl;
       cout << " If you want the sotware to ample the parameters : " << endl;
       cout << " ./ResearchP model_name parameters_file_name major_reaction_for_MCA type_of_elasticities_for_MCA number_of_iterations" << endl;
       cout << " If your parameters file contains sampled parameters : " << endl;
       cout << " ./ResearchP model_name parameters_file_name major_reaction_for_MCA type_of_elasticities_for_MCA" << endl;
       exit(1);
   }

    //Source SBML file : saving the name given as input
    string sourceSMBLFileName = argv[1];

    //Source SBML file : saving the name given as input
    string parameterDistributionFileName = argv[2];

    //Saving the major reaction and type of elasticities (scaled or unscaled, unscaled by default) to focus Metabolic Control Analysis results
    string majorReactionMCAName = argv[3];
    string elasticity = argv[4];
    ifstream inputParamFile(parameterDistributionFileName);
    string line;
    if(elasticity !="unscaled" && elasticity!="scaled")
    {
        cout << "Type of elasticity given not compatible, unscaled chosen by default." << endl;
        elasticity= "scaled";
    }

    //Saving the number of iterations given as input
    bool sampledParam = true;
    int nbrIterations;
    if(argc > 5)
    {
        nbrIterations = atoi(argv[5]);
        if(nbrIterations > 0)
        {
            sampledParam = false;

        }
    }

    cout << "Copasi steady-state and metabolic control analysis using sampled kinetic parameters are running." << endl;
    //ITERATIONS

    //Necessary if parameters need to be sampled
    map<string, lognormal_distribution<double> > parametersMapLogNormal;
    map<string, lognormal_distribution<double> >::iterator itLogNormal;
    map<string, normal_distribution<double> > parametersMapNormal;
    map<string, normal_distribution<double> >::iterator itNormal;
    map<string, uniform_real_distribution<double> > parametersMapUniform;
    map<string, uniform_real_distribution<double> >::iterator itUniform;
    vector< pair<string, string> > parametersMapCalculated;
    unsigned seed;
    std::default_random_engine generator;
    map <string, double> paramValuesMap;
    map<string, double >::iterator itParam;

    //Necessary if parameters already sampled
    vector<string> parametersNames;
    vector<string>::iterator itParamNames;
    int indexParam = 0;
    int vectorSize=1000000;
    int paramNo = 200;
    //vector<vector<double>> parametersValues(paramNo,vector<double>(vectorSize));
    vector<double> parametersValues[100000];


    if(!sampledParam)
    {
        cout << "Your parameters file contains values to create parameters distributions." << endl;
        vector<pair<string, struct parameterDistribution>> parametersMap = readParameterFile(parameterDistributionFileName);

        //Creation of the parameters distributions
        for(int i=0 ; i< parametersMap.size(); i++)
        {
            if(parametersMap[i].second.distribution_type == "lognormal")
				parametersMapLogNormal[parametersMap[i].first] = lognormal_distribution<double>(parametersMap[i].second.param1, parametersMap[i].second.param2);
            else if(parametersMap[i].second.distribution_type == "normal")
                parametersMapNormal[parametersMap[i].first] = normal_distribution<double>(parametersMap[i].second.param1, parametersMap[i].second.param2);
			else if(parametersMap[i].second.distribution_type == "uniform")
				parametersMapUniform[parametersMap[i].first] = uniform_real_distribution<double>(parametersMap[i].second.param1, parametersMap[i].second.param2);
			else
				parametersMapCalculated.push_back(make_pair(parametersMap[i].first ,parametersMap[i].second.distribution_type));
        }
    }
    else
    {
       cout << "Your parameters file contains sampled parameters." << endl;
       //Necessary to read the parameters file

       //Analysis of the header
        getline(inputParamFile,line);
        long unsigned posTab;
        posTab = line.find("\t");
        while(posTab != string::npos)
        {
            parametersNames.push_back(line.substr(0,posTab));
            line = line.substr(posTab+1);
            posTab = line.find("\t");
        }

        //Reading and analysis
        while(getline(inputParamFile, line))
        {
            posTab = line.find("\t");
            while(posTab != string::npos)
            {
                parametersValues[indexParam].push_back(atof(line.substr(0,posTab).c_str()));
                line = line.substr(posTab+1);
                posTab = line.find("\t");
                indexParam++;
            }
            indexParam = 0;
       }
        nbrIterations =  parametersValues[0].size();
        //nbrIterations =  10;
        cout<<"No. of Iteration: "<<nbrIterations<<endl;
    }

    cout << "Running ..." << endl;
    // Necessary to write output SBML files
    string outputFileNameSBML;
    string outputFileNameTXT;
    string outputFileNameCPS;
    stringstream sstm;
    size_t position;

    //Necessary to access to the kinetics parameters
    CCopasiRootContainer::init(0, NULL);
    CCopasiDataModel* pDataModel = CCopasiRootContainer::addDatamodel();
    assert(CCopasiRootContainer::getDatamodelList()->size() == 1);
    bool successImport = pDataModel->importSBML(sourceSMBLFileName, NULL);
    cout<<"Success Import:"<<successImport<<endl;
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
    double currentParamValue;

    //Necessary to do Copasi analysis
    CCopasiVectorN< CCopasiTask >taskList;
    //CSteadyStateTask* pTaskS;
    CSteadyStateMethod* pMethodS;
    CMCATask* pTaskMCA;
    CMCAMethod* pMethodMCA;
    CSteadyStateMethod::ReturnCode rc;

    //Time Series
    CTrajectoryTask* pTaskT;
    CTrajectoryMethod* pMethodT;
    CTimeSeries pTimeData;

    //Writing in output files
    ifstream fileRead;
    ofstream outputSteady;
    ofstream outputSSConcentration("resultSteadyStateSpeciesConcentration.txt");
    ofstream outputSSConcentrationRate("resultSteadyStateSpeciesConcentrationRate.txt");
    ofstream outputSSTransitionTime("resultSteadyStateTransitionTime.txt");
    ofstream outputSSReaction("resultSteadyStateReactionFlux.txt");
    ofstream outputTSConcentration("resultTimeSeriesConcentration.txt");
    ofstream outputTSFlux("resultTimeSeriesFlux.txt");
    ofstream outputMCA("resultMCAFlux.txt");
    ofstream outputMCACon("resultMCAConc.txt");
    ofstream outputMCAEla("resultMCAEla.txt");
    ofstream outputParam("resultParameters.txt");
    ofstream outputTimeSeries;
    ofstream logofs("Output_log.log");
    ofstream resultStatus("resultStatus.txt");
    ofstream outputSSMConcentration("resultSSConcentration.txt");
    ofstream outputSSMReaction("resultSSFlux.txt");
    //ofstream ofs("output_log.log");
    int indexMCAReaction;
    int foundSS=0;
    int stepDuration=atoi(argv[6]);
    int stepSize=atoi(argv[7]);
    int stepNumber=atoi(argv[8]);
    CMatrix<C_FLOAT64> mcaResults;
    CMatrix<C_FLOAT64> mcaConResults;
    CMatrix<C_FLOAT64> mcaElaResults;
    string listParamNames = "";
    string reactionNames = "";
    string speciesNames = "";
    string speciesConcentrations = "";
    string speciesConcentrationsRates = "";
    string speciesSSTransitionTime = "";
    std::ostringstream os;
    string test;
    string suc = "Success";
    //string pOutputLogFilename = "Output_log";
    string timePoint;
    std::size_t SSfound;
    //float timeSeriesData[50][1000][10000];


    for (int i = 0; i<nbrIterations; i++)
    {
        //new files for each iterations : SBML file with new parameters values, CPS (Copasi) file and text file = Copasi's report
        sstm.clear();
        sstm.str(string());
        position = sourceSMBLFileName.find('.');
        sstm << sourceSMBLFileName.substr(0,position)+"_v" << i;
        outputFileNameSBML = sstm.str()+".xml";
        outputFileNameTXT = sstm.str()+".txt";
        outputFileNameCPS = sstm.str()+".cps";

        if(!sampledParam)
        {
            //Calculation of the parameters values
			//Those wih lognormal distribution
			for(itLogNormal=parametersMapLogNormal.begin() ; itLogNormal!=parametersMapLogNormal.end() ; ++itLogNormal)
			{
				seed = std::chrono::system_clock::now().time_since_epoch().count();
                generator = std::default_random_engine(seed);
				paramValuesMap[itLogNormal->first] = itLogNormal->second(generator);
			}
			//Those with normal distribution
			for(itNormal=parametersMapNormal.begin() ; itNormal!=parametersMapNormal.end() ; ++itNormal)
			{
				seed = std::chrono::system_clock::now().time_since_epoch().count();
                generator = std::default_random_engine(seed);
				paramValuesMap[itNormal->first] = itNormal->second(generator);
			}
			//Those with uniform distribution
			for(itUniform=parametersMapUniform.begin() ; itUniform!=parametersMapUniform.end() ; ++itUniform)
			{
				seed = std::chrono::system_clock::now().time_since_epoch().count();
                generator = std::default_random_engine(seed);
				paramValuesMap[itUniform->first] = itUniform->second(generator);
			}
			//Those which need to be calculated
			for(int indexCalc=0; indexCalc < parametersMapCalculated.size(); indexCalc++)
			{
				paramValuesMap[parametersMapCalculated[indexCalc].first] = calcExpr(paramValuesMap, parametersMapCalculated[indexCalc].second);
			}
        }
        string vectorParamName;
        indexMCAReaction = -1;
        //Modification of the model
        for(int j=0; j<nbrReactions; j++)
        {
            currentReaction = cModel->getReactions()[j];
            nbrParam = currentReaction->getParameters().size();

            if(currentReaction->getObjectName() == majorReactionMCAName)
            {
                indexMCAReaction=j;
                cout<<currentReaction->getObjectName()<<" \t"<<majorReactionMCAName<<"\t"<<indexMCAReaction<<endl;
            }
            for(int k=0; k<nbrParam; k++)
            {
                currentParamName = currentReaction->getSBMLId()+"_"+currentReaction->getParameters().getName(k);
                if(i==0)
                    listParamNames += currentReaction->getSBMLId()+"_"+currentReaction->getParameters().getName(k) + '\t';
                    //listParamNames += currentReaction->getParameters().getName(k) + '\t';

                if(!sampledParam)
                {
                    //Randon pick of the parameter new value
                    itParam = paramValuesMap.find(currentParamName.c_str());
                    if(itParam!=paramValuesMap.end())
                        currentParamValue = itParam->second;
                }
                else
                {
                    for(int index=0; index<parametersNames.size(); index++)
                    {
                        //vectorParamName = currentReaction->getSBMLId()+"_"+parametersNames[index];
                        vectorParamName = parametersNames[index];
                        //cout<<"Current parameter: "<<currentParamName<<", ParametersNames: "<<vectorParamName<<endl;
                        if(vectorParamName==currentParamName)
                        {
                            currentParamValue = parametersValues[index][i];
                            //cout<<"Current parameter value: "<<currentParamValue<<endl;
                        }
                        vectorParamName="";
                    }
                }

                //Modification of the parameter value in the new SBML file
                currentReaction->setParameterValue(currentReaction->getParameters().getName(k), currentParamValue, true);
                os << currentParamValue << '\t';
            }
        }

        //Result parameters file
        if(i==0)
            outputParam << listParamNames << '\n';
        outputParam << os.str() << '\n';
        os.str("");
        os.clear();
        //Saving the new version in a SBML file
        pDataModel->exportSBML(outputFileNameSBML, true, 2, 1, false, true, NULL);

        try
        {

            //Steady-StateANALYSIS
            taskList = *pDataModel->getTaskList();

            CSteadyStateTask* pTaskS = new CSteadyStateTask();
            /*
            pTaskS->setMethodType(CCopasiMethod::Newton);
            pTaskS->getProblem()->setModel(pDataModel->getModel());

            pTaskS->setScheduled(true);

            //**** specify problem ****
            //CSteadyStateProblem* pProblem = dynamic_cast<CSteadyStateProblem*>(pTaskS->getProblem());

            //**** method settings ****
            CSteadyStateMethod* pMethod = dynamic_cast<CSteadyStateMethod*>(pTaskS->getMethod());
            pMethod->getParameter("Maximum duration for forward integration")->setValue(1.0e5);
            pMethod->getParameter("Use Back Integration")->setValue(false);
            pMethod->getParameter("Maximum duration for backward integration")->setValue(1.0e3);

            //***** run the task *********
            CCopasiVectorN< CCopasiTask > & TaskList = * pDataModel->getTaskList();

            TaskList.remove("Steady-State");
            TaskList.add(pTaskS, true);

            // Run the task
            pTaskS->initialize(CCopasiTask::OUTPUT_SE, pDataModel, NULL);
            pTaskS->process(true);
            pDataModel->finish();
            pTaskS->restore();

            //status output
            rc = pTaskS->getResult();
            */




            taskList.remove("Steady-State");
            pTaskS = (CSteadyStateTask*)pDataModel->addTask(CCopasiTask::steadyState);
            pMethodS = (CSteadyStateMethod*)pTaskS->getMethod();
            pMethodS->setValue("Resolution", 0.000001);
            pTaskS->initialize(CCopasiTask::OUTPUT_UI, pDataModel, NULL);
            pTaskS->process(true);
            pDataModel->finish();
            pTaskS->restore();

            rc = pTaskS->getResult();
            //cout<<"Status Code: "<<pTaskS->getResult()<<" Found value is "<<CSteadyStateMethod::found<<endl;
            std::string iter="\n Iteration No. ";
            std::string result_iter;
            std::stringstream log_sstm;
            log_sstm << iter << i<<"==========================\n\n";
            result_iter = log_sstm.str();
            log_sstm.clear();

            std::string methodlog = dynamic_cast<CSteadyStateMethod*>(pTaskS->getMethod())->getMethodLog();
            //std::string logname = pOutputLogFilename + std::string(".log");
            //std::ofstream logofs(logname.c_str());
            logofs<<result_iter;
            logofs << methodlog;
            //cout<<methodlog;
            //logofs.close();
            //ofs << result_iter;
            //ofs << methodlog;
            SSfound=0;
            SSfound = methodlog.find(suc);
            if (rc == CSteadyStateMethod::found)
            {
                cout << "Steady-State found for Iteration: "<<i<<'\n';
                resultStatus<<1<<endl;
            }

            else
            {
                cout<<"No Steady-state was found for Iteration: "<<i<<'\n';
                resultStatus<<0<<endl;
            }


            //ofs.close();

            /*Results sort : to simplify the use of data, Copasi's report is analysis to provide 4 files each containing an array. Each file can then be used with Excel or R.The four files contain :
            * steady state species concentration
            * steady state species concentration rate
            * steady state reaction flux
        */
        }
        catch (CCopasiException Exception)
        {
            std::cerr << Exception.getMessage().getText() << std::endl;
        }

        //cerr << "Setting Resolution Error: "<<CCopasiMessage::getAllMessageText() << endl;
        //CCopasiMessage::clearDeque();

        if(i==0)
        {
            for(int j=0; j<nbrReactions; j++)
            {
                reactionNames += pDataModel->getModel()->getReactions()[j]->getObjectName()+ '\t';
            }
            outputSSReaction << reactionNames <<'\n';
            outputMCA << reactionNames << '\n';
            //outputMCAEla << reactionNames << '\n';

            for(int j=0; j<nbrSpecies; j++)
            {
                speciesNames += pDataModel->getModel()->getMetabolites()[j]->getObjectName() + '\t';
            }
            outputSSConcentration << speciesNames <<'\n';
            outputSSConcentrationRate << speciesNames <<'\n';
            outputMCAEla << speciesNames << '\n';
        }
        if(rc == CSteadyStateMethod::foundNegative)
        {
            outputSSReaction << "Negative steady-state found" <<'\n';
            outputSSConcentration << "Negative steady-state found" <<'\n';
            outputSSConcentrationRate << "Negative steady-state found" <<'\n';
        }
        else if(rc == CSteadyStateMethod::notFound)
        {
            outputSSReaction << "Steady-state not found" <<'\n';
            outputSSConcentration << "Steady-state not found" <<'\n';
            outputSSConcentrationRate << "Steady-state not found" <<'\n';
        }
        else if(rc == CSteadyStateMethod::foundEquilibrium)
        {
            outputSSReaction << "Steady-state with Equilibrium (Flux =0)" <<'\n';
            outputSSConcentration << "Steady-state with Equilibrium (Flux =0)" <<'\n';
            outputSSConcentrationRate << "Steady-state with Equilibrium (Flux =0)" <<'\n';
        }
        //if(SSfound!=std::string::npos)
        else if(rc == CSteadyStateMethod::found)
        {
            foundSS++;
            for(int j=0; j<nbrReactions; j++)
            {
                if(i==0)
                {
                    reactionNames += pDataModel->getModel()->getReactions()[j]->getObjectName()+ '\t';
                }
                os << pDataModel->getModel()->getReactions()[j]->getFlux() << '\t';
            }
            outputSSReaction << os.str() <<'\n';
            outputSSMReaction<<os.str() <<'\n';
            os.str("");
            os.clear();

            for(int j=0; j<nbrSpecies; j++)
            {
                os << pDataModel->getModel()->getMetabolites()[j]->getConcentration();
                speciesConcentrations += os.str() + '\t';
                os.str("");
                os.clear();
                os << pDataModel->getModel()->getMetabolites()[j]->getConcentrationRate();
                speciesConcentrationsRates += os.str() + '\t';
                os.str("");
                os.clear();
                os << pDataModel->getModel()->getMetabolites()[j]->getTransitionTime();
                //cout<<pDataModel->getModel()->getMetabolites()[j]->getTransitionTime();
                speciesSSTransitionTime += os.str() + '\t';
                os.str("");
                os.clear();

            }
            outputSSConcentration << speciesConcentrations <<'\n';
            outputSSMConcentration << speciesConcentrations <<'\n';
            outputSSConcentrationRate << speciesConcentrationsRates <<'\n';
            outputSSTransitionTime << speciesSSTransitionTime <<'\n';
            speciesConcentrations = "";
            speciesConcentrationsRates = "";
            speciesSSTransitionTime = "";
        }
        //else
        //{
        //    cout<<"Steady-State cannot be calculated"<<endl;
        //}
         //write method log to separate file

        //Time Series
        if(rc == CSteadyStateMethod::found)
        {
            taskList.remove("Time Course");
            pTaskT = (CTrajectoryTask*)pDataModel->addTask(CCopasiTask::timeCourse);
            pMethodT = (CTrajectoryMethod*)pTaskT->getMethod();
            CTrajectoryProblem* pProblem = dynamic_cast<CTrajectoryProblem*>(pTaskT->getProblem());
            pProblem->setDuration(stepDuration);
            pProblem->setStepNumber(stepNumber);
            pProblem->setStepSize(stepSize);
            pTaskT->initialize(CCopasiTask::OUTPUT_UI, pDataModel, NULL);
            //pMethodT->step(0.1);
            //pTimeEle->compile(CCopasiRootContainer::getDatamodelList()->, pDataModel);

            pTaskT->processStart(true);
            pTaskT->processStep(1);
            pTaskT->process(true);
            pDataModel->finish();
            pTaskT->restore();
            //cout<<"Initial Concetration: "<<speciesNames<<pDataModel->getModel()->getMetabolites()[1]->getInitialConcentration()<<endl;
            os.str("");
            os.clear();
            timePoint = "";
            //Writing the time points
            timePoint+="Time\t";
            for(int l=0;l<=stepNumber;l++)
            {
                os << pTaskT->getTimeSeries().getConcentrationData(l,0);
                timePoint += os.str() + '\t';
                os.str("");
                os.clear();
            }

            outputTSConcentration<< timePoint <<'\n';

                for(int j=0; j<nbrSpecies; j++)
                {
                    speciesConcentrations = "";

                    speciesNames="";
                    speciesNames = pDataModel->getModel()->getMetabolites()[j]->getObjectName() + '\t';
                    speciesConcentrations +=speciesNames;
                    if(i==0)
                    {
                        timePoint="";
                        timePoint = "Time\t";
                        for(int l=0;l<=stepNumber;l++)
                        {
                            os << pTaskT->getTimeSeries().getConcentrationData(l,0);
                            timePoint += os.str() + '\t';
                            os.str("");
                            os.clear();
                        }
                        remove(speciesNames.c_str());
                        outputTimeSeries.open(speciesNames.c_str(), ios::out | ios::app);
                        outputTimeSeries<<timePoint<<endl;
                        outputTimeSeries.close();
                    }

                    for(int l=0;l<=stepNumber;l++)
                    {
                        os << pTaskT->getTimeSeries().getConcentrationData(l,j+1);
                        speciesConcentrations += os.str() + '\t';
                        os.str("");
                        os.clear();
                        //timeSeriesData[i][j][l]= pTaskT->getTimeSeries().getConcentrationData(l,j+1);
                    }
                    outputTimeSeries.open(speciesNames.c_str(), ios::out | ios::app);
                    outputTimeSeries<<speciesConcentrations<<endl;
                    outputTimeSeries.close();
                    //outputTSConcentration<< speciesConcentrations <<'\n';
                }
                speciesConcentrations = "";


        }




        //MCA
        if(indexMCAReaction!=-1 && rc == CSteadyStateMethod::found)
        {
            taskList.remove("Metabolic Control Analysis");
            pTaskMCA = (CMCATask*)pDataModel->addTask(CCopasiTask::mca);
            CMCAMethod* pMethodMCA = dynamic_cast<CMCAMethod *>(pTaskMCA->getMethod());
            pTaskMCA->initialize(CCopasiTask::OUTPUT_UI, pDataModel, NULL);
            pMethodMCA->setSteadyStateResolution(1e-18);
            pTaskMCA->process(true);

            //pDataModel->finish();
            //pMethodMCA=(CMCAMethod*)pTaskMCA->getMethod();



            //pTaskMCA->restore();
            cout<<"Calculated Concentration: "<<pMethodMCA->calculateUnscaledConcentrationCC()<<endl;


            //Results : MCA coefficient of the selected reaction (chosen in input) over all the species
            if(elasticity == "unscaled")
            {
                mcaResults = pMethodMCA->getUnscaledFluxCC();
                mcaConResults = pMethodMCA->getUnscaledConcentrationCC();
                mcaElaResults = pMethodMCA->getUnscaledElasticities();
                //cout<<"Unscaled Flux: "<<pMethodMCA->getUnscaledFluxCC()<<endl;
                //cout<<"Unscaled Concentration"<<pMethodMCA->getUnscaledConcentrationCC()<<endl;
                //cerr << "Setting Resolution Error: "<<CCopasiMessage::getAllMessageText() << endl;
                //CCopasiMessage::clearDeque();
            }
            else
            {
                mcaResults = pMethodMCA->getScaledFluxCC();
                mcaConResults = pMethodMCA->getScaledConcentrationCC();
                mcaElaResults = pMethodMCA->getScaledElasticities();
                //cout<<"Scaled Flux"<<pMethodMCA->getScaledFluxCC()<<endl;
                //cout<<"Scaled Concentration"<<pMethodMCA->getScaledConcentrationCC()<<endl;
                //cerr << "Setting Resolution Error: "<<CCopasiMessage::getAllMessageText() << endl;
                //CCopasiMessage::clearDeque();
            }
            for(int j=0; j<nbrReactions; j++)
            {
                os << mcaResults(indexMCAReaction, j);
                //cout<<mcaResults;
                outputMCA<< os.str() << '\t';
                os.str("");
                os.clear();

                os << mcaConResults;
                //cout<<mcaConResults;
                outputMCACon<< os.str() << '\t';
                os.str("");
                os.clear();


            }
            for(int j=0;j<nbrSpecies;j++)
            {
                os << mcaElaResults(indexMCAReaction, j);
                //cout<<mcaElaResults;
                outputMCAEla<< os.str() << '\t';
                os.str("");
                os.clear();
            }
            speciesNames = "";
            for(int j=0; j<nbrSpecies; j++)
            {
                speciesNames += pDataModel->getModel()->getMetabolites()[j]->getObjectName() + '\t';
            }
            outputMCA << '\n';
            outputMCAEla << '\n';
            outputMCACon << '\n';


        }
        else
        {
            cout << "MCA not performed" << endl;
        }

    }




    cout<<"Total Steady-State found: "<<foundSS<<" out of "<<nbrIterations<<endl;

    CCopasiRootContainer::destroy();
    outputSSConcentration.close();
    outputSSConcentrationRate.close();
    outputSSReaction.close();
    outputParam.close();
    outputMCA.close();
    cout << "End of the program." << endl;

    return 0;
}
