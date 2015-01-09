// Begin CVS Header
//   $Source: /Volumes/Home/Users/shoops/cvs/copasi_dev/copasi/oscillation/COscillationTask.h,v $
//   $Revision: 1.3 $
//   $Name:  $
//   $Author: shoops $
//   $Date: 2011/03/07 19:32:04 $
// End CVS Header

// Copyright (C) 2011 - 2010 by Pedro Mendes, Virginia Tech Intellectual
// Properties, Inc., University of Heidelberg, and The University
// of Manchester.
// All rights reserved.

// Copyright (C) 2008 by Pedro Mendes, Virginia Tech Intellectual
// Properties, Inc., EML Research, gGmbH, University of Heidelberg,
// and The University of Manchester.
// All rights reserved.

/**
 * COscillationTask class.
 *
 * This class implements an oscillation detection task
 */

#ifndef COPASI_COscillationTask
#define COPASI_COscillationTask

#include "utilities/CCopasiTask.h"

class COscillationProblem;
class COscillationMethod;
class CReport;

class COscillationTask : public CCopasiTask
{
  //Attributes

public:
  /**
   * The methods which can be selected for preforming this task.
   */
  static unsigned C_INT32 ValidMethods[];

  /**
   * default constructor
   * @param const CCopasiTask::Type & type (default: optimization)
   * @param const CCopasiContainer * pParent (default: NULL)
   */
  //COscillationTask(const CCopasiTask::Type & type = CCopasiTask::oscillation,
  //                 const CCopasiContainer * pParent = NULL);
  COscillationTask(const CCopasiTask::Type & type ,
                   const CCopasiContainer * pParent);

  /**
   * Copy constructor
   * @param const COscillationTask & src
   */
  COscillationTask(const COscillationTask & src, const CCopasiContainer * pParent = NULL);

  /**
   * Destructor
   */
  ~COscillationTask();

  /**
   * cleanup()
   */
  void cleanup();

  /**
   * Set the call back of the task
   * @param CProcessReport * pCallBack
   * @result bool succes
   */
  virtual bool setCallBack(CProcessReport * pCallBack);

  /**
   * Initialize the task. If an ostream is given this ostream is used
   * instead of the target specified in the report. This allows nested
   * tasks to share the same output device.
   * @param const OutputFlag & of
   * @param COutputHandler * pOutputHandler
   * @param std::ostream * pOstream (default: NULL)
   * @return bool success
   */
  virtual bool initialize(const OutputFlag & of,
                          COutputHandler * pOutputHandler,
                          std::ostream * pOstream);

  /**
   * Process the task with or without initializing to the initial state.
   * @param const bool & useInitialValues
   * @return bool success
   */
  virtual bool process(const bool & useInitialValues);

  /**
   * Set the method type applied to solve the task
   * @param const CCopasiMethod::SubType & type
   * @return bool success
   */
  virtual bool setMethodType(const int & type);

  /**
   * Create a method of the specified type to solve the task.
   * It is the duty of the caller to release the CCopasiMethod.
   * @param const CCopasiMethod::SubType & type
   * @return CCopasiMethod *
   */
  virtual CCopasiMethod * createMethod(const int & type) const;
};
#endif // COPASI_COscillationTask
