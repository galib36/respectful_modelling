// Begin CVS Header
//   $Source: /Volumes/Home/Users/shoops/cvs/copasi_dev/copasi/function/CMassAction.h,v $
//   $Revision: 1.29 $
//   $Name:  $
//   $Author: shoops $
//   $Date: 2012/05/16 23:11:32 $
// End CVS Header

// Copyright (C) 2012 - 2010 by Pedro Mendes, Virginia Tech Intellectual
// Properties, Inc., University of Heidelberg, and The University
// of Manchester.
// All rights reserved.

// Copyright (C) 2001 - 2007 by Pedro Mendes, Virginia Tech Intellectual
// Properties, Inc. and EML Research, gGmbH.
// All rights reserved.

/**
 * CMassAction
 *
 * Created for Copasi by Stefan Hoops
 * (C) Stefan Hoops 2001
 */

#ifndef COPASI_CMassAction
#define COPASI_CMassAction

#include "function/CFunction.h"

/** @dia:pos 62.0314,41.4407 */
/** @dia:route CFunction; v,58.2,35.55,39.4413,82.281,41.4407 */
class CMassAction : public CFunction
{
  // Static Attributes
public:
  static const char * Infix[];

  // Operations
public:
  /**
   * Default constructor
   * @param const std::string & name (default: "NoName")
   * @param const CCopasiContainer * pParent (default: NULL)
   */
  CMassAction(const std::string & name = "NoName",
              const CCopasiContainer * pParent = NULL);

  /**
   * Copy constructor
   * @param "const CFunction &" src
   * @param const CCopasiContainer * pParent (default: NULL)
   */
  CMassAction(const CFunction & src,
              const CCopasiContainer * pParent = NULL);

  /**
   * Specified constructor
   * @param const CCopasiContainer * pParent (default: NULL)
   */
  CMassAction(const TriLogic & reversible,
              const CCopasiContainer * pParent = NULL);

  /**
   *  Destructor
   */
  ~CMassAction();

  /**
   *  Calculates the value of the function
   *  @param const CCallParameters<C_FLOAT64> & callParameters
   *  @return const C_FLOAT64 & result
   */
  virtual const C_FLOAT64 &
  calcValue(const CCallParameters<C_FLOAT64> & callParameters);

  /**
   * Returns whether the function depends on the given parameter
   * &param const void * parameter
   * &param const CCallParameters<C_FLOAT64> & callParameters
   * &param bool dependsOn
   */
  virtual bool dependsOn(const C_FLOAT64 * parameter,
                         const CCallParameters<C_FLOAT64> & callParameters) const;

  /**
   * Set the infix description of the tree and compile it.
   * @param const std::string & infix
   * @return bool success
   */
  virtual bool setInfix(const std::string & infix);

  /**
   * Compile the evaluation tree.
   */
  virtual bool compile();

  virtual std::string writeMathML(const std::vector< std::vector< std::string > > & variables,
                                  bool expand = true, bool fullExpand = true) const;

};

#endif // COPASI_CMassAction
