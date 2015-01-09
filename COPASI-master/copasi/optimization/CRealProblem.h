/* Begin CVS Header
   $Source: /Volumes/Home/Users/shoops/cvs/copasi_dev/copasi/optimization/CRealProblem.h,v $
   $Revision: 1.12 $
   $Name:  $
   $Author: shoops $
   $Date: 2012/04/23 21:11:20 $
   End CVS Header */

// Copyright (C) 2012 by Pedro Mendes, Virginia Tech Intellectual
// Properties, Inc., University of Heidelberg, and The University
// of Manchester.
// All rights reserved.

// Copyright (C) 2001 - 2007 by Pedro Mendes, Virginia Tech Intellectual
// Properties, Inc. and EML Research, gGmbH.
// All rights reserved.

/**
 *  File name: CRealProblem.h
 *
 *  Programmer: Yongqun He
 *  Contact email: yohe@vt.edu
 *  Purpose: This is the header file of the CRealProblem class.
 *           It specifies the optimization problem with its own members and
 *           functions. It's used by COptAlgorithm class and COptimization class
 *           It's a subclass of COptProblem
 *
 *  The specific purpose here is to use it to demonstrate that the system works
 *  fine. It is for my independent study project.  -- Yongqun He
 *
 *  The function I use here is a N-Dimensional Test Function:
 *            f(x) = (1/2)Sum(j=1, n)(Xj^4 - 16Xj^2 + 5Xj)
 *        where, x = [X1, X2, ... , Xj, ..., Xn]
 *        Number of global minima = 1;
 *    Global minimum found by TRUST is:
 *        [-2.90354, -2.90354, ..., -2.90354].
 *
 */

#ifndef CREALPROBLEM_H
#define CREALPROBLEM_H

#include <string>
#include <vector>

#include "optimization/COptProblem.h"

//
/** @dia:pos 3.05,34.05 */
/** @dia:route COptProblem; v,9.3288,30.75,32.7,9.3472,34.05 */
class CRealProblem : public COptProblem
{
  //data member
  //no real variables defined here, but the base class has many members
private:

  // Implementation
public:

  /**
   * Default constructor
   */
  CRealProblem();

  /**
   * Destructor
   */
  ~CRealProblem();

  /**
   * calculate function for optimization
   */
  virtual bool calculate();
};

#endif  // the end
