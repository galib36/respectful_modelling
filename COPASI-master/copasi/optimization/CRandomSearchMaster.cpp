/* Begin CVS Header
   $Source: /Volumes/Home/Users/shoops/cvs/copasi_dev/copasi/optimization/CRandomSearchMaster.cpp,v $
   $Revision: 1.4 $
   $Name:  $
   $Author: shoops $
   $Date: 2006/04/27 01:29:53 $
   End CVS Header */

// Copyright � 2005 by Pedro Mendes, Virginia Tech Intellectual
// Properties, Inc. and EML Research, gGmbH.
// All rights reserved.

/***************************************************************************
                    CRandomSearchMaster.cpp  -  Random Optimizer
                       -------------------


Programmer           : Rohan Luktuke
email                : rluktuke@vt.edu
 ***************************************************************************/

/***************************************************************************
 * This is the implementation of the Random Algorithm for Optimization.  The
 * class is inherited from the COptAlgorithm class
 ***************************************************************************/

#include "copasi.h"
#include "COptMethod.h"

CRandomSearchMaster::CRandomSearchMaster():
    COptMethod(CCopasiMethod::RandomSearchMaster)
{}

CRandomSearchMaster::CRandomSearchMaster(const CRandomSearchMaster & src):
    COptMethod(src)
{}

/**
 * Destructor
 */
CRandomSearchMaster::~CRandomSearchMaster(){}

/**
 * Optimizer Function
 * Returns: nothing
 */
C_INT32 CRandomSearchMaster::optimise()
{return 0;}
