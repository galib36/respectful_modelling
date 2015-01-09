/* Begin CVS Header
   $Source: /Volumes/Home/Users/shoops/cvs/copasi_dev/copasi/optimization/CRandomSearchMaster.h,v $
   $Revision: 1.6 $
   $Name:  $
   $Author: shoops $
   $Date: 2006/04/27 01:29:53 $
   End CVS Header */

// Copyright � 2005 by Pedro Mendes, Virginia Tech Intellectual
// Properties, Inc. and EML Research, gGmbH.
// All rights reserved.

/**
 * CRandomSearchMaster class
 */

#ifndef COPASI_CRandomSearchMaster
#define COPASI_CRandomSearchMaster

/** @dia:pos 18.05,34.05 */
/** @dia:route COptMethod; v,46.9608,16.35,33,32.3332,34.05 */
class CRandomSearchMaster : public COptMethod
  {
    friend COptMethod * COptMethod::createMethod(CCopasiMethod::SubType subType);

    // Operations
  private:
    /**
     * Default Constructor
     */
    CRandomSearchMaster();

  public:
    /**
     * Copy Constructor
     * @param const CRandomSearchMaster & src
     */
    CRandomSearchMaster(const CRandomSearchMaster & src);

    /**
     * Destructor
     */
    virtual ~CRandomSearchMaster();

    /**
     * Execute the optimization algorithm calling simulation routine
     * when needed. It is noted that this procedure can give feedback
     * of its progress by the callback function set with SetCallback.
     */
    virtual bool optimise();
  };

#endif  // COPASI_CRandomSearchMaster
