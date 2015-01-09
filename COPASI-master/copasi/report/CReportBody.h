/* Begin CVS Header
   $Source: /Volumes/Home/Users/shoops/cvs/copasi_dev/copasi/report/CReportBody.h,v $
   $Revision: 1.8 $
   $Name:  $
   $Author: gauges $
   $Date: 2006/10/15 08:00:56 $
   End CVS Header */

// Copyright � 2005 by Pedro Mendes, Virginia Tech Intellectual
// Properties, Inc. and EML Research, gGmbH.
// All rights reserved.

/****************************************************************************
 *
 * Created: Tue Aug 11 2003
 *      by: Liang Xu
 *
 *  Created for Copasi by Liang Xu
 ****************************************************************************/

#ifndef COPASI_CREPORT_BODY
#define COPASI_CREPORT_BODY

#include <string>
#include <vector>
#include "report/CCopasiObjectName.h"

class CReportTable;
class CCopasiObjectName;

class CReportBody
  {
  private:
    std::vector<CCopasiObjectName> mReportObjectName;
    CReportTable *mpReportTable;

  public:
    /**
        *  Default constructor.
        */
    CReportBody();

    /**
    cleanup
    */
    void cleanup();

    /**
    *  Destructor.
    */
    ~CReportBody();

    /**
    gets the Table tag  in the Body Tag
    */
    CReportTable* getTable();

    /**
    sets the Table tag in the Body tag
    */
    void setTable(CReportTable *table);
  };

class CReportTable
  {
  private:
    std::vector<CCopasiObjectName>* mpObjectNameVector;
    std::string mSeparator;
    bool mPrintTitle;

  public:
    /**
        *  Default constructor.
        */
    CReportTable();

    /**
    cleanup
    */
    void cleanup();

    /**
       *  Destructor.
       */
    ~CReportTable();

    /**
    gets the Object tag of the Table tag
    */
    inline std::vector<CCopasiObjectName>* getObjectNameVector()
    {
      return mpObjectNameVector;
    }

    /**
       sets the Object tag of the Table tag
    */
    inline void setObjectNameVector(std::vector<CCopasiObjectName> *objectVector)
    {
      mpObjectNameVector = objectVector;
    }

    /**
    gets the Separator attribute of the Table tag
    */
    inline const std::string& getSeparator()
    {
      return mSeparator;
    }

    /**
       sets the separator attribute of the Table tag
    */
    inline void setSeparator(const std::string separator)
    {
      mSeparator = separator;
    }

    /**
    gets the printTitle of the Table tag
    */
    bool printTitle()
    {
      return mPrintTitle;
    }

    /**
       sets the print Title of the Table tag
    */
    void setPrintTitle(bool printTitle)
    {
      mPrintTitle = printTitle;
    }
  };

#endif
