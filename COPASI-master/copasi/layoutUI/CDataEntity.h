// Begin CVS Header
//   $Source: /Volumes/Home/Users/shoops/cvs/copasi_dev/copasi/layoutUI/CDataEntity.h,v $
//   $Revision: 1.10 $
//   $Name:  $
//   $Author: shoops $
//   $Date: 2011/03/07 19:29:16 $
// End CVS Header

// Copyright (C) 2011 - 2010 by Pedro Mendes, Virginia Tech Intellectual
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

#ifndef CDATAENTITY_H
#define CDATAENTITY_H

#include <string>
#include <map>

#include "copasi.h"

class CDataEntity
{
private:
  C_FLOAT64 mTime;
  std::map<std::string, C_FLOAT64> mSpeciesValueMap;
  std::map<std::string, C_FLOAT64> mOrigValueMap;

public:
  CDataEntity();
  CDataEntity(C_FLOAT64 t);
  CDataEntity(const CDataEntity & src);

  void putValueForSpecies(std::string nodeKey, C_FLOAT64 value);
  void putOrigValueForSpecies(std::string nodeKey, C_FLOAT64 value);

  C_FLOAT64 getValueForSpecies(std::string nodeKey);
  C_FLOAT64 getOrigValueForSpecies(std::string nodeKey);
  size_t getNumberOfElements() {return mSpeciesValueMap.size();}
  C_FLOAT64 getTime() {return mTime;}
  friend std::ostream & operator<<(std::ostream &os, const CDataEntity & gn);
};
#endif
