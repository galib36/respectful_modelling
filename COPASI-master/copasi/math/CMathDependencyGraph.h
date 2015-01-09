// Copyright (C) 2011 - 2013 by Pedro Mendes, Virginia Tech Intellectual
// Properties, Inc., University of Heidelberg, and The University
// of Manchester.
// All rights reserved.

#ifndef COPASI_CMathDependencyGraph
#define COPASI_CMathDependencyGraph

#include <map>
#include <set>
#include <vector>

#include "copasi/report/CCopasiObject.h"
#include "copasi/math/CMathEnum.h"

class CMathDependencyNode;

class CMathDependencyGraph
{
public:
  typedef std::map< const CObjectInterface *, CMathDependencyNode * > NodeMap;
  typedef NodeMap::iterator iterator;
  typedef NodeMap::const_iterator const_iterator;

  // Operations
  /**
   * Constructor
   */
  CMathDependencyGraph();

  /**
   * Destructor
   */
  ~CMathDependencyGraph();

  /**
   * Clear the whole dependency tree
   */
  void clear();

  /**
   * Add an object and all its prerequisites to the tree
   * @param const CObjectInterface * pObject
   * @return std::map< const CObjectInterface *, CMathDependencyNode * >::iterator itObject
   */
  iterator addObject(const CObjectInterface * pObject);

  /**
   * Construct a update sequence for the given context
   * @param const CMath::SimulationContextFlag & context
   * @param const CObjectInterface::ObjectSet & changedObjects
   * @param const CObjectInterface::ObjectSet & requestedObjects
   * @param CObjectInterface::UpdateSequence & updateSequence)
   * @return bool success
   */
  bool getUpdateSequence(const CMath::SimulationContextFlag & context,
                         const CObjectInterface::ObjectSet & changedObjects,
                         const CObjectInterface::ObjectSet & requestedObjects,
                         CObjectInterface::UpdateSequence & updateSequence);

  void exportDOTFormat(std::ostream & os, const std::string & name) const;

private:
  std::string getDOTNodeId(const CObjectInterface * pObject) const;

  // Attributes
  NodeMap mObjects2Nodes;

  mutable std::map< const CObjectInterface *, size_t > mObject2Index;
};

#endif // COPASI_CMathDependencyGraph
