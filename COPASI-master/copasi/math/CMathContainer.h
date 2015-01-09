// Copyright (C) 2011 - 2013 by Pedro Mendes, Virginia Tech Intellectual
// Properties, Inc., University of Heidelberg, and The University
// of Manchester.
// All rights reserved.

#ifndef COPASI_CMathContainer
#define COPASI_CMathContainer

#include <map>

#include "copasi/report/CCopasiContainer.h"

#include "copasi/math/CMathObject.h"
#include "copasi/math/CMathEvent.h"
#include "copasi/math/CMathDependencyGraph.h"
#include "copasi/math/CMathReaction.h"

#include "copasi/utilities/CVector.h"

#include "copasi/model/CModelParameter.h"

class CModelEntity;
class CReaction;
class CMoiety;

template < class CType > class CCopasiVector;

class CMathContainer: public CCopasiContainer
{
private:

  struct sDiscontinuous
  {
  public:
    CMathEventN * pEvent;
    CMathObject * pDiscontinuous;
  };

  /**
   * Default Constructor
   */
  CMathContainer();

public:

  /**
   * Specific Constructor
   * @param CModel & model
   */
  CMathContainer(CModel & model);

  /**
   * Copy constructor
   * @param const CMathContainer & src
   */
  CMathContainer(const CMathContainer & src);

  /**
   * Destructor
   */
  virtual ~CMathContainer();

  /**
   * Retrieve the values of all mathematical objects
   * @return const CVectorCore< C_FLOAT64 > & values
   */
  const CVectorCore< C_FLOAT64 > & getValues() const;

  /**
   * Retrieve the values of all mathematical objects
   * @return CVectorCore< C_FLOAT64 > & values
   */
  CVectorCore< C_FLOAT64 > & getValues();

  /**
   * Set the values of all mathematical objects
   * @param const CVectorCore< C_FLOAT64 > & values
   */
  void setValues(const CVectorCore< C_FLOAT64 > & values);

  /**
   * Retrieves the initial state values, i.e., all initial values of objects of
   * simulation type Fixed, Time, ODE, Dependent, and Independent. It includes only
   * extensive values for species.
   * @return const CVectorCore< C_FLOAT64 > & initialState
   */
  const CVectorCore< C_FLOAT64 > & getInitialState() const;

  /**
   * Retrieves the initial state values, i.e., all initial values of objects of
   * simulation type Fixed, Time, ODE, Dependent, and Independent. It includes only
   * extensive values for species.
   * @return CVectorCore< C_FLOAT64 > & initialState
   */
  CVectorCore< C_FLOAT64 > & getInitialState();

  /**
   * Set the initial state values, i.e., all initial values of objects of
   * simulation type Fixed, Time, ODE, Dependent, and Independent. It includes only
   * extensive values for species.
   * @param const CVectorCore< C_FLOAT64 > & initialState
   */
  void setInitialState(const CVectorCore< C_FLOAT64 > & initialState);

  /**
   * Retrieves the state values, i.e., all initial values of objects of
   * simulation type Time, ODE, Dependent, and Independent. It includes only
   * extensive values for species.
   * @return const CVectorCore< C_FLOAT64 > & state
   */
  const CVectorCore< C_FLOAT64 > & getState() const;

  /**
   * Retrieves the state values, i.e., all initial values of objects of
   * simulation type Time, ODE, Dependent, and Independent. It includes only
   * extensive values for species.
   * @return const CVectorCore< C_FLOAT64 > & state
   */
  CVectorCore< C_FLOAT64 > & getState();

  /**
   * Set the state values, i.e., all initial values of objects of
   * simulation type Time, ODE, Dependent, and Independent. It includes only
   * extensive values for species.
   * @param const CVectorCore< C_FLOAT64 > & state
   */
  void setState(const CVectorCore< C_FLOAT64 > & state);

  /**
   * Retrieves the state values, i.e., all initial values of objects of
   * simulation type Time, ODE, and Dependent. It includes only
   * extensive values for species.
   * @return const CVectorCore< C_FLOAT64 > & stateReduced
   */
  const CVectorCore< C_FLOAT64 > & getStateReduced() const;

  /**
   * Retrieves the reduced state values, i.e., all initial values of objects of
   * simulation type Time, ODE, and Dependent. It includes only
   * extensive values for species.
   * @return const CVectorCore< C_FLOAT64 > & stateReduced
   */
  CVectorCore< C_FLOAT64 > & getStateReduced();

  /**
   * Set the reduced state values, i.e., all initial values of objects of
   * simulation type Time, ODE, Dependent, and Independent. It includes only
   * extensive values for species.
   * @param const CVectorCore< C_FLOAT64 > & stateReduced
   */
  void setStateReduced(const CVectorCore< C_FLOAT64 > & stateReduced);

  /**
   * Calculate all dependent initial values based on initial extensive
   * or intensive values
   * @param const CModelParameter::Framework & framework
   */
  void updateInitialValues(const CModelParameter::Framework & framework);

  /**
   * Set the transient values to the initial values and calculate all
   * dependent values. Please call updateInitialValues before calling this
   * method to assure that all values are up to date.
   */
  void applyInitialValues();

  /**
   * Calculate all values required for simulation based on the current state
   * @param const bool & useMoieties
   */
  void updateSimulatedValues(const bool & useMoieties);

  /**
   * Apply the given update sequence to the mathematical objects in the container
   * @param const CObjectInterface::UpdateSequence & updateSequence
   */
  void applyUpdateSequence(const CObjectInterface::UpdateSequence & updateSequence);

  /**
   * Fetch the initial state from the associated model
   */
  void fetchInitialState();

  /**
   * Push the initial state to the associated model
   */
  void pushInitialState();

  /**
   * Fetch the state from the associated model
   */
  void fetchState();

  /**
   * Push the state to the associated model
   */
  void pushState();

  /**
   * Retrieve the CN of the math container
   * The math container provides values for the numerical values of model objects.
   * For the CN mechanism to work properly it has to pretend to be the model.
   * @return CCopasiObjectName
   */
  virtual CCopasiObjectName getCN() const;

  /**
   * Retrieve a descendant object by its CN.
   */
  virtual const CObjectInterface * getObject(const CCopasiObjectName & cn) const;

  /**
   * Retrieve a pointer to corresponding the mathematical object
   * for the given abject if it exists otherwise NULL
   * @param const CObjectInterface * pObject
   * @return CMathObject * pMathObject
   */
  CMathObject * getMathObject(const CObjectInterface * pObject) const;

  /**
   * Retrieve a pointer to the mathematical object for data value pointer
   * @param const C_FLOAT64 * pDataValue
   * @return CMathObject * pMathObject
   */
  CMathObject * getMathObject(const C_FLOAT64 * pDataValue) const;

  /**
   * Retrieve a pointer to the mathematical object for the given CN
   * @param const CCopasiObjectName & cn
   * @return CMathObject * pMathObject
   */
  CMathObject * getMathObject(const CCopasiObjectName & cn) const;

  /**
   * Retrieve the model for the container
   * @return const CModel & model
   */
  const CModel & getModel() const;

  /**
   * Copy a node and all its children. Nodes are converted to suite the math container,
   * i.e., objects nodes point to math object or numbers, function calls are expanded
   * and discrete nodes are replaced by object node pointing to newly created math objects,
   * which are and targets of automatically generated events.
   * @param const CEvaluationNode * pSrc
   * @param const bool & replaceDiscontinuousNodes
   * @return CEvaluationNode * pCopy
   */
  CEvaluationNode * copyBranch(const CEvaluationNode * pSrc,
                               const bool & replaceDiscontinuousNodes);

  /**
   * Copy a node and all its children. Nodes are converted to suite the math container,
   * i.e., objects nodes point to math object or numbers, function calls are expanded
   * and discrete nodes are replaced by object node pointing to newly created math objects,
   * which are and targets of automatically generated events.
   * @param const CEvaluationNode * pSrc
   * @param const CMath::Variables & variables
   * @param const bool & replaceDiscontinuousNodes
   * @return CEvaluationNode * pCopy
   */
  CEvaluationNode * copyBranch(const CEvaluationNode * pSrc,
                               const CMath::Variables< CEvaluationNode * > & variables,
                               const bool & replaceDiscontinuousNodes);

  /**
   * Replace a discontinuous node by an object node pointing to newly created math objects,
   * which are targets of automatically generated events.
   * @param const CEvaluationNode * pSrc
   * @param const std::vector< CEvaluationNode * > & children
   * @return CEvaluationNode * pCopy
   */
  CEvaluationNode * replaceDiscontinuousNode(const CEvaluationNode * pSrc,
      const std::vector< CEvaluationNode * > & children);

  /**
   * Retrieve the pointer to the corresponding initial value
   * @param const C_FLOAT64 * pValue
   * @return C_FLOAT64 * pInitialvalue
   */
  C_FLOAT64 * getInitialValuePointer(const C_FLOAT64 * pValue) const;

private:
  /**
   * Initialize the mathematical model
   */
  void init();

  /**
   * Allocate the memory for objects and values
   */
  void allocate();

  /**
   * Initialize the pointers
   * @param sPointers & pointers
   */
  void initializePointers(CMath::sPointers & pointers);

#ifdef COPASI_DEBUG
  /**
   * Print pointer state
   * @param sPointers & pointers
   */
  void printPointers(CMath::sPointers & pointers);
#endif // COAPSI_DEBUG

  /**
   * Initialize the pointers used for the conversion of
   * discontinuous nodes.
   */
  void initializeDiscontinuousCreationPointer();

  /**
   * Initialize the objects
   * @param sPointers & pointers
   */
  void initializeObjects(CMath::sPointers & pointers);

  /**
   * Initialize the events
   * @param sPointers & pointers
   */
  void initializeEvents(CMath::sPointers & pointers);

  /**
   * Compile the objects
   */
  bool compileObjects();

  /**
   * Compile the events
   */
  bool compileEvents();

  /**
   * Convert nodes which generate discrete changes into events
   * and compiles them.
   */
  bool compileDiscreteChanges();

  /**
   * Create a node based on the given pointer to the math object
   * @param const CObjectInterface * pMathObject
   * @return CEvaluationNode * pNode
   */
  CEvaluationNode * createNodeFromObject(const CObjectInterface * pMathObject);

  /**
   * Create a node based on the given pointer to a data value
   * @param const C_FLOAT64 * pDataValue
   * @return CEvaluationNode * pNode
   */
  CEvaluationNode * createNodeFromValue(const C_FLOAT64 * pDataValue);

  /**
   * Create Dependency Graphs
   */
  void createDependencyGraphs();

  /**
   * Create the update sequences needed to synchronize the initial values
   */
  void createSynchronizeInitialValuesSequence();

  /**
   * Create the update sequence used when applying the initial state
   */
  void createApplyInitialValuesSequence();

  /**
   * Create the update sequences used to calculate all values required for simulation
   */
  void createUpdateSimulationValuesSequence();

  /**
   * Determine the entity type of an entity
   * @param const CModelEntity * pEntity
   * @return CMath::EntityType entityType
   */
  static CMath::EntityType getEntityType(const CModelEntity * pEntity);

  /**
   * Initialize several mathematical objects of a common simulation type and
   * advance relevant pointers
   * @param const std::vector<const CModelEntity*> & entities
   * @param const CMath::SimulationType & simulationType
   * @param sPointers & pointers
   */
  void initializeMathObjects(const std::vector<const CModelEntity*> & entities,
                             const CMath::SimulationType & simulationType,
                             CMath::sPointers & pointers);

  /**
   * Initialize several mathematical objects for local reaction parameters and
   * advance relevant pointers
   * @param const std::vector<const CCopasiObject *> & parameters
   * @param CMathContainer::sPointers & p
   */
  void initializeMathObjects(const std::vector<const CCopasiObject *> & parameters,
                             CMath::sPointers & p);

  /**
   * Initialize several mathematical objects for local reaction parameters and
   * advance relevant pointers
   * @param const CCopasiVector< CReaction > & reactions
   * @param CMathContainer::sPointers & p
   */
  void initializeMathObjects(const CCopasiVector< CReaction > & reactions,
                             CMath::sPointers & p);

  /**
   * Initialize several mathematical objects for local reaction parameters and
   * advance relevant pointers
   * @param const CCopasiVector< CMoiety > & moieties
   * @param CMathContainer::sPointers & p
   */
  void initializeMathObjects(const CCopasiVector< CMoiety > & moieties,
                             CMath::sPointers & p);

  /**
   * Determine whether on object has calculation dependencies.
   */
  static bool hasDependencies(const CCopasiObject * pObject);

  /**
   * Map the data object to the math object
   * @param CCopasiObject * pDataObject
   * @param CMathObject * pMathObject
   */
  void map(CCopasiObject * pDataObject, CMathObject * pMathObject);

  /**
   * Create an event of type CEvent::Discontinuity for each discontinuity in the model
   */
  void createDiscontinuityEvents();

  /**
   * Create an event of type CEvent::Discontinuity for each discontinuity the tree with root pNode
   * @param const CEvaluationTree * pTree
   */
  void createDiscontinuityEvents(const CEvaluationTree * pTree);

  /**
   * Create an event of type CEvent::Discontinuity for the given node
   * @param const CEvaluationNode * pNode
   */
  void createDiscontinuityDataEvent(const CEvaluationNode * pNode);

  /**
   * Create the infix for trigger of the event which tracks changes in a discontinuity
   * represented by the node.
   * @param const CEvaluationNode * pNode
   */
  std::string createDiscontinuityTriggerInfix(const CEvaluationNode * pNode);

  // Attributes

  /**
   * A pointer to the data model which mathematics are contained
   */
  CModel * mpModel;
  const CObjectInterface * mpAvogadro;
  const CObjectInterface * mpQuantity2NumberFactor;

  CVector< C_FLOAT64 > mValues;

  CVectorCore< C_FLOAT64 > mInitialExtensiveValues;
  CVectorCore< C_FLOAT64 > mInitialIntensiveValues;
  CVectorCore< C_FLOAT64 > mInitialExtensiveRates;
  CVectorCore< C_FLOAT64 > mInitialIntensiveRates;
  CVectorCore< C_FLOAT64 > mInitialParticleFluxes;
  CVectorCore< C_FLOAT64 > mInitialFluxes;
  CVectorCore< C_FLOAT64 > mInitialTotalMasses;
  CVectorCore< C_FLOAT64 > mInitialEventTriggers;

  CVectorCore< C_FLOAT64 > mExtensiveValues;
  CVectorCore< C_FLOAT64 > mIntensiveValues;
  CVectorCore< C_FLOAT64 > mExtensiveRates;
  CVectorCore< C_FLOAT64 > mIntensiveRates;
  CVectorCore< C_FLOAT64 > mParticleFluxes;
  CVectorCore< C_FLOAT64 > mFluxes;
  CVectorCore< C_FLOAT64 > mTotalMasses;
  CVectorCore< C_FLOAT64 > mEventTriggers;

  CVectorCore< C_FLOAT64 > mEventDelays;
  CVectorCore< C_FLOAT64 > mEventPriorities;
  CVectorCore< C_FLOAT64 > mEventAssignments;
  CVectorCore< C_FLOAT64 > mEventRoots;
  CVectorCore< C_FLOAT64 > mEventRootStates;
  CVectorCore< C_FLOAT64 > mPropensities;
  CVectorCore< C_FLOAT64 > mDependentMasses;
  CVectorCore< C_FLOAT64 > mDiscontinuous;

  size_t mFixedCount;
  size_t mEventTargetCount;
  size_t mODECount;
  size_t mIndependentCount;
  size_t mDependentCount;
  size_t mAssignmentCount;

  /**
   * The initial state contains also all fixed values
   */
  CVectorCore< C_FLOAT64 > mInitialState;

  /**
   * The state contains values of type Time, ODE, Independent, and Dependent
   */
  CVectorCore< C_FLOAT64 > mState;

  /**
   * The state contains values of type Time, ODE, Independent
   */
  CVectorCore< C_FLOAT64 > mStateReduced;

  /**
   * Dependency graph for initial value calculations
   */
  CMathDependencyGraph mInitialDependencies;

  /**
   * Dependency graph for transient value calculations
   */
  CMathDependencyGraph mTransientDependencies;

  /**
   * The sequence of updates needed to synchronize the initial values based
   * on extensive values, i.e., species amounts
   */
  CObjectInterface::UpdateSequence mSynchronizeInitialValuesSequenceExtensive;

  /**
   * The sequence of updates needed to synchronize the initial values based
   * on intensive values, i.e., species concentrations.
   */
  CObjectInterface::UpdateSequence mSynchronizeInitialValuesSequenceIntensive;

  /**
   * The sequence of updates needed to apply the initial values
   */
  CObjectInterface::UpdateSequence mApplyInitialValuesSequence;

  /**
   * The sequence of updates needed to calculate all simulation required values based
   * on the assumption that all state values may have changed
   */
  CObjectInterface::UpdateSequence mSimulationValuesSequence;

  /**
   * The sequence of updates needed to calculate all simulation required values based
   * on the assumption that all state values may have changed
   */
  CObjectInterface::UpdateSequence mSimulationValuesSequenceReduced;

  /**
   * The set of objects which determine the initial state of the model based on extensive
   * values
   */
  CObjectInterface::ObjectSet mInitialStateValueExtensive;

  /**
   * The set of objects which determine the initial state of the model based on intensive
   * values
   */
  CObjectInterface::ObjectSet mInitialStateValueIntensive;

  /**
   * The set of objects which determine the transient state of the model
   */
  CObjectInterface::ObjectSet mStateValues;

  /**
   * The set of objects which determine the transient state of the reduced model
   */
  CObjectInterface::ObjectSet mReducedStateValues;

  /**
   * The objects which are required to be up to date for simulation of the full model,
   * i.e., the right hand side of ODEs, rates of species determined by reaction,
   * and event roots.
   */
  CObjectInterface::ObjectSet mSimulationRequiredValues;

  /**
   * The objects which are required to be up to date for simulation of the reduced model,
   * i.e., the right hand side of ODEs, rates of independent species determined by reaction,
   * and event roots.
   */
  // CObjectInterface::ObjectSet mSimulationRequiredValuesReduced;

  /**
   * A vector containing all math objects.
   */
  CVector< CMathObject > mObjects;

  /**
   * A vector containing all math events.
   */
  CVector< CMathEventN > mEvents;

  /**
   * A vector containing all math reactions.
   */
  CVector< CMathReaction > mReactions;

  /**
   * Structure of pointers used for creating discontinuities.
   */
  sDiscontinuous mCreateDiscontinuousPointer;

  /**
   * A map from data objects to math objects
   */
  std::map< CCopasiObject *, CMathObject * > mDataObject2MathObject;

  /**
   * A map from data objects values to math objects
   */
  std::map< C_FLOAT64 *, CMathObject * > mDataValue2MathObject;

  /**
   * A vector of data events for discontinuities
   */
  CCopasiVector< CEvent > mDiscontinuityEvents;

  /**
   * A map from the infix of the expression of a discontinuity to the object
   * representing it.
   */
  std::map< std::string, CMathObject * > mDiscontinuityInfix2Object;

  /**
   * A map from the discontinuity trigger infix of to event
   * representing it.
   */
  std::map< std::string,  CMathEventN * > mTriggerInfix2Event;
};

#endif // COPASI_CMathContainer
