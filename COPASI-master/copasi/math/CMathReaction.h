// Copyright (C) 2013 by Pedro Mendes, Virginia Tech Intellectual
// Properties, Inc., University of Heidelberg, and The University
// of Manchester.
// All rights reserved.

#ifndef COPASI_CMathReaction
#define COPASI_CMathReaction

#include <utility>

#include "copasi/report/CCopasiObject.h"
#include "copasi/utilities/CVector.h"

class CReaction;
class CMathObject;
class CMathContainer;

class CMathReaction
{
public:
  /**
   * Constructor
   */
  CMathReaction();

  /**
   * Destructor
   */
  ~CMathReaction();

  /**
   * Initialize the reaction from the model reaction in the given container context
   * @param const CReaction * pReaction
   * @param CMathContainer & container
   */
  void initialize(const CReaction * pReaction, CMathContainer & container);

  /**
   * Copy an existing object
   * @param const CMathReaction & src
   * @param CMathContainer & container
   * @param const size_t & valueOffset
   * @param const size_t & objectOffset
   */
  void copy(const CMathReaction & src, CMathContainer & container, const size_t & valueOffset, const size_t & objectOffset);

  /**
   * Fire the reaction count times. Note count must not necessarily be an integer
   * @param const C_FLOAT64 & count
   */
  void fire(const C_FLOAT64 & count);

  /**
   * Retrieve a pointer to the mathematical object for the particle flux.
   * @return const CMathObject * pParticleFluxObject
   */
  const CMathObject * getParticleFluxObject() const;

  /**
   * Retrieve a pointer to the mathematical object for the flux.
   * @return const CMathObject * pFluxObject
   */
  const CMathObject * getFluxObject() const;

  /**
   * Retrieve a pointer to the mathematical object for the propensity.
   * @return const CMathObject * pPropensityObject
   */
  const CMathObject * getPropensityObject() const;

  /**
   * Retrieve the set of modified species
   * @return const CObjectInterface::ObjectSet & modifiedSpecies
   */
  const CObjectInterface::ObjectSet & getModifiedSpecies() const;

  /**
   * Retrieve a pointer to the model reaction
   * @return const CReaction * pModelReaction
   */
  const CReaction * getModelReaction() const;

private:
  /**
   * A pointer to model reaction.
   */
  const CReaction * mpReaction;

  /**
   * A pointer to the mathematical particle flux.
   */
  CMathObject * mpParticleFlux;

  /**
   * A pointer to the mathematical flux.
   */
  CMathObject * mpFlux;

  /**
   * A pointer to the mathematical propensity.
   */
  CMathObject * mpPropensity;

  /**
   * The set of modified species
   */
  CObjectInterface::ObjectSet mModifiedSpecies;

  /**
   * Information for updating the species particle numbers when the
   * reaction fires.
   */
  CVector< std::pair< C_FLOAT64, C_FLOAT64 * > > mStepUpdates;
};

#endif // COPASI_CMathReaction
