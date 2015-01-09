// Copyright (C) 2011 - 2014 by Pedro Mendes, Virginia Tech Intellectual
// Properties, Inc., University of Heidelberg, and The University
// of Manchester.
// All rights reserved.

#include <limits>

#include "copasi.h"

#include "CMathEvent.h"
#include "CMathContainer.h"
#include "CMathExpression.h"

#include "model/CMathModel.h"
#include "model/CEvent.h"

#include "function/CFunction.h"
#include "report/CCopasiRootContainer.h"
#include "utilities/CNodeIterator.h"

CMathEventN::CAssignment::CAssignment():
  mpTarget(NULL),
  mpAssignment(NULL)
{}

CMathEventN::CAssignment::~CAssignment()
{}

void CMathEventN::CAssignment::initialize(CMath::sPointers & pointers)
{
  // Initialize the assignment object
  mpAssignment = pointers.pEventAssignmentsObject;
  CMathObject::initialize(pointers.pEventAssignmentsObject, pointers.pEventAssignments,
                          CMath::EventAssignment, CMath::Event, CMath::SimulationTypeUndefined,
                          false, false, NULL);
}

void CMathEventN::CAssignment::copy(const CMathEventN::CAssignment & src,
                                    CMathContainer & /* container */,
                                    const size_t & /* valueOffset */,
                                    const size_t & objectOffset)
{
  mpTarget = src.mpTarget + objectOffset;
  mpAssignment = src.mpAssignment + objectOffset;
}

bool CMathEventN::CAssignment::compile(CEventAssignment * pDataAssignment,
                                       CMathContainer & container)
{
  // A compiled pDataAssignment is prerequisite.
  bool success = true;

  // Determine the target object
  mpTarget = container.getMathObject(pDataAssignment->getTargetObject());

  if (mpTarget != NULL &&
      mpTarget->getSimulationType() == CMath::Fixed)
    {
      mpTarget->setSimulationType(CMath::EventTarget);
    }

  std::vector< CCopasiContainer * > ListOfContainer;

  // Compile the assignment object in the model context
  CExpression AssignmentExpression("AssignmentExpression", &container);
  success &= AssignmentExpression.setInfix(pDataAssignment->getExpression());
  success &= AssignmentExpression.compile(ListOfContainer);
  success &= mpAssignment->setExpression(AssignmentExpression, container);

  return success;
}

void CMathEventN::CAssignment::setTarget(CMathObject * pTarget)
{
  mpTarget = pTarget;
}

void CMathEventN::CAssignment::setExpression(CMathObject * pExpression)
{
  mpAssignment = pExpression;
}

CMathEventN::CTrigger::CRoot::CRoot():
  mpRoot(NULL),
  mpRootState(NULL),
  mEquality(false),
  mDiscrete(false)
{}

CMathEventN::CTrigger::CRoot::~CRoot()
{}

void CMathEventN::CTrigger::CRoot::initialize(CMath::sPointers & pointers)
{
  // Initialize the root object!
  mpRoot = pointers.pEventRootsObject;
  *pointers.pEventRoots = 1.0;
  CMathObject::initialize(pointers.pEventRootsObject, pointers.pEventRoots,
                          CMath::EventRoot, CMath::Event, CMath::SimulationTypeUndefined,
                          false, false, NULL);

  // Initialize the root state object!
  mpRootState = pointers.pEventRootStatesObject;
  *pointers.pEventRootStates = 1.0;
  CMathObject::initialize(pointers.pEventRootStatesObject, pointers.pEventRootStates,
                          CMath::EventRootState, CMath::Event, CMath::SimulationTypeUndefined,
                          false, false, NULL);
}

void CMathEventN::CTrigger::CRoot::copy(const CMathEventN::CTrigger::CRoot & src,
                                        CMathContainer & /* container */,
                                        const size_t & /* valueOffset */,
                                        const size_t & objectOffset)
{
  mpRoot = src.mpRoot + objectOffset;
  mpRootState = src.mpRootState + objectOffset;
  mEquality = src.mEquality;
  mDiscrete = src.mDiscrete;
}

bool CMathEventN::CTrigger::CRoot::compile(CEvaluationNode * pRootNode,
    const bool & equality,
    CMathContainer & container)
{
  mEquality = equality;

  std::vector< CCopasiContainer * > ListOfContainer;
  ListOfContainer.push_back(const_cast< CMathContainer * >(&container));

  CMathExpression * pExpression = new CMathExpression("RootExpression", container);
  bool success = static_cast< CEvaluationTree * >(pExpression)->setRoot(pRootNode);
  success &= mpRoot->setExpressionPtr(pExpression);

  // Compile the root state object
  CMathExpression * pStateExpression = new CMathExpression("RootStateExpression", container);

  CEvaluationNode * pStateExpressionNode = new CEvaluationNodeLogical(CEvaluationNodeLogical::GT, "GT");
  pStateExpressionNode->addChild(new CEvaluationNodeObject((C_FLOAT64 *) mpRoot->getValuePointer()));
  pStateExpressionNode->addChild(new CEvaluationNodeNumber(CEvaluationNodeNumber::DOUBLE, "0.5"));

  success &= static_cast< CEvaluationTree * >(pStateExpression)->setRoot(pStateExpressionNode);
  success &= mpRootState->setExpressionPtr(pStateExpression);

  return success;
}

// static
CEvaluationNode * CMathEventN::CTrigger::CRoot::createTriggerExpressionNode() const
{
  return  new CEvaluationNodeObject((C_FLOAT64 *) mpRootState->getValuePointer());
}

CMathEventN::CTrigger::CTrigger():
  mpTrigger(NULL),
  mpInitialTrigger(NULL),
  mRoots(),
  mDualAction(false)
{}

CMathEventN::CTrigger::~CTrigger()
{}

void CMathEventN::CTrigger::allocate(const CEvent * pDataEvent,
                                     const CMathContainer & container)
{
  // Determine the number of roots.
  CMath::Variables< size_t > Variables;

  std::vector< CCopasiContainer * > Container;
  Container.push_back(const_cast< CMathContainer * >(&container));

  CExpression Trigger("EventTrigger", &container);
  Trigger.setIsBoolean(true);

  if (Trigger.setInfix(pDataEvent->getTriggerExpression()))
    {
      Trigger.compile(Container);
      mRoots.resize(countRoots(Trigger.getRoot(), Variables));
    }
  else
    {
      CFunction TriggerFunction("EventTrigger", &container);

      if (TriggerFunction.setInfix(pDataEvent->getTriggerExpression()))
        {
          TriggerFunction.compile();
          mRoots.resize(countRoots(TriggerFunction.getRoot(), Variables));
        }
    }

  return;
}

void CMathEventN::CTrigger::allocateDiscontinuous(const size_t & nRoots,
    const CMathContainer & /* container */)
{
  mRoots.resize(nRoots);
}

void CMathEventN::CTrigger::initialize(CMath::sPointers & pointers)
{
  // Initialize trigger object.
  mpTrigger = pointers.pEventTriggersObject;
  *pointers.pEventTriggers = 1.0;
  CMathObject::initialize(pointers.pEventTriggersObject, pointers.pEventTriggers,
                          CMath::EventTrigger, CMath::Event, CMath::SimulationTypeUndefined,
                          false, false, NULL);

  // Initialize initial trigger object.
  mpInitialTrigger = pointers.pInitialEventTriggersObject;
  *pointers.pInitialEventTriggers = 1.0;
  CMathObject::initialize(pointers.pInitialEventTriggersObject, pointers.pInitialEventTriggers,
                          CMath::EventTrigger, CMath::Event, CMath::SimulationTypeUndefined,
                          false, true, NULL);

  // Initialize root object.
  CRoot * pRoot = mRoots.array();
  CRoot * pRootEnd = pRoot + mRoots.size();

  for (; pRoot != pRootEnd; ++pRoot)
    {
      pRoot->initialize(pointers);
    }
}

void CMathEventN::CTrigger::copy(const CMathEventN::CTrigger & src,
                                 CMathContainer & container,
                                 const size_t & valueOffset,
                                 const size_t & objectOffset)
{
  mpTrigger = src.mpTrigger + objectOffset;
  mpInitialTrigger = src.mpInitialTrigger + objectOffset;

  mRoots.resize(src.mRoots.size());
  CRoot * pRoot = mRoots.array();
  CRoot * pRootEnd = pRoot + mRoots.size();
  const CRoot * pRootSrc = src.mRoots.array();

  for (; pRoot != pRootEnd; ++pRoot, ++pRootSrc)
    {
      pRoot->copy(*pRootSrc, container, valueOffset, objectOffset);
    }

  mDualAction = src.mDualAction;
}

bool CMathEventN::CTrigger::compile(CEvent * pDataEvent,
                                    CMathContainer & container)
{
  bool success = true;

  CMath::Variables< CEvaluationNode * > Variables;

  std::vector< CCopasiContainer * > ListOfContainer;
  ListOfContainer.push_back(&container);

  CExpression DataTrigger("DataTrigger", &container);
  DataTrigger.setIsBoolean(true);

  if (pDataEvent != NULL)
    {
      DataTrigger.setInfix(pDataEvent->getTriggerExpression());
      mDualAction = false;
    }
  else
    {
      DataTrigger.setInfix(mpTrigger->getExpressionPtr()->getInfix());
    }

  success &= DataTrigger.compile();

  CEvaluationNode * pTriggerRoot = NULL;
  CRoot * pRoot = mRoots.array();

  pTriggerRoot = compile(DataTrigger.getRoot(), Variables, pRoot, container);

  assert(pRoot == mRoots.array() + mRoots.size());

  CMathExpression * pTrigger = new CMathExpression("EventTrigger", container);
  success &= static_cast< CEvaluationTree * >(pTrigger)->setRoot(pTriggerRoot);

  success &= mpTrigger->setExpressionPtr(pTrigger);

  return success;
}

const CVector< CMathEventN::CTrigger::CRoot > & CMathEventN::CTrigger::getRoots() const
{
  return mRoots;
}

void CMathEventN::CTrigger::setExpression(const std::string & infix,
    CMathContainer & container)
{
  assert(mpTrigger != NULL);

  mpTrigger->setExpression(infix, true, container);

  compile(NULL, container);

  std::cout << *mpTrigger;
  std::cout << *mpInitialTrigger;

  CRoot * pRoot = mRoots.array();
  CRoot * pRootEnd = pRoot + mRoots.size();

  for (; pRoot != pRootEnd; ++pRoot)
    {
      std::cout << *pRoot->mpRoot;
      std::cout << *pRoot->mpRootState;
    }
}

// static
size_t CMathEventN::CTrigger::countRoots(const CEvaluationNode * pNode,
    const CMath::Variables< size_t > & variables)
{
  size_t RootCount = 0;

  // TODO CRITICAL We only need to count in boolean functions see compile for details.
  CNodeContextIterator< const CEvaluationNode, std::vector< size_t > > itNode(pNode);
  itNode.setProcessingModes(CNodeIteratorMode::Before | CNodeIteratorMode::After);

  while (itNode.next() != itNode.end())
    {
      if (*itNode == NULL)
        {
          continue;
        }

      switch (itNode.processingMode())
        {
          case CNodeIteratorMode::Before:

            // Variables return always false we need to dig deeper.
            if (CEvaluationNode::type(itNode->getType()) == CEvaluationNode::VARIABLE ||
                itNode->isBoolean())
              {
                continue;
              }

            // We found a non boolean node which does not create a root.
            itNode.skipChildren();
            RootCount = 0;

            break;

          case CNodeIteratorMode::After:
          {
            // We do not need to check whether the root is boolean as non boolean nodes are
            // already processed
            const CEvaluationNode::Type & Type = itNode->getType();

            switch (CEvaluationNode::type(Type))
              {
                case CEvaluationNode::LOGICAL:

                  switch ((int) CEvaluationNode::subType(Type))
                    {
                      case CEvaluationNodeLogical::EQ:
                      case CEvaluationNodeLogical::NE:
                        RootCount = countRootsEQ(*itNode, itNode.context());
                        break;

                      case CEvaluationNodeLogical::LE:
                      case CEvaluationNodeLogical::LT:
                      case CEvaluationNodeLogical::GE:
                      case CEvaluationNodeLogical::GT:
                        RootCount = 1;
                        break;

                      default:
                        RootCount = countRootsDefault(itNode.context());
                        break;
                    }

                  break;

                case CEvaluationNode::CALL:

                  switch ((int) CEvaluationNode::subType(Type))
                    {
                      case CEvaluationNodeCall::FUNCTION:
                      case CEvaluationNodeCall::EXPRESSION:
                        RootCount = countRootsFUNCTION(*itNode, itNode.context());
                        break;

                      default:
                        RootCount = countRootsDefault(itNode.context());
                        break;
                    }

                  break;

                case CEvaluationNode::VARIABLE:

                  switch ((int) CEvaluationNode::subType(Type))
                    {
                      case CEvaluationNodeVariable::ANY:
                        RootCount = countRootsVARIABLE(*itNode, variables);
                        break;

                      default:
                        RootCount = countRootsDefault(itNode.context());
                        break;
                    }

                  break;

                default:
                  RootCount = countRootsDefault(itNode.context());
                  break;
              }
          }
          break;

          default:
            break;
        }

      if (itNode.parentContextPtr() != NULL)
        {
          itNode.parentContextPtr()->push_back(RootCount);
        }
    }

  return RootCount;
}

// static
size_t CMathEventN::CTrigger::countRootsDefault(const std::vector< size_t > & children)
{
  size_t RootCount = 0;

  std::vector< size_t >::const_iterator it = children.begin();
  std::vector< size_t >::const_iterator end = children.end();

  for (; it != end; ++it)
    {
      RootCount += *it;
    }

  return RootCount;
}

// static
size_t CMathEventN::CTrigger::countRootsEQ(const CEvaluationNode * pNode,
    const std::vector< size_t > & children)
{
  size_t nRoots = children[0] + children[1];

  // Equality can be determined between Boolean and double values.
  if (static_cast<const CEvaluationNode *>(pNode->getChild())->isBoolean())
    {
      nRoots = children[0] + children[1];
    }
  else
    {
      nRoots = 2;
    }

  return nRoots;
}

// static
size_t CMathEventN::CTrigger::countRootsFUNCTION(const CEvaluationNode * pNode,
    const std::vector< size_t > & children)
{
  const CEvaluationNode * pTreeRoot =
    static_cast< const CEvaluationNodeCall * >(pNode)->getCalledTree()->getRoot();

  size_t nRoots = countRoots(pTreeRoot, children);

  return nRoots;
}

// static
size_t CMathEventN::CTrigger::countRootsVARIABLE(const CEvaluationNode * pNode,
    const CMath::Variables< size_t > & variables)
{
  size_t Index =
    static_cast< const CEvaluationNodeVariable * >(pNode)->getIndex();

  if (Index == C_INVALID_INDEX)
    {
      return 0;
    }

  return variables[Index];
}

// static
CEvaluationNode * CMathEventN::CTrigger::compile(const CEvaluationNode * pTriggerNode,
    const CMath::Variables< CEvaluationNode * > & variables,
    CMathEventN::CTrigger::CRoot *& pRoot,
    CMathContainer & container)
{
  CEvaluationNode * pNode = NULL;

  CNodeContextIterator< const CEvaluationNode, std::vector< CEvaluationNode * > > itNode(pTriggerNode);
  itNode.setProcessingModes(CNodeIteratorMode::Before | CNodeIteratorMode::After);

  while (itNode.next() != itNode.end())
    {
      if (*itNode == NULL)
        {
          continue;
        }

      switch (itNode.processingMode())
        {
          case CNodeIteratorMode::Before:

            // Variables return always false we need to dig deeper.
            if (CEvaluationNode::type(itNode->getType()) == CEvaluationNode::VARIABLE)
              {
                size_t Index =
                  static_cast< const CEvaluationNodeVariable * >(*itNode)->getIndex();

                if (Index != C_INVALID_INDEX &&
                    Index < variables.size() &&
                    variables[Index]->isBoolean())
                  {
                    continue;
                  }

                // We found a non boolean node which we simply copy.
                itNode.skipChildren();

                // Since a variable may be referred to multiple times we need to copy it.
                if (Index != C_INVALID_INDEX)
                  {
                    pNode = variables[Index]->copyBranch();
                  }
                else
                  {
                    // Variables must not appear in mathematical expressions.
                    // We create an constant node with the variable name and value NaN.
                    pNode = new CEvaluationNodeConstant(CEvaluationNodeConstant::_NaN, itNode->getData());
                  }
              }
            else if (!itNode->isBoolean())
              {
                // We found a non boolean node which we simply copy.
                itNode.skipChildren();
                pNode = container.copyBranch(*itNode, variables, true);
              }
            else
              {
                continue;
              }

            break;

          case CNodeIteratorMode::After:

            // We do not need to check whether the root is boolean as non boolean nodes are
            // already processed
            switch ((int) itNode->getType())
              {
                case(CEvaluationNode::LOGICAL | CEvaluationNodeLogical::AND):
                case(CEvaluationNode::LOGICAL | CEvaluationNodeLogical::OR):
                case(CEvaluationNode::LOGICAL | CEvaluationNodeLogical::XOR):
                  pNode = compileAND(*itNode, itNode.context(), variables, pRoot, container);
                  break;

                case(CEvaluationNode::LOGICAL | CEvaluationNodeLogical::EQ):
                  pNode = compileEQ(*itNode, itNode.context(), variables, pRoot, container);
                  break;

                case(CEvaluationNode::LOGICAL | CEvaluationNodeLogical::NE):
                  pNode = compileNE(*itNode, itNode.context(), variables, pRoot, container);
                  break;

                case(CEvaluationNode::LOGICAL | CEvaluationNodeLogical::LE):
                case(CEvaluationNode::LOGICAL | CEvaluationNodeLogical::LT):
                case(CEvaluationNode::LOGICAL | CEvaluationNodeLogical::GE):
                case(CEvaluationNode::LOGICAL | CEvaluationNodeLogical::GT):
                  pNode = compileLE(*itNode, itNode.context(), variables, pRoot, container);
                  break;

                case(CEvaluationNode::FUNCTION | CEvaluationNodeFunction::NOT):
                  pNode = compileNOT(*itNode, itNode.context(), variables, pRoot, container);
                  break;

                case(CEvaluationNode::CALL | CEvaluationNodeCall::FUNCTION):
                case(CEvaluationNode::CALL | CEvaluationNodeCall::EXPRESSION):
                  pNode = compileFUNCTION(*itNode, itNode.context(), variables, pRoot, container);
                  break;

                case(CEvaluationNode::VARIABLE | CEvaluationNodeVariable::ANY):
                  pNode = compileVARIABLE(*itNode, itNode.context(), variables, pRoot, container);
                  break;

                case(CEvaluationNode::CONSTANT | CEvaluationNodeConstant::TRUE):
                case(CEvaluationNode::CONSTANT | CEvaluationNodeConstant::FALSE):
                default:
                  pNode = itNode->copyNode(itNode.context());
                  break;
              }

            break;

          default:
            // This will never happen
            break;
        }

      if (itNode.parentContextPtr() != NULL)
        {
          itNode.parentContextPtr()->push_back(pNode);
        }
    }

  return pNode;
}

// static
CEvaluationNode * CMathEventN::CTrigger::compileAND(const CEvaluationNode * pTriggerNode,
    const std::vector< CEvaluationNode * > & children,
    const CMath::Variables< CEvaluationNode * > & /* variables */,
    CMathEventN::CTrigger::CRoot *& /* pRoot */,
    CMathContainer & /* container */)
{
  CEvaluationNode * pNode = NULL;

  switch ((int) CEvaluationNode::subType(pTriggerNode->getType()))
    {
      case CEvaluationNodeLogical::AND:
        pNode = new CEvaluationNodeLogical(CEvaluationNodeLogical::AND, "AND");
        break;

      case CEvaluationNodeLogical::OR:
        pNode = new CEvaluationNodeLogical(CEvaluationNodeLogical::OR, "OR");
        break;

      case CEvaluationNodeLogical::XOR:
        pNode = new CEvaluationNodeLogical(CEvaluationNodeLogical::XOR, "XOR");
        break;

      default:
        break;
    }

  if (pNode != NULL)
    {
      pNode->addChild(children[0]);
      pNode->addChild(children[1]);
    }

  return pNode;
}

// static
CEvaluationNode * CMathEventN::CTrigger::compileEQ(const CEvaluationNode * pTriggerNode,
    const std::vector< CEvaluationNode * > & children,
    const CMath::Variables< CEvaluationNode * > & variables,
    CMathEventN::CTrigger::CRoot *& pRoot,
    CMathContainer & container)
{
  CEvaluationNode * pNode = NULL;

  // Equality can be determined between Boolean and double values.
  if (!static_cast< const CEvaluationNode * >(pTriggerNode->getChild())->isBoolean())
    {
      // We treat x EQ y as (x GE y) AND (y GE x)
      pNode = new CEvaluationNodeLogical(CEvaluationNodeLogical::AND, "AND");

      CEvaluationNodeLogical GELeft(CEvaluationNodeLogical::GE, "GE");
      CEvaluationNode * pGELeft = compileLE(&GELeft, children, variables, pRoot, container);
      pNode->addChild(pGELeft);

      // We need to duplicate and reverse the order for the right
      std::vector< CEvaluationNode * > RightChildren;
      RightChildren.push_back(children[1]->copyBranch());
      RightChildren.push_back(children[0]->copyBranch());

      CEvaluationNodeLogical GERight(CEvaluationNodeLogical::GE, "GE");
      CEvaluationNode * pGERight = compileLE(&GERight, RightChildren, variables, pRoot, container);
      pNode->addChild(pGERight);
    }
  else
    {
      pNode = new CEvaluationNodeLogical(CEvaluationNodeLogical::EQ, "EQ");
      pNode->addChild(children[0]);
      pNode->addChild(children[1]);
    }

  return pNode;
}

// static
CEvaluationNode * CMathEventN::CTrigger::compileNE(const CEvaluationNode * /* pTriggerNode */,
    const std::vector< CEvaluationNode * > & children,
    const CMath::Variables< CEvaluationNode * > & variables,
    CMathEventN::CTrigger::CRoot *& pRoot,
    CMathContainer & container)
{
  CEvaluationNode * pNode = NULL;

  // We treat this as NOT and EQ.
  // For this we create a modified copy of the current node.

  CEvaluationNode * pNotNode = new CEvaluationNodeFunction(CEvaluationNodeFunction::NOT, "NOT");

  CEvaluationNodeLogical EqNode(CEvaluationNodeLogical::EQ, "EQ");

  EqNode.addChild(children[0]);
  EqNode.addChild(children[1]);

  CEvaluationNode * pEqNode = compileEQ(&EqNode, children, variables, pRoot, container);
  pNotNode->addChild(pEqNode);

  // We need to remove the children since the ownership has been transferred to pEqNode.
  EqNode.removeChild(children[0]);
  EqNode.removeChild(children[1]);

  return pNode;
}

// static
CEvaluationNode * CMathEventN::CTrigger::compileLE(const CEvaluationNode * pTriggerNode,
    const std::vector< CEvaluationNode * > & children,
    const CMath::Variables< CEvaluationNode * > & /* variables */,
    CMathEventN::CTrigger::CRoot *& pRoot,
    CMathContainer & container)
{
  CEvaluationNode * pNode = NULL;

  // We need to compile the root finding structure
  // Create a root expression
  CEvaluationNode * pRootNode = new CEvaluationNodeOperator(CEvaluationNodeOperator::MINUS, "-");

  bool Equality = false;

  // We need to create a copy the left and right data nodes with the variables being replaced.
  switch ((int) CEvaluationNode::subType(pTriggerNode->getType()))
    {
      case CEvaluationNodeLogical::LE:
        pRootNode->addChild(children[1]);
        pRootNode->addChild(children[0]);
        Equality = true;
        break;

      case CEvaluationNodeLogical::LT:
        pRootNode->addChild(children[1]);
        pRootNode->addChild(children[0]);
        Equality = false;
        break;

      case CEvaluationNodeLogical::GE:
        pRootNode->addChild(children[0]);
        pRootNode->addChild(children[1]);
        Equality = true;
        break;

      case CEvaluationNodeLogical::GT:
        pRootNode->addChild(children[0]);
        pRootNode->addChild(children[1]);
        Equality = false;
        break;
    }

  pRoot->compile(pRootNode, Equality, container);
  pNode = pRoot->createTriggerExpressionNode();
  pRoot++;

  // We do not need to delete pRootNode as CRoot::compile takes car of it.

  return pNode;
}

// static
CEvaluationNode * CMathEventN::CTrigger::compileNOT(const CEvaluationNode * /* pTriggerNode */,
    const std::vector< CEvaluationNode * > & children,
    const CMath::Variables< CEvaluationNode * > & /* variables */,
    CMathEventN::CTrigger::CRoot *& /* pRoot */,
    CMathContainer & /* container */)
{
  CEvaluationNode * pNode = NULL;

  pNode = new CEvaluationNodeFunction(CEvaluationNodeFunction::NOT, "NOT");
  pNode->addChild(children[0]);

  return pNode;
}

// static
CEvaluationNode * CMathEventN::CTrigger::compileFUNCTION(const CEvaluationNode * pTriggerNode,
    const std::vector< CEvaluationNode * > & children,
    const CMath::Variables< CEvaluationNode * > & /* variables */,
    CMathEventN::CTrigger::CRoot *& pRoot,
    CMathContainer & container)
{
  const CEvaluationNode * pCalledNode =
    static_cast< const CEvaluationNodeCall * >(pTriggerNode)->getCalledTree()->getRoot();

  CEvaluationNode * pNode = compile(pCalledNode, children, pRoot, container);

  // We need to delete the children as the variables have been copied in place.
  std::vector< CEvaluationNode * >::const_iterator it = children.begin();
  std::vector< CEvaluationNode * >::const_iterator end = children.end();

  for (; it != end; ++it)
    {
      delete *it;
    }

  return pNode;
}

// static
CEvaluationNode * CMathEventN::CTrigger::compileVARIABLE(const CEvaluationNode * pTriggerNode,
    const std::vector< CEvaluationNode * > & /* children */,
    const CMath::Variables< CEvaluationNode * > & variables,
    CMathEventN::CTrigger::CRoot *& /* pRoot */,
    CMathContainer & /* container */)
{
  // We need to mimic the process in CMathContainer::copyBranch;
  size_t Index =
    static_cast< const CEvaluationNodeVariable * >(pTriggerNode)->getIndex();

  if (Index != C_INVALID_INDEX &&
      Index < variables.size())
    {
      // Since a variable may be referred to multiple times we need to copy it.
      return variables[Index]->copyBranch();
    }
  else
    {
      // Variables must not appear in mathematical expressions.
      // We create a constant node with the variable name and value NaN.
      return new CEvaluationNodeConstant(CEvaluationNodeConstant::_NaN, pTriggerNode->getData());
    }
}

// static
void CMathEventN::allocate(CMathEventN * pEvent,
                           const CEvent * pDataEvent,
                           const CMathContainer & container)
{
  pEvent->mTrigger.allocate(pDataEvent, container);
  pEvent->mAssignments.resize(pDataEvent->getAssignments().size());
}

CMathEventN::CMathEventN():
  mType(CEvent::Assignment),
  mTrigger(),
  mAssignments(),
  mpDelay(NULL),
  mpPriority(NULL),
  mFireAtInitialTime(false),
  mPersistentTrigger(false),
  mDelayAssignment(true)
{}

/**
 * Destructor
 */
CMathEventN::~CMathEventN()
{}

const CMathEventN::CTrigger & CMathEventN::getTrigger() const
{
  return mTrigger;
}

void CMathEventN::initialize(CMath::sPointers & pointers)
{
  mTrigger.initialize(pointers);

  CAssignment * pAssignment = mAssignments.array();
  CAssignment * pAssignmentEnd = pAssignment + mAssignments.size();

  // Initialize assignments.
  for (; pAssignment != pAssignmentEnd; ++pAssignment)
    {
      pAssignment->initialize(pointers);
    }

  // Initialize delay object.
  mpDelay = pointers.pEventDelaysObject;
  CMathObject::initialize(pointers.pEventDelaysObject, pointers.pEventDelays,
                          CMath::EventDelay, CMath::Event, CMath::SimulationTypeUndefined,
                          false, false, NULL);

  // Initialize priority object.
  mpPriority = pointers.pEventDelaysObject;
  CMathObject::initialize(pointers.pEventPrioritiesObject, pointers.pEventPriorities,
                          CMath::EventPriority, CMath::Event, CMath::SimulationTypeUndefined,
                          false, false, NULL);
}

void CMathEventN::copy(const CMathEventN & src, CMathContainer & container, const size_t & valueOffset, const size_t & objectOffset)
{
  mType = src.mType;
  mTrigger.copy(src.mTrigger, container, valueOffset, objectOffset);

  mAssignments.resize(src.mAssignments.size());
  CAssignment * pAssignment = mAssignments.array();
  CAssignment * pAssignmentEnd = pAssignment + mAssignments.size();
  const CAssignment * pAssignmentSrc = src.mAssignments.array();

  for (; pAssignment != pAssignmentEnd; ++pAssignment, ++pAssignmentSrc)
    {
      pAssignment->copy(*pAssignmentSrc, container, valueOffset, objectOffset);
    }

  mpDelay = src.mpDelay + objectOffset;
  mpPriority = src.mpPriority + objectOffset;
  mFireAtInitialTime = src.mFireAtInitialTime;
  mPersistentTrigger = src.mPersistentTrigger;
  mDelayAssignment = src.mDelayAssignment;
}

bool CMathEventN::compile(CEvent * pDataEvent,
                          CMathContainer & container)
{
  bool success = true;

  mFireAtInitialTime = pDataEvent->getFireAtInitialTime();
  mPersistentTrigger = pDataEvent->getPersistentTrigger();
  mDelayAssignment = pDataEvent->getDelayAssignment();

  // Compile Trigger
  success &= mTrigger.compile(pDataEvent, container);

  CAssignment * pAssignment = mAssignments.array();
  CAssignment * pAssignmentEnd = pAssignment + mAssignments.size();
  CCopasiVector< CEventAssignment >::const_iterator itAssignment = pDataEvent->getAssignments().begin();

  // Compile assignments.
  for (; pAssignment != pAssignmentEnd; ++pAssignment, ++itAssignment)
    {
      success &= pAssignment->compile(*itAssignment, container);
    }

  std::vector< CCopasiContainer * > ListOfContainer;
  // ListOfContainer.push_back(const_cast< CMathContainer * >(&container));

  // Compile the delay object.
  CExpression DelayExpression("DelayExpression", &container);
  success &= DelayExpression.setInfix(pDataEvent->getDelayExpression());
  success &= DelayExpression.compile(ListOfContainer);
  success &= mpDelay->setExpression(DelayExpression, container);

  // Compile the priority object.
  CExpression PriorityExpression("PriorityExpression", &container);
  success &= PriorityExpression.setInfix(pDataEvent->getPriorityExpression());
  success &= PriorityExpression.compile(ListOfContainer);
  success &= mpPriority->setExpression(PriorityExpression, container);

  return success;
}

bool CMathEventN::compile(CMathContainer & container)
{
  bool success = true;

  // Compile Trigger
  success &= mTrigger.compile(NULL, container);

  // Compile assignments.
  // Nothing to do since the target and expression objects are already compiled

  // The delay object is already compiled.
  // The priority object is already compiled.

  return success;
}

void CMathEventN::setTriggerExpression(const std::string & infix, CMathContainer & container)
{
  mTrigger.setExpression(infix, container);
}

void CMathEventN::setDelayExpression(const std::string & infix, CMathContainer & container)
{
  assert(mpDelay != NULL);

  mpDelay->setExpression(infix, false, container);
}

void CMathEventN::setPriorityExpression(const std::string & infix, CMathContainer & container)
{
  assert(mpPriority != NULL);

  mpPriority->setExpression(infix, false, container);
}

void CMathEventN::addAssignment(CMathObject * pTarget, CMathObject * pExpression)
{
  size_t OldSize = mAssignments.size();
  mAssignments.resize(OldSize + 1, true);

  CAssignment & Assignment = mAssignments[OldSize];
  Assignment.setTarget(pTarget);
  Assignment.setExpression(pExpression);
}

const CVector< CMathEventN::CAssignment > & CMathEventN::getAssignments() const
{
  return mAssignments;
}

CMathEvent::CAssignment::CAssignment(const CCopasiContainer * pParent) :
  CCopasiContainer("MathEventAssignment", pParent),
  mpTarget(NULL),
  mExpression("Expression", this)
{}

CMathEvent::CAssignment::CAssignment(const CMathEvent::CAssignment & src,
                                     const CCopasiContainer * pParent) :
  CCopasiContainer(src, pParent),
  mpTarget(src.mpTarget),
  mExpression(src.mExpression, this)
{}

CMathEvent::CAssignment::~CAssignment()
{}

bool CMathEvent::CAssignment::compile(const CEventAssignment * pAssignment,
                                      std::vector< CCopasiContainer * > listOfContainer)
{
  if (pAssignment == NULL)
    return false;

  bool success = true;

  mpTarget = NULL;

  if (pAssignment->getTargetObject() != NULL)
    {
      mpTarget = (C_FLOAT64 *) pAssignment->getTargetObject()->getValuePointer();
    }

  if (mpTarget == NULL)
    success = false;

  success &= mExpression.setInfix(pAssignment->getExpression());
  success &= mExpression.compile(listOfContainer);

  setDirectDependencies(mExpression.getDirectDependencies());

  return success;
}

CMathEvent::CMathEvent(const CCopasiContainer * pParent) :
  CCopasiContainer("MathEvent", pParent, "MathEvent"),
  mTrigger(this),
  mOrder(false),
  mHaveDelay(false),
  mDelay("DelayExpression", this),
  mDelayAssignment(true),
  mAssignments("ListOfMathEventAssignment", this),
  mDelayValueRefreshes(),
  mAssignmentValueRefreshes(),
  mDependentValueRefreshes(),
  mType(CEvent::Assignment)
{}

CMathEvent::CMathEvent(const CMathEvent & src,
                       const CCopasiContainer * pParent) :
  CCopasiContainer(src, pParent),
  mTrigger(src.mTrigger, this),
  mOrder(src.mOrder),
  mHaveDelay(src.mHaveDelay),
  mDelay(src.mDelay, this),
  mDelayAssignment(src.mDelayAssignment),
  mAssignments(src.mAssignments, this),
  mDelayValueRefreshes(src.mDelayValueRefreshes),
  mAssignmentValueRefreshes(src.mAssignmentValueRefreshes),
  mDependentValueRefreshes(src.mDependentValueRefreshes),
  mType(src.mType)
{}

CMathEvent::~CMathEvent()
{}

bool CMathEvent::compile(const CEvent * pEvent,
                         std::vector< CCopasiContainer * > listOfContainer)
{
  // A CMathEvent must be part of CMathModel to be compiled.
  CMathModel * pMathModel = dynamic_cast< CMathModel *>(getObjectAncestor("CMathModel"));

  if (pMathModel == NULL)
    return false;

  bool success = true;

  success &= mTrigger.compile(pEvent->getTriggerExpressionPtr(), listOfContainer);

  success &= mDelay.setInfix(pEvent->getDelayExpression());
  success &= mDelay.compile(listOfContainer);

  mHaveDelay = (mDelay.getInfix() != "");

  mType = pEvent->getType();

  // Build the list of refresh calls needed to assure that the delay expression
  // can be calculated.

  mDelayValueRefreshes = pMathModel->buildRequiredRefreshList(mDelay.getDirectDependencies());

  mDelayAssignment = pEvent->getDelayAssignment();

  mAssignments.clear();

  CCopasiVectorN< CEventAssignment >::const_iterator it = pEvent->getAssignments().begin();
  CCopasiVectorN< CEventAssignment >::const_iterator end = pEvent->getAssignments().end();

  std::set< const CCopasiObject * > Assignments;
  std::set< const CCopasiObject * > Targets;

  for (; it != end; ++it)
    {
      CAssignment * pAssignment = new CAssignment();

      mAssignments.add(pAssignment, true);
      success &= pAssignment->compile(*it, listOfContainer);

      Assignments.insert(pAssignment);
      Targets.insert((*it)->getTargetObject());
    }

  // Build the list of refresh calls needed to assure that the assignment expressions
  // can be calculated.
  mAssignmentValueRefreshes = pMathModel->buildRequiredRefreshList(Assignments);

  // Build the list of refresh calls needed to assure that all dependent model values
  // are updated after the assignments are executed.
  mDependentValueRefreshes = pMathModel->buildDependendRefreshList(Targets);

  return success;
}

void CMathEvent::fire(const C_FLOAT64 & time,
                      const bool & equality,
                      CProcessQueue & processQueue)
{
  if (mDelayAssignment)
    {
      processQueue.addAssignment(getAssignmentTime(time), equality, getTargetValues(), this);
    }
  else
    {
      processQueue.addCalculation(getCalculationTime(time), equality, this);
    }
}

CVector< C_FLOAT64 > CMathEvent::getTargetValues()
{
  applyValueRefreshes();

  CVector< C_FLOAT64 > Values(mAssignments.size());
  C_FLOAT64 * pValue = Values.array();
  CCopasiVector< CAssignment >::iterator itAssignment = mAssignments.begin();
  CCopasiVector< CAssignment >::iterator endAssignment = mAssignments.end();

  for (; itAssignment != endAssignment; ++itAssignment, ++pValue)
    {
      *pValue = (*itAssignment)->mExpression.calcValue();
    }

  return Values;
}

bool CMathEvent::setTargetValues(const CVector< C_FLOAT64 > & values)
{
  bool StateChanged = false;

  const C_FLOAT64 * pValue = values.array();
  CCopasiVector< CAssignment >::iterator itAssignment = mAssignments.begin();
  CCopasiVector< CAssignment >::iterator endAssignment = mAssignments.end();

  for (; itAssignment != endAssignment; ++itAssignment, ++pValue)
    {
      if (*(*itAssignment)->mpTarget != *pValue)
        {
          StateChanged = true;
          *(*itAssignment)->mpTarget = *pValue;
        }
    }

  if (StateChanged)
    {
      applyDependentRefreshes();
    }

  return StateChanged;
}

bool CMathEvent::executeAssignment()
{
  return setTargetValues(getTargetValues());
}

void CMathEvent::applyDelayRefreshes()
{
  std::vector< Refresh * >::const_iterator itRefresh = mDelayValueRefreshes.begin();
  std::vector< Refresh * >::const_iterator endRefresh = mDelayValueRefreshes.end();

  while (itRefresh != endRefresh)
    (**itRefresh++)();

  return;
}

void CMathEvent::applyValueRefreshes()
{
  std::vector< Refresh * >::const_iterator itRefresh = mAssignmentValueRefreshes.begin();
  std::vector< Refresh * >::const_iterator endRefresh = mAssignmentValueRefreshes.end();

  while (itRefresh != endRefresh)
    (**itRefresh++)();

  return;
}

void CMathEvent::applyDependentRefreshes()
{
  std::vector< Refresh * >::const_iterator itRefresh = mDependentValueRefreshes.begin();
  std::vector< Refresh * >::const_iterator endRefresh = mDependentValueRefreshes.end();

  while (itRefresh != endRefresh)
    (**itRefresh++)();

  return;
}

CMathTrigger & CMathEvent::getMathTrigger()
{
  return mTrigger;
}

const size_t & CMathEvent::getOrder() const
{
  return mOrder;
}

const bool & CMathEvent::delayAssignment() const
{
  return mDelayAssignment;
}

C_FLOAT64 CMathEvent::getCalculationTime(const C_FLOAT64 & currentTime)
{
  if (mDelayAssignment)
    {
      return currentTime;
    }

  return calculateDelayedTime(currentTime);
}

C_FLOAT64 CMathEvent::getAssignmentTime(const C_FLOAT64 & currentTime)
{
  if (!mDelayAssignment)
    {
      return currentTime;
    }

  return calculateDelayedTime(currentTime);
}

C_FLOAT64 CMathEvent::calculateDelayedTime(const C_FLOAT64 & currentTime)
{
  if (mDelay.getInfix() == "")
    {
      return currentTime;
    }

  // We make sure everything is up to date.
  applyDelayRefreshes();

  C_FLOAT64 DelayedTime = currentTime + mDelay.calcValue();

  // Events are only allowed in forward integration. Thus the ExecutionTime
  // must not be less than the time.
  if (DelayedTime - currentTime < 0.0)
    {
      // We allow small numerical errors.
      C_FLOAT64 Scale =
        (fabs(DelayedTime) + fabs(currentTime)) * 50.0 * std::numeric_limits< C_FLOAT64 >::epsilon();

      // Both are approximately zero
      if (Scale < 100.0 * std::numeric_limits< C_FLOAT64 >::min())
        {
          DelayedTime = currentTime;
        }
      // The difference is small compared to the scale
      else if (fabs(DelayedTime - currentTime) < Scale)
        {
          DelayedTime = currentTime;
        }
      // The execution time is definitely in the past
      else
        {
          // Create an error message and throw an exception.
          CCopasiMessage(CCopasiMessage::EXCEPTION, MCMathModel + 2, DelayedTime, currentTime);
        }
    }

  return DelayedTime;
}

const CEvent::Type & CMathEvent::getType() const
{
  return mType;
}
