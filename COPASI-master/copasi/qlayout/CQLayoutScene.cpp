// Copyright (C) 2013 by Pedro Mendes, Virginia Tech Intellectual
// Properties, Inc., University of Heidelberg, and The University
// of Manchester.
// All rights reserved.

#include <QtCore/QCoreApplication>
#include <QtGui/QGraphicsScene>
#include <QtGui/QGraphicsItem>
#include <QtGui/QPainter>
#include <QtGui/QPrinter>
#include <QtGui/QImage>
#include <QtGui/QGraphicsEffect>

#include "copasi.h"

#include "qlayout/CQLayoutScene.h"
#include "qlayout/CQCopasiEffect.h"
#include "qlayout/CQLabelGraphicsItem.h"
#include "qlayout/CQStyledGraphicsItem.h"
#include "qlayout/CQConnectionGraphicsItem.h"
#include "qlayout/CQRenderConverter.h"
#include "layout/CLayout.h"
#include "layout/CLGlyphs.h"
#include "layout/CLText.h"
#include "layout/CLReactionGlyph.h"
#include "layout/CLRenderResolver.h"
#include "layout/CLGlobalRenderInformation.h"
#include "layout/CListOfLayouts.h"
#include "layout/CLLocalRenderInformation.h"
#include "layout/CLDefaultStyles.h"
#include "layout/CCopasiSpringLayout.h"

#include "CopasiDataModel/CCopasiDataModel.h"
#include "copasi/report/CCopasiRootContainer.h"
#include "copasi/report/CKeyFactory.h"

CQLayoutScene::CQLayoutScene(CLayout* layout, CCopasiDataModel* model, CLRenderInformationBase* renderInformation)
  : QGraphicsScene()
  , mpLayout(layout)
  , mpRender(renderInformation)
  , mpResolver(NULL)
{
  initializeResolver(model, renderInformation);
  connect(this, SIGNAL(recreateNeeded()), this, SLOT(recreate()), Qt::QueuedConnection);
}

void CQLayoutScene::setLayout(CLayout *layout, CCopasiDataModel* model, CLRenderInformationBase* renderInformation)
{
  mpLayout = layout;
  setRenderInformation(model, renderInformation);
}

void CQLayoutScene::setRenderInformation(CCopasiDataModel* model, CLRenderInformationBase* renderInformation)
{
  initializeResolver(model, renderInformation);
}

const CLayout* CQLayoutScene::getCurrentLayout() const
{
  return mpLayout;
}

CLayout* CQLayoutScene::getCurrentLayout()
{
  return mpLayout;
}

const CLRenderInformationBase* CQLayoutScene::getCurrentRenderInfo() const
{
  return mpRender;
}

void CQLayoutScene::saveToFile(const std::string& fileName, const std::string& fileType /*= "pdf"*/)
{
  if (fileType == "pdf")
    {
      QPrinter printer(QPrinter::HighResolution);
      printer.setOutputFormat(QPrinter::PdfFormat);
      printer.setOutputFileName(fileName.c_str());
      QPainter painter(&printer);
      painter.setRenderHints(
        QPainter::Antialiasing | QPainter::HighQualityAntialiasing | QPainter::SmoothPixmapTransform);
      render(&painter, QRect(), itemsBoundingRect());
      painter.end();
    }
  else
    {
      const int scale = 2;
      QImage image(QSize(width()*scale, height()*scale), QImage::Format_ARGB32);
      QPainter painter(&image);
      painter.setRenderHints(
        QPainter::Antialiasing | QPainter::HighQualityAntialiasing | QPainter::SmoothPixmapTransform);
      render(&painter, image.rect(), itemsBoundingRect());
      painter.end();
      image.save(fileName.c_str(), fileType.c_str());
    }
}

void CQLayoutScene::initializeResolver(CCopasiDataModel* model, CLRenderInformationBase* renderInformation)
{
  if (model == NULL)
    return;

  if (renderInformation == NULL)
    {
      if (mpLayout != NULL && mpLayout->getListOfLocalRenderInformationObjects().size() > 0)
        mpRender = mpLayout->getListOfLocalRenderInformationObjects()[0];
      else if (model->getListOfLayouts()->getListOfGlobalRenderInformationObjects().size() > 0)
        mpRender = model->getListOfLayouts()->getListOfGlobalRenderInformationObjects()[0];
      else
        mpRender = getDefaultStyle(0);
    }
  else
    mpRender = renderInformation;

  if (mpLayout == NULL || mpRender == NULL)
    return;

  CLLocalRenderInformation* local = dynamic_cast<CLLocalRenderInformation*>(mpRender);

  if (local != NULL)
    mpResolver = QSharedPointer<CLRenderResolver>(new CLRenderResolver(*local, mpLayout->getListOfLocalRenderInformationObjects(),   model->getListOfLayouts()->getListOfGlobalRenderInformationObjects()));
  else
    mpResolver = QSharedPointer<CLRenderResolver>(new CLRenderResolver(*dynamic_cast<CLGlobalRenderInformation*>(mpRender), model->getListOfLayouts()->getListOfGlobalRenderInformationObjects()));
}

void CQLayoutScene::setResolver(CLRenderResolver* resolver)
{
  mpResolver = QSharedPointer<CLRenderResolver>(resolver);
}

const CLRenderResolver* CQLayoutScene::getResolver() const
{
  return mpResolver.data();
}

CQLayoutScene::~CQLayoutScene()
{
}

void CQLayoutScene::recreate()
{
  fillFromLayout(mpLayout);
  invalidate();
}

void CQLayoutScene::addGlyph(const CLGraphicalObject* go)
{
  if (go == NULL) return;

  const CLGlyphWithCurve* curveGlyph = dynamic_cast<const CLGlyphWithCurve*>(go);
  const CLReactionGlyph* reaction = dynamic_cast<const CLReactionGlyph*>(go);
  const CLTextGlyph* text = dynamic_cast<const CLTextGlyph*>(go);
  const CLGeneralGlyph* general = dynamic_cast<const CLGeneralGlyph*>(go);
  QGraphicsItem *item = NULL;

  if (curveGlyph != NULL)
    {
      if (curveGlyph->getCurve().getNumCurveSegments() > 0 || reaction != NULL || general != NULL)
        item = new CQConnectionGraphicsItem(curveGlyph,
                                            mpResolver == NULL ? NULL : mpResolver.data());
    }
  else if (text != NULL)
    {
      item = new CQLabelGraphicsItem(text, mpResolver == NULL ? NULL : mpResolver.data());
    }
  else
    {
      item = new CQStyledGraphicsItem(go, mpResolver == NULL ? NULL : mpResolver.data());
    }

  if (item != NULL)
    {
      CCopasiObject* obj = go->getModelObject();

      if (obj != NULL && text == NULL)
        {
          item->setData(COPASI_OBJECT_CN, QString(obj->getCN().c_str()));
          mItems[obj->getCN()] = item;
        }

      addItem(item);
    }

  if (general != NULL)
    {
      const CCopasiVector<CLGraphicalObject> & subGlyphs = general->getListOfSubglyphs();
      CCopasiVector<CLGraphicalObject>::const_iterator it = subGlyphs.begin();

      while (it != subGlyphs.end())
        {
          addGlyph(*it);
          ++it;
        }
    }
}

QGraphicsItem* CQLayoutScene::getItemFor(const std::string& cn)
{
  return mItems[cn];
}

void CQLayoutScene::fillFromLayout(const CLayout* layout)
{
  if (layout == NULL) return;

  clear();
  mItems.clear();

  if (mpRender != NULL && mpResolver != NULL)
    {
      CQRenderConverter::setBackground(this, mpRender->getBackgroundColor(), mpResolver.data());
    }

  const CCopasiVector<CLCompartmentGlyph> & comps = layout->getListOfCompartmentGlyphs();

  CCopasiVector<CLCompartmentGlyph>::const_iterator itComp = comps.begin();

  while (itComp != comps.end())
    {
      addGlyph(*itComp);
      ++itComp;
    }

  const CCopasiVector<CLReactionGlyph> & reactions = layout->getListOfReactionGlyphs();

  CCopasiVector<CLReactionGlyph>::const_iterator itReactions = reactions.begin();

  while (itReactions != reactions.end())
    {
      addGlyph(*itReactions);
      ++itReactions;
    }

  const CCopasiVector<CLMetabGlyph> & species = layout->getListOfMetaboliteGlyphs();

  CCopasiVector<CLMetabGlyph>::const_iterator itSpecies = species.begin();

  while (itSpecies != species.end())
    {
      addGlyph(*itSpecies);
      ++itSpecies;
    }

  const CCopasiVector<CLTextGlyph> & texts = layout->getListOfTextGlyphs();

  CCopasiVector<CLTextGlyph>::const_iterator itTexts = texts.begin();

  while (itTexts != texts.end())
    {
      addGlyph(*itTexts);
      ++itTexts;
    }

  const CCopasiVector<CLGeneralGlyph> & list = layout->getListOfGeneralGlyphs();

  CCopasiVector<CLGeneralGlyph>::const_iterator itList = list.begin();

  while (itList != list.end())
    {
      addGlyph(*itList);
      ++itList;
    }
}

CLGraphicalObject* getTextForItem(const CLayout* layout, const CLGraphicalObject* obj)
{
  const CCopasiVector<CLTextGlyph> & texts = layout->getListOfTextGlyphs();
  CCopasiVector<CLTextGlyph>::const_iterator it = texts.begin();

  while (it != texts.end())
    {
      if ((*it)->getGraphicalObjectKey() == obj->getKey())
        return *it;

      ++it;
    }

  return NULL;
}

#include <model/CCompartment.h>
#include <model/CMetab.h>

CLGraphicalObject* getMetabGlyphForKey(const CLayout* layout, const CMetab* metab)
{
  const CCopasiVector<CLMetabGlyph> & metabs = layout->getListOfMetaboliteGlyphs();
  CCopasiVector<CLMetabGlyph>::const_iterator it = metabs.begin();

  while (it != metabs.end())
    {
      if ((*it)->getModelObjectKey() == metab->getKey())
        return *it;

      ++it;
    }

  return NULL;
}

void moveObject(CLGraphicalObject* obj, const CLPoint& delta, CLayout* layout)
{
  if (obj == NULL) return;

  // move object
  obj->moveBy(delta);

  // move its label
  CLGraphicalObject* text = getTextForItem(layout, obj);

  if (text != NULL)
    text->moveBy(delta);

  // move species within compartments as well
  CLCompartmentGlyph* lcomp = dynamic_cast<CLCompartmentGlyph*>(obj);

  if (lcomp == NULL)
    return;

  CCompartment*  comp = dynamic_cast<CCompartment*>(lcomp ->getModelObject());

  if (comp == NULL)
    return;

  CCopasiVectorNS < CMetab > & metabs = comp->getMetabolites();
  CCopasiVectorNS < CMetab >::const_iterator it = metabs.begin();

  while (it != metabs.end())
    {
      moveObject(getMetabGlyphForKey(layout, (*it)), delta, layout);
      ++it;
    }
}

void CQLayoutScene::updatePosition(const QString& key, const QPointF& newPos)
{
  CKeyFactory* kf = CCopasiRootContainer::getKeyFactory();

  if (kf == NULL) return;

  CLGraphicalObject* obj = dynamic_cast<CLGraphicalObject*>(kf->get(key.toStdString()));

  if (obj == NULL) return;

  CLPoint delta(newPos.x(), newPos.y());
  moveObject(obj, delta, mpLayout);

  // restore lines
  CCopasiSpringLayout::Parameters p;
  CCopasiSpringLayout l(mpLayout, &p);
  l.finalizeState();

  emit recreateNeeded();
}
