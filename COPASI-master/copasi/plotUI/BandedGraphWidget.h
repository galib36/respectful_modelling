// Copyright (C) 2011 - 2013 by Pedro Mendes, Virginia Tech Intellectual
// Properties, Inc., University of Heidelberg, and The University
// of Manchester.
// All rights reserved.

#ifndef BANDEDGRAPHWIDGET_H
#define BANDEDGRAPHWIDGET_H

#include <QtCore/QVariant>

#include "copasi/plotUI/ui_BandedGraphWidget.h"
#include "copasi/plotUI/CQPlotEditWidget.h"

class CCopasiObject;
class CModel;
class CPlotItem;

class BandedGraphWidget : public CQPlotEditWidget, public Ui::BandedGraphWidget
{
  Q_OBJECT

public:
  BandedGraphWidget(QWidget* parent = 0, const char* name = 0, Qt::WindowFlags fl = 0);
  ~BandedGraphWidget();

  virtual bool LoadFromCurveSpec(const CPlotItem * curve);
  virtual bool SaveToCurveSpec(CPlotItem * curve, const CPlotItem *original = NULL) const;

  /**
   * In multiple edit mode, we don't want to edit name & channels
   */
  virtual void setMultipleEditMode(bool mode);

protected:
  const CCopasiObject* mpObjectYone;
  const CCopasiObject* mpObjectYtwo;
  const CCopasiObject* mpObjectX;

protected slots:

  virtual void buttonPressedX();
  virtual void buttonPressedY() {buttonPressedYone();};
  void buttonPressedYone();
  void buttonPressedYtwo();
};

#endif // BANDEDGRAPHWIDGET_H
