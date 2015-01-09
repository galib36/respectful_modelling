// Copyright (C) 2010 - 2013 by Pedro Mendes, Virginia Tech Intellectual
// Properties, Inc., University of Heidelberg, and The University
// of Manchester.
// All rights reserved.

// Copyright (C) 2009 by Pedro Mendes, Virginia Tech Intellectual
// Properties, Inc., EML Research, gGmbH, University of Heidelberg,
// and The University of Manchester.
// All rights reserved.

#include <QtCore/QVariant>
#include <QtGui/QCheckBox>
#include <QtGui/QLabel>
#include <QtGui/QLineEdit>
#include <QtGui/QPushButton>
#include <QtGui/QRadioButton>
#include <QtGui/QLayout>
#include <QtGui/QToolTip>
#include <QtGui/QWhatsThis>
#include <QtGui/QMessageBox>
#include <QtGui/QToolButton>
#include <QtGui/QImage>

#include <QtGui/QFrame>
#include <QtGui/QFileDialog>
#include <QtGui/QTableWidget>

#include <algorithm>

#include "SensWidgetComboBox.h"
#include "qtUtilities.h"

SensWidgetComboBox::SensWidgetComboBox(QWidget * parent, const char * name)
  : QComboBox(parent)
{
  setObjectName(QString::fromUtf8(name));
}

void SensWidgetComboBox::fillFromList(const std::vector<CObjectLists::ListType> & list)
{
  //store old selection
  CObjectLists::ListType oldItem = getCurrentObjectList();

  mIndexTable = list;

  //fill combobox
  clear();
  std::vector<CObjectLists::ListType>::const_iterator it, itEnd = mIndexTable.end();

  for (it = mIndexTable.begin(); it != itEnd; ++it)
    insertItem(this->count(), FROM_UTF8(CObjectLists::ListTypeName[*it]));

  //restore old selection, if possible
  if (!setCurrentObjectList(oldItem))
    setCurrentIndex(0);
}

CObjectLists::ListType SensWidgetComboBox::getCurrentObjectList() const
{
  unsigned int index = currentIndex();

  if (index < mIndexTable.size())
    return mIndexTable[currentIndex()];
  else
    return CObjectLists::EMPTY_LIST;
}

bool SensWidgetComboBox::setCurrentObjectList(CObjectLists::ListType lt)
{
  std::vector<CObjectLists::ListType>::const_iterator it;
  it = std::find(mIndexTable.begin(),
                 mIndexTable.end(),
                 lt);

  if (it == mIndexTable.end()) return false;

  setCurrentIndex(it - mIndexTable.begin());
  return true;
}
