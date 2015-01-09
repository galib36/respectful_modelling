// Copyright (C) 2013 by Pedro Mendes, Virginia Tech Intellectual
// Properties, Inc., University of Heidelberg, and The University
// of Manchester.
// All rights reserved.

#ifndef QANIMATION_CONTROLS_H
#define QANIMATION_CONTROLS_H

#include <QtGui/QWidget>

#include <qlayout/ui_CQAnimationControls.h>

class QTimer;

class CQAnimationControls : public QWidget, public Ui::CQAnimationControls
{
  Q_OBJECT
public:
  CQAnimationControls(QWidget* parent = NULL);
  virtual ~CQAnimationControls();

  void setNumSteps(size_t numSteps);
  void setCurrentStep(size_t currentStep);
  bool isPlaying() const;

public slots:
  void slotChangeInterval(int);
  void slotBack();
  void slotForward();
  void slotTogglePlayPause();
  void slotStepBack();
  void slotStepForward();
  void slotStop();
  void slotShowStep(int);

signals:
  void play();
  void pause();
  void stop();
  void forward();
  void backward();
  void step_forward();
  void step_backward();
  void showStep(int);

protected:
  size_t mNumSteps;
  size_t mCurrentStep;
  bool mShouldLoop;
  QTimer *mpTimer;

  void updateButtons();
  void updateActions();
  void createActions();
};

#endif // QANIMATION_CONTROLS_H
