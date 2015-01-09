# -*- coding: utf-8 -*-
# Begin CVS Header 
#   $Source: /Volumes/Home/Users/shoops/cvs/copasi_dev/copasi/bindings/python/unittests/Test_CCopasiTask.py,v $ 
#   $Revision: 1.11 $ 
#   $Name:  $ 
#   $Author: shoops $ 
#   $Date: 2010/07/16 18:55:59 $ 
# End CVS Header 

# Copyright (C) 2010 by Pedro Mendes, Virginia Tech Intellectual 
# Properties, Inc., University of Heidelberg, and The University 
# of Manchester. 
# All rights reserved. 

# Copyright (C) 2008 by Pedro Mendes, Virginia Tech Intellectual 
# Properties, Inc., EML Research, gGmbH, University of Heidelberg, 
# and The University of Manchester. 
# All rights reserved. 

import COPASI
import unittest
from types import *


class Test_CCopasiTask(unittest.TestCase):
  def setUp(self):
    self.datamodel=COPASI.CCopasiRootContainer.addDatamodel()
    self.task=self.datamodel.getTask(0)


  def test_getType(self):
    t=self.task.getType()
    self.assert_(type(t)==IntType)

  def test_setType(self):
    task=COPASI.CCopasiTask("TestTask")
    task.setType(COPASI.CCopasiTask.scan)
    self.assert_(task.getType()==COPASI.CCopasiTask.scan)

  def test_getKey(self):
    key=self.task.getKey()
    self.assert_(type(key)==StringType)

  def test_setScheduled(self):
    v=False
    self.task.setScheduled(v)
    self.assert_(self.task.isScheduled()==v)
    v=True
    self.task.setScheduled(v)
    self.assert_(self.task.isScheduled()==v)

  def test_isScheduled(self):
    v=self.task.isScheduled()
    self.assert_(type(v)==BooleanType)

  def test_setUpdateModel(self):
    v=False
    self.task.setUpdateModel(v)
    self.assert_(self.task.isUpdateModel()==v)
    v=True
    self.task.setUpdateModel(v)
    self.assert_(self.task.isUpdateModel()==v)

  def test_isUpdateModel(self):
    v=self.task.isUpdateModel()
    self.assert_(type(v)==BooleanType)

  def test_getValidMethods(self):
    task=COPASI.CCopasiTask();
    self.assert_(task!=None)
    validMethods=task.getValidMethods();
    self.assert_(type(validMethods)==TupleType);
    self.assert_(len(validMethods)==0)

  def test_getReport(self):
    report=self.task.getReport()
    self.assert_(report!=None)
    self.assert_(report.__class__==COPASI.CReport)


def suite():
  tests=[
          'test_getType'
         ,'test_setType'
         ,'test_getKey'
         ,'test_isScheduled'
         ,'test_setScheduled'
         ,'test_isUpdateModel'
         ,'test_setUpdateModel'
         ,'test_getValidMethods'
         ,'test_getReport'                       
        ]
  return unittest.TestSuite(map(Test_CCopasiTask,tests))


if(__name__ == '__main__'):
    unittest.TextTestRunner(verbosity=2).run(suite())


                                      
