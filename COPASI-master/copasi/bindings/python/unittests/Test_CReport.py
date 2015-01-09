# -*- coding: utf-8 -*-
# Begin CVS Header 
#   $Source: /Volumes/Home/Users/shoops/cvs/copasi_dev/copasi/bindings/python/unittests/Test_CReport.py,v $ 
#   $Revision: 1.4 $ 
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


class Test_CReport(unittest.TestCase):
  def setUp(self):
    self.datamodel=COPASI.CCopasiRootContainer.addDatamodel()
    self.task=self.datamodel.getTask(0)
    self.report=self.task.getReport()

  def test_getReportDefinition(self):
    reportDefinition=self.report.getReportDefinition()
    self.assert_(reportDefinition!=None)
    self.assert_(reportDefinition.__class__==COPASI.CReportDefinition)

  def test_setReportDefinition(self):
    listOfReportDefinitions=self.datamodel.getReportDefinitionList()
    reportDefinition=listOfReportDefinitions.createReportDefinition("MyReportDefinition","No Comment")
    self.assert_(reportDefinition!=None)
    self.report.setReportDefinition(reportDefinition)
    self.assert_(self.report.getReportDefinition().getKey()==reportDefinition.getKey())

  def test_getTarget(self):
    target=self.report.getTarget()
    self.assert_(target!=None)
    self.assert_(type(target)==StringType)

  def test_setTarget(self):
    target="MyTaget.txt"
    self.report.setTarget(target)
    t=self.report.getTarget()
    self.assert_(t!=None)
    self.assert_(t==target)

  def test_append(self):
    append=self.report.append()
    self.assert_(type(append)==BooleanType)

  def test_setAppend(self):
    self.report.setAppend(True)
    append=self.report.append()
    self.assert_(append==True)
    self.report.setAppend(False)
    append=self.report.append()
    self.assert_(append==False)




def suite():
  tests=[
          'test_getReportDefinition'                       
         ,'test_setReportDefinition'                       
         ,'test_getTarget'                       
         ,'test_setTarget'                       
         ,'test_append'                       
         ,'test_setAppend'                       
        ]
  return unittest.TestSuite(map(Test_CReport,tests))


if(__name__ == '__main__'):
    unittest.TextTestRunner(verbosity=2).run(suite())


                                      
