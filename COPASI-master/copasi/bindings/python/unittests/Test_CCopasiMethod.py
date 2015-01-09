# -*- coding: utf-8 -*-
# Begin CVS Header 
#   $Source: /Volumes/Home/Users/shoops/cvs/copasi_dev/copasi/bindings/python/unittests/Test_CCopasiMethod.py,v $ 
#   $Revision: 1.9 $ 
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


class Test_CCopasiMethod(unittest.TestCase):
  def setUp(self):
    self.datamodel=COPASI.CCopasiRootContainer.addDatamodel()
    self.method=self.datamodel.getTask(0).getMethod()

  def test_TypeNameToEnum(self):
    typename="Enhanced Newton"
    e=COPASI.CCopasiMethod.TypeNameToEnum(typename)
    self.assert_(type(e)==IntType)
    self.assert_(e==COPASI.CCopasiMethod.Newton)

  def test_getType(self):
    t=self.method.getType()
    self.assert_(type(t)==IntType)

  def test_getSubType(self):
    t=self.method.getSubType()
    self.assert_(type(t)==IntType)


def suite():
  tests=[
          'test_TypeNameToEnum'
         ,'test_getType'
         ,'test_getSubType'
        ]
  return unittest.TestSuite(map(Test_CCopasiMethod,tests))

if(__name__ == '__main__'):
    unittest.TextTestRunner(verbosity=2).run(suite())


