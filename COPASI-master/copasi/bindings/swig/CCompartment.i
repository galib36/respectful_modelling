// Begin CVS Header 
//   $Source: /fs/turing/cvs/copasi_dev/copasi/bindings/swig/CCompartment.i,v $ 
//   $Revision: 1.8 $ 
//   $Name:  $ 
//   $Author: bergmann $ 
//   $Date: 2012/04/11 15:40:24 $ 
// End CVS Header 

// Copyright (C) 2012 - 2010 by Pedro Mendes, Virginia Tech Intellectual 
// Properties, Inc., University of Heidelberg, and The University 
// of Manchester. 
// All rights reserved. 

// Copyright (C) 2008 by Pedro Mendes, Virginia Tech Intellectual 
// Properties, Inc., EML Research, gGmbH, University of Heidelberg, 
// and The University of Manchester. 
// All rights reserved. 

// Copyright (C) 2001 - 2007 by Pedro Mendes, Virginia Tech Intellectual 
// Properties, Inc. and EML Research, gGmbH. 
// All rights reserved. 

%{

#include "model/CCompartment.h"  
#include <assert.h>  
%}

%ignore CCompartment::load;
%ignore CCompartment::getDeletedObjects;

#if (defined SWIGJAVA || defined SWIGCSHARP)
// remove some const methods to get rid of warnings
%ignore CCompartment::getMetabolites() const;

#endif // SWIGJAVA || CSHARP

%include "model/CCompartment.h"

%extend CCompartment{

    /**
     * This method is there for backwards compatibility and will be
     * removed eventually.
     */
    bool removeMetabolite(CMetab* metab)
    {
       CModel* pModel=dynamic_cast<CModel*>($self->getObjectParent()->getObjectParent());
       assert(pModel!=NULL);
       return pModel->removeMetabolite(metab->getKey());
    }
}
