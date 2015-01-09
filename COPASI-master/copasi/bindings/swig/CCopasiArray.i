// Copyright (C) 2010 - 2013 by Pedro Mendes, Virginia Tech Intellectual 
// Properties, Inc., University of Heidelberg, and The University 
// of Manchester. 
// All rights reserved. 

// Copyright (C) 2009 by Pedro Mendes, Virginia Tech Intellectual 
// Properties, Inc., EML Research, gGmbH, University of Heidelberg, 
// and The University of Manchester. 
// All rights reserved. 


%{

#include "report/CCopasiArray.h"

%}

%ignore CCopasiAbstractArray::operator[] (const index_type & index); 
%ignore CCopasiAbstractArray::operator[] (const index_type & index) const;
%ignore CCopasiArray::operator[] (const index_type & index); 
%ignore CCopasiArray::operator[] (const index_type & index) const;
%ignore CCopasiMatrixInterface::operator[] (const index_type & index); 
%ignore CCopasiMatrixInterface::operator[] (const index_type & index) const;
%ignore CCopasiVectorInterface::operator[] (const index_type & index) ;
%ignore CCopasiVectorInterface::operator[] (const index_type & index) const;
%ignore CArrayAnnotation::operator=(const CArrayAnnotation&);
%ignore CArrayAnnotation::array() const;
%ignore operator<<(std::ostream &os, const CArrayAnnotation & o);


#ifdef SWIGPERL
// we need a new typemap for the get function below
// to accept a reference to a cell array as the index argument
%typemap (in) const CCopasiAbstractArray::index_type & MY_FUNNY_INDEX {
  // this typemap is almost identical with the one that
  // is described for handling char** values in section 30.8.1
  // of the SWIG manual  

  // first we check if the argument is a reference
  if (!SvROK($input)) {
    croak("Argument $argnum is not a reference.");
  }
  // the argument has to be an array
  if (SvTYPE(SvRV($input)) != SVt_PVAV) {
    croak("Argument $argnum is not an array.");
  }

  // now we convert the array to a index type which should always be 
  // a std::vector with size_t values
  AV* tmpav = (AV*)SvRV($input);
  // av_len returns the index of the last element, not the length
  I32 len = av_len(tmpav);

  // the object that we dynamically allocate here
  // has to be deleted by a separate freearg typemap
  // (see below)
  $1 = new CCopasiAbstractArray::index_type;
  SV **tv;
  for(I32 i=0;i <= len;++i)
  {
    tv = av_fetch(tmpav,i,0);
    assert(*tv != NULL);
    $1->push_back(SvIV(*tv));
  }
};

// Since we allocate the temporary vector in the in typemap above on the heap,
// we need to make sure the object is deleted after the call.
%typemap (freearg) const CCopasiAbstractArray::index_type & MY_FUNNY_INDEX {
 delete $1;
}
#endif // SWIGPERL


%include "report/CCopasiArray.h"


%extend CCopasiAbstractArray {
   /* convert the operator[] to get methods */
   /* The argument needs a name that is unique so that I can write a specific typemap 
      (see above).    
   */ 
   virtual CCopasiAbstractArray::data_type get(const CCopasiAbstractArray::index_type & MY_FUNNY_INDEX)
   {
      return (*($self))[MY_FUNNY_INDEX];
   }
};

%template(AnnotatedFloatMatrix) CCopasiMatrixInterface<CMatrix<C_FLOAT64> >;
typedef CCopasiMatrixInterface<CMatrix<C_FLOAT64> > AnnotatedFloatMatrix;

