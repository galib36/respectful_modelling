// Begin CVS Header 
//   $Source: /Volumes/Home/Users/shoops/cvs/copasi_dev/CopasiWS/copasiws_dev/copasiws/parameterestimationws/src/org/copasi/copasiws/services/parameterestimationws/types/SetFitItemsAndMethodDocument.java,v $ 
//   $Revision: 1.1 $ 
//   $Name:  $ 
//   $Author: jdada $ 
//   $Date: 2008/11/11 12:19:51 $ 
// End CVS Header 

// Copyright (C) 2008 by Pedro Mendes, Virginia Tech Intellectual 
// Properties, Inc., EML Research, gGmbH, University of Heidelberg, 
// and The University of Manchester. 
// All rights reserved. 

/*
 * An XML document type.
 * Localname: SetFitItemsAndMethod
 * Namespace: http://parameterestimationws.services.copasiws.copasi.org/types
 * Java type: org.copasi.copasiws.services.parameterestimationws.types.SetFitItemsAndMethodDocument
 *
 * Automatically generated - do not modify.
 */
package org.copasi.copasiws.services.parameterestimationws.types;


/**
 * A document containing one SetFitItemsAndMethod(@http://parameterestimationws.services.copasiws.copasi.org/types) element.
 *
 * This is a complex type.
 */
public interface SetFitItemsAndMethodDocument extends org.apache.xmlbeans.XmlObject {
    public static final org.apache.xmlbeans.SchemaType type = (org.apache.xmlbeans.SchemaType) org.apache.xmlbeans.XmlBeans.typeSystemForClassLoader(SetFitItemsAndMethodDocument.class.getClassLoader(),
            "schemaorg_apache_xmlbeans.system.parameterestimationwsTypeSystem")
                                                                                                                           .resolveHandle("setfititemsandmethodfff2doctype");

    /**
     * Gets the "SetFitItemsAndMethod" element
     */
    org.copasi.copasiws.services.parameterestimationws.types.SetFitItemsAndMethodDocument.SetFitItemsAndMethod getSetFitItemsAndMethod();

    /**
     * Sets the "SetFitItemsAndMethod" element
     */
    void setSetFitItemsAndMethod(
        org.copasi.copasiws.services.parameterestimationws.types.SetFitItemsAndMethodDocument.SetFitItemsAndMethod setFitItemsAndMethod);

    /**
     * Appends and returns a new empty "SetFitItemsAndMethod" element
     */
    org.copasi.copasiws.services.parameterestimationws.types.SetFitItemsAndMethodDocument.SetFitItemsAndMethod addNewSetFitItemsAndMethod();

    /**
     * An XML SetFitItemsAndMethod(@http://parameterestimationws.services.copasiws.copasi.org/types).
     *
     * This is a complex type.
     */
    public interface SetFitItemsAndMethod extends org.apache.xmlbeans.XmlObject {
        public static final org.apache.xmlbeans.SchemaType type = (org.apache.xmlbeans.SchemaType) org.apache.xmlbeans.XmlBeans.typeSystemForClassLoader(SetFitItemsAndMethod.class.getClassLoader(),
                "schemaorg_apache_xmlbeans.system.parameterestimationwsTypeSystem")
                                                                                                                               .resolveHandle("setfititemsandmethodebefelemtype");

        /**
         * Gets the "resourceId" element
         */
        int getResourceId();

        /**
         * Gets (as xml) the "resourceId" element
         */
        org.apache.xmlbeans.XmlInt xgetResourceId();

        /**
         * Sets the "resourceId" element
         */
        void setResourceId(int resourceId);

        /**
         * Sets (as xml) the "resourceId" element
         */
        void xsetResourceId(org.apache.xmlbeans.XmlInt resourceId);

        /**
         * Gets the "userId" element
         */
        java.lang.String getUserId();

        /**
         * Gets (as xml) the "userId" element
         */
        org.apache.xmlbeans.XmlString xgetUserId();

        /**
         * Sets the "userId" element
         */
        void setUserId(java.lang.String userId);

        /**
         * Sets (as xml) the "userId" element
         */
        void xsetUserId(org.apache.xmlbeans.XmlString userId);

        /**
         * Gets array of all "itemToFit" elements
         */
        org.copasi.copasiws.services.parameterestimationws.types.ItemToFit[] getItemToFitArray();

        /**
         * Gets ith "itemToFit" element
         */
        org.copasi.copasiws.services.parameterestimationws.types.ItemToFit getItemToFitArray(
            int i);

        /**
         * Returns number of "itemToFit" element
         */
        int sizeOfItemToFitArray();

        /**
         * Sets array of all "itemToFit" element
         */
        void setItemToFitArray(
            org.copasi.copasiws.services.parameterestimationws.types.ItemToFit[] itemToFitArray);

        /**
         * Sets ith "itemToFit" element
         */
        void setItemToFitArray(int i,
            org.copasi.copasiws.services.parameterestimationws.types.ItemToFit itemToFit);

        /**
         * Inserts and returns a new empty value (as xml) as the ith "itemToFit" element
         */
        org.copasi.copasiws.services.parameterestimationws.types.ItemToFit insertNewItemToFit(
            int i);

        /**
         * Appends and returns a new empty value (as xml) as the last "itemToFit" element
         */
        org.copasi.copasiws.services.parameterestimationws.types.ItemToFit addNewItemToFit();

        /**
         * Removes the ith "itemToFit" element
         */
        void removeItemToFit(int i);

        /**
         * Gets the "method" element
         */
        org.copasi.copasiws.services.parameterestimationws.types.Method getMethod();

        /**
         * True if has "method" element
         */
        boolean isSetMethod();

        /**
         * Sets the "method" element
         */
        void setMethod(
            org.copasi.copasiws.services.parameterestimationws.types.Method method);

        /**
         * Appends and returns a new empty "method" element
         */
        org.copasi.copasiws.services.parameterestimationws.types.Method addNewMethod();

        /**
         * Unsets the "method" element
         */
        void unsetMethod();

        /**
         * A factory class with static methods for creating instances
         * of this type.
         */
        public static final class Factory {
            private Factory() {
            } // No instance of this class allowed

            public static org.copasi.copasiws.services.parameterestimationws.types.SetFitItemsAndMethodDocument.SetFitItemsAndMethod newInstance() {
                return (org.copasi.copasiws.services.parameterestimationws.types.SetFitItemsAndMethodDocument.SetFitItemsAndMethod) org.apache.xmlbeans.XmlBeans.getContextTypeLoader()
                                                                                                                                                                .newInstance(type,
                    null);
            }

            public static org.copasi.copasiws.services.parameterestimationws.types.SetFitItemsAndMethodDocument.SetFitItemsAndMethod newInstance(
                org.apache.xmlbeans.XmlOptions options) {
                return (org.copasi.copasiws.services.parameterestimationws.types.SetFitItemsAndMethodDocument.SetFitItemsAndMethod) org.apache.xmlbeans.XmlBeans.getContextTypeLoader()
                                                                                                                                                                .newInstance(type,
                    options);
            }
        }
    }

    /**
     * A factory class with static methods for creating instances
     * of this type.
     */
    public static final class Factory {
        private Factory() {
        } // No instance of this class allowed

        public static org.copasi.copasiws.services.parameterestimationws.types.SetFitItemsAndMethodDocument newInstance() {
            return (org.copasi.copasiws.services.parameterestimationws.types.SetFitItemsAndMethodDocument) org.apache.xmlbeans.XmlBeans.getContextTypeLoader()
                                                                                                                                       .newInstance(type,
                null);
        }

        public static org.copasi.copasiws.services.parameterestimationws.types.SetFitItemsAndMethodDocument newInstance(
            org.apache.xmlbeans.XmlOptions options) {
            return (org.copasi.copasiws.services.parameterestimationws.types.SetFitItemsAndMethodDocument) org.apache.xmlbeans.XmlBeans.getContextTypeLoader()
                                                                                                                                       .newInstance(type,
                options);
        }

        /** @param xmlAsString the string value to parse */
        public static org.copasi.copasiws.services.parameterestimationws.types.SetFitItemsAndMethodDocument parse(
            java.lang.String xmlAsString)
            throws org.apache.xmlbeans.XmlException {
            return (org.copasi.copasiws.services.parameterestimationws.types.SetFitItemsAndMethodDocument) org.apache.xmlbeans.XmlBeans.getContextTypeLoader()
                                                                                                                                       .parse(xmlAsString,
                type, null);
        }

        public static org.copasi.copasiws.services.parameterestimationws.types.SetFitItemsAndMethodDocument parse(
            java.lang.String xmlAsString, org.apache.xmlbeans.XmlOptions options)
            throws org.apache.xmlbeans.XmlException {
            return (org.copasi.copasiws.services.parameterestimationws.types.SetFitItemsAndMethodDocument) org.apache.xmlbeans.XmlBeans.getContextTypeLoader()
                                                                                                                                       .parse(xmlAsString,
                type, options);
        }

        /** @param file the file from which to load an xml document */
        public static org.copasi.copasiws.services.parameterestimationws.types.SetFitItemsAndMethodDocument parse(
            java.io.File file)
            throws org.apache.xmlbeans.XmlException, java.io.IOException {
            return (org.copasi.copasiws.services.parameterestimationws.types.SetFitItemsAndMethodDocument) org.apache.xmlbeans.XmlBeans.getContextTypeLoader()
                                                                                                                                       .parse(file,
                type, null);
        }

        public static org.copasi.copasiws.services.parameterestimationws.types.SetFitItemsAndMethodDocument parse(
            java.io.File file, org.apache.xmlbeans.XmlOptions options)
            throws org.apache.xmlbeans.XmlException, java.io.IOException {
            return (org.copasi.copasiws.services.parameterestimationws.types.SetFitItemsAndMethodDocument) org.apache.xmlbeans.XmlBeans.getContextTypeLoader()
                                                                                                                                       .parse(file,
                type, options);
        }

        public static org.copasi.copasiws.services.parameterestimationws.types.SetFitItemsAndMethodDocument parse(
            java.net.URL u)
            throws org.apache.xmlbeans.XmlException, java.io.IOException {
            return (org.copasi.copasiws.services.parameterestimationws.types.SetFitItemsAndMethodDocument) org.apache.xmlbeans.XmlBeans.getContextTypeLoader()
                                                                                                                                       .parse(u,
                type, null);
        }

        public static org.copasi.copasiws.services.parameterestimationws.types.SetFitItemsAndMethodDocument parse(
            java.net.URL u, org.apache.xmlbeans.XmlOptions options)
            throws org.apache.xmlbeans.XmlException, java.io.IOException {
            return (org.copasi.copasiws.services.parameterestimationws.types.SetFitItemsAndMethodDocument) org.apache.xmlbeans.XmlBeans.getContextTypeLoader()
                                                                                                                                       .parse(u,
                type, options);
        }

        public static org.copasi.copasiws.services.parameterestimationws.types.SetFitItemsAndMethodDocument parse(
            java.io.InputStream is)
            throws org.apache.xmlbeans.XmlException, java.io.IOException {
            return (org.copasi.copasiws.services.parameterestimationws.types.SetFitItemsAndMethodDocument) org.apache.xmlbeans.XmlBeans.getContextTypeLoader()
                                                                                                                                       .parse(is,
                type, null);
        }

        public static org.copasi.copasiws.services.parameterestimationws.types.SetFitItemsAndMethodDocument parse(
            java.io.InputStream is, org.apache.xmlbeans.XmlOptions options)
            throws org.apache.xmlbeans.XmlException, java.io.IOException {
            return (org.copasi.copasiws.services.parameterestimationws.types.SetFitItemsAndMethodDocument) org.apache.xmlbeans.XmlBeans.getContextTypeLoader()
                                                                                                                                       .parse(is,
                type, options);
        }

        public static org.copasi.copasiws.services.parameterestimationws.types.SetFitItemsAndMethodDocument parse(
            java.io.Reader r)
            throws org.apache.xmlbeans.XmlException, java.io.IOException {
            return (org.copasi.copasiws.services.parameterestimationws.types.SetFitItemsAndMethodDocument) org.apache.xmlbeans.XmlBeans.getContextTypeLoader()
                                                                                                                                       .parse(r,
                type, null);
        }

        public static org.copasi.copasiws.services.parameterestimationws.types.SetFitItemsAndMethodDocument parse(
            java.io.Reader r, org.apache.xmlbeans.XmlOptions options)
            throws org.apache.xmlbeans.XmlException, java.io.IOException {
            return (org.copasi.copasiws.services.parameterestimationws.types.SetFitItemsAndMethodDocument) org.apache.xmlbeans.XmlBeans.getContextTypeLoader()
                                                                                                                                       .parse(r,
                type, options);
        }

        public static org.copasi.copasiws.services.parameterestimationws.types.SetFitItemsAndMethodDocument parse(
            javax.xml.stream.XMLStreamReader sr)
            throws org.apache.xmlbeans.XmlException {
            return (org.copasi.copasiws.services.parameterestimationws.types.SetFitItemsAndMethodDocument) org.apache.xmlbeans.XmlBeans.getContextTypeLoader()
                                                                                                                                       .parse(sr,
                type, null);
        }

        public static org.copasi.copasiws.services.parameterestimationws.types.SetFitItemsAndMethodDocument parse(
            javax.xml.stream.XMLStreamReader sr,
            org.apache.xmlbeans.XmlOptions options)
            throws org.apache.xmlbeans.XmlException {
            return (org.copasi.copasiws.services.parameterestimationws.types.SetFitItemsAndMethodDocument) org.apache.xmlbeans.XmlBeans.getContextTypeLoader()
                                                                                                                                       .parse(sr,
                type, options);
        }

        public static org.copasi.copasiws.services.parameterestimationws.types.SetFitItemsAndMethodDocument parse(
            org.w3c.dom.Node node) throws org.apache.xmlbeans.XmlException {
            return (org.copasi.copasiws.services.parameterestimationws.types.SetFitItemsAndMethodDocument) org.apache.xmlbeans.XmlBeans.getContextTypeLoader()
                                                                                                                                       .parse(node,
                type, null);
        }

        public static org.copasi.copasiws.services.parameterestimationws.types.SetFitItemsAndMethodDocument parse(
            org.w3c.dom.Node node, org.apache.xmlbeans.XmlOptions options)
            throws org.apache.xmlbeans.XmlException {
            return (org.copasi.copasiws.services.parameterestimationws.types.SetFitItemsAndMethodDocument) org.apache.xmlbeans.XmlBeans.getContextTypeLoader()
                                                                                                                                       .parse(node,
                type, options);
        }

        /** @deprecated {@link XMLInputStream} */
        public static org.copasi.copasiws.services.parameterestimationws.types.SetFitItemsAndMethodDocument parse(
            org.apache.xmlbeans.xml.stream.XMLInputStream xis)
            throws org.apache.xmlbeans.XmlException,
                org.apache.xmlbeans.xml.stream.XMLStreamException {
            return (org.copasi.copasiws.services.parameterestimationws.types.SetFitItemsAndMethodDocument) org.apache.xmlbeans.XmlBeans.getContextTypeLoader()
                                                                                                                                       .parse(xis,
                type, null);
        }

        /** @deprecated {@link XMLInputStream} */
        public static org.copasi.copasiws.services.parameterestimationws.types.SetFitItemsAndMethodDocument parse(
            org.apache.xmlbeans.xml.stream.XMLInputStream xis,
            org.apache.xmlbeans.XmlOptions options)
            throws org.apache.xmlbeans.XmlException,
                org.apache.xmlbeans.xml.stream.XMLStreamException {
            return (org.copasi.copasiws.services.parameterestimationws.types.SetFitItemsAndMethodDocument) org.apache.xmlbeans.XmlBeans.getContextTypeLoader()
                                                                                                                                       .parse(xis,
                type, options);
        }

        /** @deprecated {@link XMLInputStream} */
        public static org.apache.xmlbeans.xml.stream.XMLInputStream newValidatingXMLInputStream(
            org.apache.xmlbeans.xml.stream.XMLInputStream xis)
            throws org.apache.xmlbeans.XmlException,
                org.apache.xmlbeans.xml.stream.XMLStreamException {
            return org.apache.xmlbeans.XmlBeans.getContextTypeLoader()
                                               .newValidatingXMLInputStream(xis,
                type, null);
        }

        /** @deprecated {@link XMLInputStream} */
        public static org.apache.xmlbeans.xml.stream.XMLInputStream newValidatingXMLInputStream(
            org.apache.xmlbeans.xml.stream.XMLInputStream xis,
            org.apache.xmlbeans.XmlOptions options)
            throws org.apache.xmlbeans.XmlException,
                org.apache.xmlbeans.xml.stream.XMLStreamException {
            return org.apache.xmlbeans.XmlBeans.getContextTypeLoader()
                                               .newValidatingXMLInputStream(xis,
                type, options);
        }
    }
}
