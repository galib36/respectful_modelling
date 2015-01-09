/* Begin CVS Header
   $Source: /Volumes/Home/Users/shoops/cvs/copasi_dev/copasi/utilities/CWriteConfig.h,v $
   $Revision: 1.10 $
   $Name:  $
   $Author: shoops $
   $Date: 2006/04/27 01:32:43 $
   End CVS Header */

// Copyright � 2005 by Pedro Mendes, Virginia Tech Intellectual
// Properties, Inc. and EML Research, gGmbH.
// All rights reserved.

/**
 *  CWriteConfig class. A more elaborate class description.
 *
 *  New Class based on pmutils read functionality
 * (C) Stefan Hoops 2001
 */

#ifndef COPASI_CWriteConfig
#define COPASI_CWriteConfig

#include <sstream>

class CWriteConfig
  {
  public:
    /**
     *  Default consructor.
     *  This creates a configuration buffer without assigning a filename.
     *  It is currently useless.
     */
    CWriteConfig();

    /**
     *  Specified consructor.
     *  This opens the configuration file using the filename specified
     *  as the argument.
     *  @param name name of the confguration file.
     *  @param mode output mode, defaults to creating a new file.
     */
    CWriteConfig(const std::string& name,
                 std::ios_base::openmode mode = std::ios_base::out);

    /**
     *  Destructor.
     *  The destructor calls the method flush().
     */
    ~CWriteConfig();

    /**
     *  Flush the output buffer to the configuration file.
     *  @return mFail
     *  @see mFail
     */
    C_INT32 flush();

    /**
     *  Returns the failure status.
     *  @return mFail
     *  @see mFail
     */
    C_INT32 fail();

    /**
     *  Writes a variable to the output file.
     *  @param name name of the variable to be written.
     *  @param type type of the variable to be written.
     *  @param *pout pointer to the location where the variable
     *               is stored.
     *  @return mFail
     *  @see mFail
     */
    C_INT32 setVariable(const std::string & name,
                        const std::string & type,
                        const void * pout);
    C_INT32 setVariable(const std::string & name,
                        const std::string & type,
                        const void * pout1,
                        const void * pout2);

    /**
     *  Save various default values
     */
    void setDefaults();

  private:
    /**
     *  Commits all information to the configuration file.
     *  This method is called by the destructor.
     */
    C_INT32 commit();

    /**
     *  Writes a version number to the file.
     *  This method is called by one of the constructors.
     */
    void writeVersion(void);

    /**
     *  Name of the configuration file.
     */
    std::string mFileName;               // Config File Name

    /**
     *  Output file
     */
    std::stringstream mBuffer;

    /**
     *  Current line number in the configuration file
     */
    unsigned C_INT32 mLineNumber;             // Current Line Number

    /**
     *  Open mode
     */
    std::ios_base::openmode mOpenMode;

    /**
     *  Failure status:
     *  0 = no error
     *  !0 = error
     */
    C_INT32 mFail;                   // Failure State
  };
#endif // COPASI_CWriteConfig
