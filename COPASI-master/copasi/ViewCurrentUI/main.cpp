// Copyright (C) 2012 - 2013 by Pedro Mendes, Virginia Tech Intellectual
// Properties, Inc., University of Heidelberg, and The University
// of Manchester.
// All rights reserved.

#include <stdexcept>

#include <QtGui/QApplication>
#include <QtCore/QString>

#define COPASI_MAIN

#include "copasi.h"

#include "UI/CQCopasiApplication.h"
#include "UI/copasiui3window.h"
#include "UI/DataModelGUI.h"
#include "UI/CQMessageBox.h"

#include "CopasiDataModel/CCopasiDataModel.h"
#include "report/CCopasiRootContainer.h"
#include "function/CFunctionDB.h"
#include "function/CFunction.h"
#include "commandline/COptionParser.h"
#include "commandline/COptions.h"

#ifdef COPASI_SBW_INTEGRATION
// SBW includes
#include "SBW/SBW.h"
#include <sstream>
#endif // COPASI_SBW_INTEGRATION

#ifdef Darwin
# include <QtCore/QDir>
# include "UI/qtUtilities.h"
#endif // Darwin

#include <worker.h>
#include <arguments.h>

int main(int argc, char *argv[])
{
  CQCopasiApplication a(argc, argv);

  // parse args
  Arguments args(argc, argv);

  if (args.handleCommandLine())
    {
      return 0;
    }

  a.setAttribute(Qt::AA_DontShowIconsInMenus, false);

  Q_INIT_RESOURCE(copasi);

  // Parse the commandline options
  try
    {
      // Create the root container but don't pass in any args past the executable name
      CCopasiRootContainer::init(1, argv, true);
    }
  catch (copasi::option_error & msg)
    {
      CQMessageBox::critical(NULL, "Initialization Error",
                             msg.what(),
                             QMessageBox::Ok , QMessageBox::Ok);

      return 1;
    }

#ifdef Darwin
  std::string PluginDir;

  COptions::getValue("CopasiDir", PluginDir);
  PluginDir += "/Contents/plugins";

  QApplication::setLibraryPaths(QStringList(FROM_UTF8(PluginDir)));
#endif // Darwin

  // Create the global data model.
  CCopasiRootContainer::addDatamodel();

  // instantiate model and apply all changes
  args.prepareModel();

  // Create the main application window.
  CopasiUI3Window *pWindow = CopasiUI3Window::create();
  Worker* pWorker = NULL;
#ifdef COPASI_SBW_INTEGRATION

  if (COptions::compareValue("SBWRegister", true))
    goto finish;

#endif // COPASI_SBW_INTEGRATION

  if (pWindow != NULL)
    {
      a.setMainWindow(pWindow);

      // pass control to the worker
      if (args.isValid())
        pWorker = new Worker(pWindow, &args);

      a.exec();
    }

finish:

  try // To suppress any access violations during destruction works only under Windows
    {
      if (pWorker != NULL)
        {
          delete pWorker;
          pWorker = NULL;
        }

      CCopasiRootContainer::destroy();
    }
  catch (...)
    {}

  return 0;
}
