#include "GaudiKernel/DeclareFactoryEntries.h"

#include "SimplePIDSvc/SimplePIDSvc.h"

DECLARE_SERVICE_FACTORY( SimplePIDSvc )

DECLARE_FACTORY_ENTRIES( SimplePIDSvc ) { 
  DECLARE_SERVICE( SimplePIDSvc );
}
