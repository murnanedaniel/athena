#include "GaudiKernel/DeclareFactoryEntries.h"
#include "DerivationFrameworkTau/TauSelectionWrapper.h"
#include "DerivationFrameworkTau/TauTruthMatchingWrapper.h"
#include "DerivationFrameworkTau/TauPVRefitTool.h"
#include "DerivationFrameworkTau/TauPVTrkSelectionTool.h"
#include "DerivationFrameworkTau/TauOverlappingElectronLLHDecoratorWrapper.h"
#include "DerivationFrameworkTau/TauPFOCalHitDecorator.h"
#include "DerivationFrameworkTau/DiTauMassDecorator.h"
#include "DerivationFrameworkTau/DiTauProngDecorator.h"

using namespace DerivationFramework;

DECLARE_TOOL_FACTORY( TauSelectionWrapper )
DECLARE_TOOL_FACTORY( TauTruthMatchingWrapper )
DECLARE_TOOL_FACTORY( TauPVRefitTool )
DECLARE_TOOL_FACTORY( TauPVTrkSelectionTool )
DECLARE_TOOL_FACTORY( TauOverlappingElectronLLHDecoratorWrapper )
DECLARE_TOOL_FACTORY( TauPFOCalHitDecorator )
DECLARE_TOOL_FACTORY( DiTauMassDecorator )
DECLARE_TOOL_FACTORY( DiTauProngDecorator )

DECLARE_FACTORY_ENTRIES( DerivationFrameworkTau ) {
   DECLARE_TOOL( TauSelectionWrapper )
   DECLARE_TOOL( TauTruthMatchingWrapper )
   DECLARE_TOOL( TauPVRefitTool )
   DECLARE_TOOL( TauPVTrkSelectionTool )
   DECLARE_TOOL( TauOverlappingElectronLLHDecoratorWrapper )
   DECLARE_TOOL( TauPFOCalHitDecorator )
   DECLARE_TOOL( DiTauMassDecorator )
   DECLARE_TOOL( DiTauProngDecorator )
}
