#include "GaudiKernel/DeclareFactoryEntries.h"

#include "FTK_RecTools/FTK_VertexFinderTool.h"
#include "FTK_RecTools/FTK_PixelClusterOnTrackTool.h"
#include "FTK_RecTools/FTK_SCTClusterOnTrackTool.h"

DECLARE_TOOL_FACTORY(FTK_VertexFinderTool)
DECLARE_TOOL_FACTORY(FTK_PixelClusterOnTrackTool)
DECLARE_TOOL_FACTORY(FTK_SCTClusterOnTrackTool)

DECLARE_FACTORY_ENTRIES(FTK_RecTools) {
  DECLARE_TOOL(FTK_VertexFinderTool)
  DECLARE_TOOL(FTK_PixelClusterOnTrackTool)
  DECLARE_TOOL(FTK_SCTClusterOnTrackTool)
}
