from PowhegControl import PowhegConfig_Dijet

# Use the Powheg_Dijet configuration
transform_runArgs = runArgs if 'runArgs' in dir() else None
transform_opts = opts if 'opts' in dir() else None
PowhegConfig = PowhegConfig_Dijet( runArgs=transform_runArgs, opts=transform_opts )

# if 'runArgs' in dir() :
#   PowhegConfig = PowhegConfig_Dijet(runArgs)
# else :
#   PowhegConfig = PowhegConfig_Dijet()
