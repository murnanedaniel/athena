from Hto4lControl import Hto4lPowhegDefault

# Use the Powheg parameters for configuration
transform_runArgs   = runArgs if 'runArgs' in dir() else None
transform_opts      = opts if 'opts' in dir() else None
hto4lConfig4e    = Hto4lPowhegDefault( runArgs=transform_runArgs, opts=transform_opts )
hto4lConfig4mu   = Hto4lPowhegDefault( runArgs=transform_runArgs, opts=transform_opts )
hto4lConfig2e2mu = Hto4lPowhegDefault( runArgs=transform_runArgs, opts=transform_opts )
