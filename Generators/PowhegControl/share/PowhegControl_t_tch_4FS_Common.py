import PowhegControl
transform_runArgs = runArgs if "runArgs" in dir() else None
transform_opts = opts if "opts" in dir() else None
PowhegConfig = PowhegControl.PowhegControl(process_name="t_tch_4FS", run_args=transform_runArgs, run_opts=transform_opts)
