all tests (with exception of GUI) should run in a Docker/Singularity image using possibly 
the [swrap](https://github.com/flow123d/swrap) tool.

`archive_mlmc.sh` - for given directory compress or extract sampling subdirectory `output`
`meshes` - tests for mesh creation functions
`test-data` - various reference input data for the tests
`transport_vis.pvsm` - paraview pipeline for visualization of the transport model
`transport_vis.py` - same pipeline but in Python (could be used to modify view parametricaly)
