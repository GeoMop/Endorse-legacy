# Fine scale local models

Model of single container position with EDZ and possible heterogenities, no fractures.
=> from application of several boundary conditions determination of ekvivalent properties at neigbourhood of given points

# Fine scale

Fully resolved center storage borehole and its EDZ, meshing governed by a distance field from the borehole surface. 
Some fractures.
Side boreholes not resolved but replaced by ekvivalent properties 

# First coarsening

Replace by model with tunnels replaced by eqivalent properties, elements about 5m, only major fractures.

# Second coarsening

Possible replacement of major fractures. More advanced, possible usage of analytical determination of the equivalent properties around fractures.

Neccessary tools:
- mark elements of the coarse mesh containing the fine scale models -> assing particular regions -> calculate the anisotropic tensors on these regions, pass them to Flow123d
