# Input data from other simulations and research

**config.yaml** 
Example main config file with commented fields.

**large_model.msh** 

result of a hg model of a catchment scale at one test locality for the repository

**conc_flux_UOS.csv** 
expected concentration flux (?? or concentration) at the interface bentonite - rock
for times from 10 years up to 1e6 years, max at 39 years 
half of max at about 100years, timestep 1year up to 1000 years, after just few times with 
negligable conc.
!! Seems to be very short peak, what is the model for this.
What will be the time of leak? We should also model random leak times and compare it 

**samples_params.csv**
Set of parameters from Bayes inversion. For testing purpose.

**tunnel_mesh_cut_healed.msh**
2D model of TSX experiment. Outer boundery square [-50,50]^2, 
tunnel is cutted out as an ellipse 1.8m vertical semiaxis, 2.2m horizontal semiaxis
EDZ region elipse: 8m vertical, 10m horizontal (its boundary is present, but the region is not resolved)    

**rectangle_2x5.msh**
Test 2d mesh , no particular reason. 
