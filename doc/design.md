# Main script
Two main scripts:

## bayes
Perform bayes inversion for parameters of the EDZ model.
Actions:

run - run inversion, collect accepted samples
    parameters should constrain: number of total samples, duration, processes

process - make plots, reports, file with samples

## mlmc
Perform MLMC forward uncertainty propagation.

Common arguments:
"<variants>" "<positions>"

<variants> - space separated varioant names, accept glob patterns
<positions> - space separated position indices, accept ranges "begin:end:step"

run - run MLMC sampling, limited either by target variance or by time and processes
process - make plots

plot indicators "<varaints>" "<positions>"
    
Other plot commands can be added.
