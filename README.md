# CTBNtoolbox
A toolbox for exact and approximate inference in continuous-time Bayesian networks.
This toolbox contains MATLAB implementations of methods for exact and approximate inference for Continuous-time Bayesian networks (CTBNs) as described in
"D. Linzner and H. Koeppl, Cluster Variational Approximations for Structure Learning of Continuous-Time Bayesian Networks from Incomplete Data. Advances in Neural Information Processing Systems 31, 7890--7900,2018".

We provide a couple of test scripts in the folder SCRIPTS for posterior inference and structure learning. With the exception of the IRMA test script (in SynthGRN_IRMA folder) these experiments are performed on synthetic data. 
For the IRMA test script, data is available at:

http://www.cell.com/supplemental/S0092-8674(09)00156-1 

Reference: I. Cantone et al., A Yeast Synthetic Network for In Vivo Assessment of Reverse-Engineering and Modeling Approaches, Cell, vol. 137, no. 1, pp. 172--181, 2009.

We provide a pre-processing script in SynthGRN_IRMA, implementing the observation model in the supplementary material of our manuscript.

Our structure learning method takes in data of the form of two cell-arrays:
likelihood of latent state given the observation model:
Z=cell{number of trajectory}(state,time-point,node index)
Time-points of measurement: 
TZ=cell{number of trajectory}(time-point)

IMPORTANT: Prior-rates for structure learning have to match time-scales of data provided. Otherwise our method will reject samples with the warning
"Warning: Could not process data-point x of node y".
We suggest re-scaling data to ~1 Transition per time unit and assuming a rate prior with alpha/beta ~ 1 

This toolbox depends in parts on the MATLAB statistics toolbox.
For parallel computing MATLABs parallel computing toolbox is needed.

Acknowledgements: Dominik Linzner is funded by the European Union's Horizon 2020 research and innovation programme under grant agreement 668858 (PrECISE).

-Dominik Linzner (dominik.linzner@bcs.tu-darmstadt.de)
