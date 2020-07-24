%Structure learning from time-course data
%Performs an exhaustive graph-search with families of up to k parents
%input: Simulation ID
%function will request input:
%Number of nodes: Number of variables
%Number of ODE-solver-iterations (default 20)
%Number of EM-iterations (default 10)
%Number of parents: maximum numbers of parents searched
%Number of sweeps: Number of family searches per nodes (default 3)

%function will require a file of the name 'DATA.mat' that contains the data
%Z: Noisy state observations of dimensions 
%Z{number of trajectories}(size of local state-space,number of time points, number of nodes)
%TZ: Time points of observations of dimensions
%TZ{number of trajectories}(number of time points)

%function will require a file of the name 'priot.mat' that contains the
%rate prior. We added a default prior for IRMA-dataset:
%tau=ones(2^4,5,2)*10;
%M_T=ones(2^4,5,2,2)*5;
%In general (uniformative prior):
%tau=ones(2^MAX_PAR,L,size of local state-space)*tau_0;
%M_T=ones(2^MAX_PAR,L,size of local state-space,size of local state-space)*alpha_0;

%We suggest re-scaling time-scales of data to ~1 transition per time-unit
%and tau_0=10 and alpha_0=5;

%If data-points are at a timescale very different as indicated from rate prior 
%function will neglect those data-points and throw a warning
%some data-points may also be unreachable for a given graph, a few warnings
%may be normal during search

exitcode = CTBNCVM_klearn_exhaustive_sweep('IRMA_test');