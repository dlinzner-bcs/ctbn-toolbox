%set experiment parameters
%
L=3; %number of nodes of random graph

%exhaustive search in small graphs will take approx 10mins
%CAREFULL! graph search scales super-exponentially in L
%for large graphs (~10nodes) only greedy searches(by restricting possible families)
%are possible
%parallelized scoring allows network inference of large graphs 

N_TRAJ=5; %number of synthetic trajectories
SIGMA=0.1; %noise in synthetic trajectories
MAX_PAR=2; %maximum number of parents in synthetic experiments
MAX_SWEEPS=3; %maximum number of sweeps

%set number of parallel computers
%delete(gcp('nocreate'))
%parpool(min(16,L));

%%
%Generate synthetic data
steps=5; %number of sampled transitions
dt=0.005; %simulation time step

%define ground-truth graph A (preset tree)
%%initial net
A=-ones(L);
for i=1:L/2-1 
    A(i,2*i) =1;
    A(i,2*i+1) =1;  
end
if mod(L,2)==0
    A(L/2,L)=1;
else
   A(floor(L/2),L)=1;
   A(floor(L/2),L-1)=1;
end

%A(L,1)=1;
A=(A+1)'/2;

%parameters of ground-truth (kinetic ising-model)
ta=1; %rate scale
beta=0.6; % inverse temperature

%draw samples from ground-truth graph
node0=createLibOfNodes(L,A, -ones(1,L), beta,ta);
steps0=5; %number of sampled states
dt0=0.01; %granularity of sample trajectories

%draw sample trajectories
[time0,DATA0] = CTBN_TRAJ(node0,N_TRAJ,steps,steps0,dt0,dt);
%add noise
[DATAC,~] = corrupted_observation_gaussian(DATA0,SIGMA);
%%
%%define prior over rates
tau=ones(2^MAX_PAR,L,2)*10;
M_T=ones(2^MAX_PAR,L,2,2)*5;

%%initial random net
B=rand(L,L);
for i=1:L
    B(i,i)=0;
end
B=B==max(B);
B=B';

[node]=createLibOfNodesGammaPrior(L,B, -ones(1,L),tau,M_T);
disp('Ground-truth adjacency matrix:')
disp(node2net(node0))
disp('Initial adjacency matrix:')
disp(node2net(node))

%%
%hill-climb parameters
thresh=0.01; %convergence citerion
H=10;        %number of EM-iterations
M=10;       %number of ODE-iterations

%F output is marginal log-likelihood
%T_scores is graph scores
for i=1:MAX_SWEEPS
[NET,T_SCORES,node] = SWEEPS_GREEDY_HILLCLIMB(MAX_PAR,L,node,DATAC,time0,dt,M,H,tau,M_T,thresh);
disp('Estimated adjacency matrix:')
disp(node2net(node)) %print current best graph after sweep
end