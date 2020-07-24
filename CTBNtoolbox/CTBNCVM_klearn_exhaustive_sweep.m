function exitcode = CTBNCVM_klearn_exhaustive_sweep(simID)
%feel free to tailor function to your needs..
rng shuffle
fprintf('running\n');

L=input('Number of Nodes:');
M=input('Number of ODE-solver-iterations:');
H=input('Number of EM-iterations:');
MAX_PAR=input('Number of Parents:');
MAX_SWEEPS=input('Number of sweeps:');
MAX_PROCESSORS=input('Number of processors (for parallel computing):');

load('DATA.mat');
delete(gcp('nocreate'))
parpool(min(MAX_PROCESSORS,nchoosek(L,MAX_PAR)));

%%load prior over rates
load('prior.mat');

%%initial random net
B=rand(L,L);
for i=1:L
    B(i,i)=0;
end
B=B==max(B);
B=B';
[node]=createLibOfNodesGammaPrior(L,B, -ones(1,L),tau,M_T);
disp('Initial adjacency matrix:')
disp(node2net(node))

%hill-climb parameters
thresh=0.01; %convergence citerium
dt=0.005;   %time-step of ODE solver

%F output is marginal log-likelihood
%T_scores is graph scores
for i=1:MAX_SWEEPS
    [NET,T_SCORES,node] = SWEEPS_GREEDY_HILLCLIMB(MAX_PAR,L,node,Z,time,dt,M,H,tau,M_T,thresh);
    disp('Estimated adjacency matrix:')
    disp(node2net(node)) %print current best graph after sweep
    name=sprintf('%s_INFERRED_NET_K%d_SWEEP%d_LEARN_L%d_dt%3g.mat',simID,MAX_PAR,i,L,dt);
    save(name,'T_SCORES','NET','node','Z','time');
end
exitcode=0;

end