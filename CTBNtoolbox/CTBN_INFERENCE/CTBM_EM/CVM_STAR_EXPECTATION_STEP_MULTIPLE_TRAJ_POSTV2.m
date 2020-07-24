function [MU,RHO,TK,F] = CVM_STAR_EXPECTATION_STEP_MULTIPLE_TRAJ_POSTV2(node,DATA0,time0,dt,M)

d=size(DATA0{1});
L=d(3);
NO_TRAJ=length(time0);
F=0;

%set initial state
X0=zeros(2,L);
X0(1,:)=1/2;
X0(2,:)=1/2;

%set initial time
t0=21;

for j=1:NO_TRAJ
    timeVector=time0{j};
    statesMatrix=DATA0{j};
    
    [mu,rho,f,~] = P_CVMCTBN_STAR_POSTV3(node,2,dt,M,X0,t0,statesMatrix,timeVector,0.05);
    
    TK{j}=linspace(0,timeVector(end),ceil(timeVector(end)/dt));
    MU{j}=mu;
    RHO{j}=rho;
    F=F+f(end);
end

end

