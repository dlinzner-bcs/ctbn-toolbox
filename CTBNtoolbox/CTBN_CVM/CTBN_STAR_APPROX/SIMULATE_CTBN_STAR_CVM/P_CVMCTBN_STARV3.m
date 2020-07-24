function [MU,RHO,F,m] = P_CVMCTBN_STARV3(node,D,T,dt,M,e0,eT,thresh)
%%%%%%%%%%%%%CLUSTER VARIATIONAL STAR-APPROXIMATION FOR CTBNs%%%%%%%%%%%%%%
%%%%%%INPUT:%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%node is NODE STRUCT. DEFINING CTBN
%D is local dim. (only D=2 possible a.t.m.)
%T is SIMULATION TIME
%dt is time-step
%M is max. number of iterations for iterative ODE solver
%e0 is initial condition
%eT is condition at end of sim. interval
%%%%%%OUTPUT:%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%mu(m,t,i,d) is MARGINAL PROB. IN m'th IERATION OF NODE i TO BE IN STATE d
%AT TIME t
%rho(m,t,i,d) is LAGRANGE MULTIPLIER " " " " " " " " ".
%qt(t,i,d) is AVERAGED RATE """""""
%pst(t,i,d) is AVERAGED CHILD RESPONSE " " " " ".

%%%%%%%%%%%READ OUR SIM. PARAM.%%%%%%%%%%%%%%%%%%

L=length(node);
tau=ceil(T/dt);
p=zeros(tau,L,2);

xt=linspace(0,T,tau);
xt2=linspace(0,T-dt,tau-1);

%%%%%%%%%%SET INITIAL CONDITIONS%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%%%CREATE ANSATZ TRAJECTORY BY PICKING RND CIM%%%
mu_init(:,:,:)=p;
mu_init(1,1:L,1)=e0(1,:);
mu_init(1,1:L,2)=e0(2,:);

TSPAN = linspace(0,T-dt,tau-1);
for i=1:L
    %%%PICK RANDOM CIM FOR ANSATZ SOLUTION
    q=zeros(tau,2);
    Q_i=node(i).cellOfCondRM;
    c=node(i).allPosStatesOfParents;
    sc=size(c);
    if sc(1)>0
        Q=Q_i{randi(sc(1))};
    else
        Q=Q_i{1};
    end
    q(:,1)=Q(1,2).*ones(tau,1);
    q(:,2)=Q(2,1).*ones(tau,1);
    
    [Tt Y] = ode15s(@(t,y) CVM_CTBN_mu_fast(D,t, y, xt2, q(1:tau-1,:), xt2, ones(tau-1,2)), TSPAN, e0(:,i));
    mu_init(1:tau-1,i,:)=Y(1:tau-1,:);
end

rho_init(:,:,:)=p;
rho_init(end,1:L,1)=eT(1,:);
rho_init(end,1:L,2)=eT(2,:);

rho=zeros(M,tau,L,2);
mu=zeros(M,tau,L,2);

for m=1:M
    mu(m,:,:,:)=squeeze(mu_init(:,:,:));
    rho(m,end,:,:)=squeeze(rho_init(end,:,:));
end

%%%%%%%%%%ACTUAL SIMULATION OF CTBN%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
psi0=zeros(tau,D);
TSPAN_R = linspace(T,0,tau);
TSPAN_M = linspace(0,T-dt,tau-1);
%%%%%%%%%%%%%%%%%%%%%%%%ITERATE%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
for m=1:M-1
    %%%%%%%%%%UPDATE LAGRANGE MULT. GIVEN MARG. PROBS%%%%%%%%%%%%%%%%%%%%%%
    for i=1:L
        [q,~] = P_CVMCTBN_EFF_RATES_V2(node,i,squeeze(mu(m,:,:,:)),D);
        if m>1
            psi0 = P_CVMCTBN_COUP_V2(node,i,squeeze(mu(m,:,:,:)),squeeze(rho(m-1,:,:,:)),D);
        end
        %%%%%%%%%%ACTUAL TIME EVOLUTION OF LAGR. MULT. (BACKWARDS IN TIME)%
        [Tt Y] = ode15s(@(t,y) CVM_CTBN_rho_fast(D,t, y, xt, q, xt, psi0), TSPAN_R, rho(m,end,i,:));
        rho(m,tau:-1:1,i,:)=Y;
        %%%%%%%%%%UPDATE MARG. PROBS GIVEN LAGRANGE MULT. %%%%%%%%%%%%%%%%%
        %%%%%%%%%%ACTUAL TIME EVOLUTION OF MARG. PROB. (FORWARDS IN TIME)%%
        [Tt Y] = ode15s(@(t,y) CVM_CTBN_mu_fast(D,t, y, xt2, q(1:tau-1,:), xt2, squeeze(rho(m,1:tau-1,i,:))), TSPAN_M, mu(m,1,i,:));
        mu(m+1,1:tau-1,i,:)=Y(1:tau-1,:);
    end
    %%%%%CALCULATE LIKELIHOOD LOWER BOUND TO CHECK FOR CONVERGENCE
    if m>1
        [F(m),~,~] = CVM_CTBN_likelihood_star(node,squeeze(mu(m,:,:,:)),squeeze(rho(m-1,:,:,:)),dt);
    end
    if m>5
        stdF=std(F(end:-1:end-4));
        if stdF<thresh
            break
        end
        
    end
    
end
%%%%CONVERGED SOLUTIONS
MU=squeeze(mu(m,:,:,:));
RHO=squeeze(rho(m-1,:,:,:));

end



