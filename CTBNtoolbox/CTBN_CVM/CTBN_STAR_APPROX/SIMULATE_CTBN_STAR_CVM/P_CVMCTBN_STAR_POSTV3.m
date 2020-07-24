function [MU,RHO,F,m] = P_CVMCTBN_STAR_POSTV3(node,D,dt,M,X0,t0,Z,TZ,thresh)
%%%%%%%%%%%%%CLUSTER VARIATIONAL STAR-APPROXIMATION FOR CTBNs WITH ERROR MODEL%%%%%
%%%%%%INPUT:%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%node is NODE STRUCT. DEFINING CTBN
%D is local dim. (only D=2 possible a.t.m.)
%T is SIMULATION TIME
%dt is time-step
%M is max. number of iterations for iterative ODE solver
%X0 is initial condition 
%Z is noisy observation of state 
%TZ is time of observation
%%%%%%OUTPUT:%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%mu(m,t,i,d) is MARGINAL PROB. IN m'th IERATION OF NODE i TO BE IN STATE d 
%AT TIME t
%rho(m,t,i,d) is LAGRANGE MULTIPLIER " " " " " " " " ".
%qt(t,i,d) is AVERAGED RATE """""""
%pst(t,i,d) is AVERAGED CHILD RESPONSE " " " " ".

%%%%%%%%%%%READ OUR SIM. PARAM.%%%%%%%%%%%%%%%%%%

L=length(node);
T=TZ(end)+dt*t0;
tau=ceil(TZ(end)/dt)+t0;
p=zeros(tau,L,2);
e0=squeeze(Z(:,1,:));

xt=linspace(0,T,tau);
%xt2=linspace(0,T-dt,tau-1);

%%%%%%%%%%SET INITIAL CONDITIONS%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%%%CREATE ANSATZ TRAJECTORY BY PICKING RND CIM%%%
mu_init(:,:,:)=p;
mu_init(1,1:L,1)=X0(1,:);
mu_init(1,1:L,2)=X0(2,:);

TSPAN = linspace(0,T,tau);
for i=1:L
    %%%PICK RANDOM CIM FOR ANSATZ SOLUTION
    q=zeros(tau,2);
    Q_i=node(i).cellOfCondRM;
    c=node(i).allPosStatesOfParents;
    sc=size(c);
    if isempty(c)
       Q=Q_i{1}; 
    else
       Q=Q_i{randi(sc(1))};
    end
    q(:,1)=0*Q(1,2).*ones(tau,1);
    q(:,2)=0*Q(2,1).*ones(tau,1);
    
    [Tt Y] = ode15s(@(t,y) CVM_CTBN_mu_fast(D,t, y, xt, q(1:tau,:), xt, ones(tau,2)), TSPAN, X0(:,i));
    mu_init(1:tau,i,:)=Y(1:tau,:);
end


rho=zeros(M,tau,L,2);
mu=zeros(M,tau,L,2);

for m=1:M
    mu(m,:,:,:)=squeeze(mu_init(:,:,:));
end

%%%%%%%%%%ACTUAL SIMULATION OF CTBN%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
psi0=zeros(tau,D);
TSPAN_R = linspace(T,0,tau); 
TSPAN_M = linspace(0,T,tau);

%options=odeset('RelTol',1e-12,'AbsTol',1e-13);
%%%%%%%%%%%%%%%%%%%%%%%%ITERATE%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
for m=1:M-1
    %%%%%%%%%%UPDATE LAGRANGE MULT. GIVEN MARG. PROBS%%%%%%%%%%%%%%%%%%%%%%
    for i=1:L
        [q,~] = P_CVMCTBN_EFF_RATES_V2(node,i,squeeze(mu(m,:,:,:)),D);
        
        if m>1
            psi0 = P_CVMCTBN_COUP_V2(node,i,squeeze(mu(m,:,:,:)),squeeze(rho(m-1,:,:,:)),D);
        end
        
        rhoT=ones(1,2);
        TSPAN_Y=linspace(T,TZ(end),t0); 
        %%%%%%%%%%ACTUAL TIME EVOLUTION OF LAGR. MULT. (BACKWARDS IN TIME)%
        %options = odeset('Jacobian',@(t,y)J_CVM_CTBN_mu(D,t, xt, q(1:tau,:), xt, squeeze(rho(m,1:tau,i,:))));
        [Tt Y] = ode15s(@(t,y) CVM_CTBN_rho_fast(D,t, y, xt, q, xt, psi0), TSPAN_Y, rhoT');
        R{length(TZ)+1}=Y;
        
     for k=length(TZ):-1:1
         
        rhoT=Y(end,:).*(squeeze(Z(:,k,i))');
        if k>1
          TSPAN_Y=linspace(TZ(k),TZ(k-1),ceil((TZ(k)-TZ(k-1))/dt)); 
        else
          TSPAN_Y=linspace(TZ(1),0,ceil(TZ(1)/dt));     
        end
       
        %options = odeset('Jacobian',@(t,y)J_CVM_CTBN_mu(D,t, xt, q(1:tau,:), xt, squeeze(rho(m,1:tau,i,:))));
        [Tt Y] = ode15s(@(t,y) CVM_CTBN_rho_fast(D,t, y, xt, q, xt, psi0), TSPAN_Y, rhoT');
        R{k}=Y;
        
         [msglast, msgidlast] = lastwarn;
         if isempty(msglast)==0
             Z(:,k,i)=1; %remove data-point
             warning('') %clear last warning
             msglast=[];
             sprintf('Warning: Could not process data-point %d of node %d',k,i)
         end
        
     end
     
     A=[];
     B=[];
     for k=length(TZ)+1:-1:1

         A=[A ;R{k}(:,1)];
         B=[B ;R{k}(:,2)];
   
     end
     rho(m,tau:-1:1,i,1)=A(1:tau);
     rho(m,tau:-1:1,i,2)=B(1:tau);
        %%%%%%%%%%UPDATE MARG. PROBS GIVEN LAGRANGE MULT. %%%%%%%%%%%%%%%%%
        %%%%%%%%%%ACTUAL TIME EVOLUTION OF MARG. PROB. (FORWARDS IN TIME)%%
      %  options = odeset('Jacobian',@(t,y)J_CVM_CTBN_mu(D,t, xt, q(1:tau,:), xt, squeeze(rho(m,1:tau,i,:))));
        [Tt Y] = ode15s(@(t,y) CVM_CTBN_mu_fast(D,t, y, xt, q(1:tau,:), xt, squeeze(rho(m,1:tau,i,:))), TSPAN_M, mu(m,1,i,:));
        mu(m+1,1:tau,i,:)=Y(1:tau,:);
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
MU=squeeze(mu(m,1:tau-t0,:,:));
RHO=squeeze(rho(m-1,1:tau-t0,:,:));

end



