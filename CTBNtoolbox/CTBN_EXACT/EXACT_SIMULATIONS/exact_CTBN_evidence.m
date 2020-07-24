function [mu,rho,gam] = exact_CTBN_evidence(Q,T,dt,e0,eT,D,L)
%exact calculation global markov process

tau=ceil(T/dt);

p=zeros(1,D^L);
mu_init=p;

for t=1:tau
    mu_init(t,:)=e0;
end

rho_init=p;
rho_init(end,:)=eT;

rho=zeros(tau,D^L);
mu=zeros(tau,D^L);

mu=squeeze(mu_init);
rho(end,:)=sparse(squeeze(rho_init(end,:)));

drho=@(y) -Q*y;
TSPAN=linspace(T,0,ceil(T/dt));
[TY,Y]=ode15s(@(t,y) drho(y), TSPAN, sparse(eT));

rho(tau:-1:1,:)=Y;


A=zeros(tau-1,D^L,D^L);
R=sparse(zeros(D^L,D^L));
for t=1:tau-1
    irho(t,:)=1./rho(t,:);
    R=squeeze(rho(t,:)'*irho(t,:));
   % A(t,:,:)=sparse(Q)'.*R;  
    A(t,:,:)=sparse(Q)'.*R;
end
for t=1:tau-1
    for i=1:D^L
        A(t,i,i)=0;
        A(t,i,i)=-sum(A(t,:,i));
    end
end

xt2=linspace(0,T-dt,tau-1);
[TY Y] = ode15s(@(t,y) dmu(D,L,t, y, xt2, A), xt2, sparse(e0)); % Solve ODE

mu(1:tau-1,:)=Y(1:tau-1,:);

gam=zeros(tau,D^L,D^L);
for t=1:tau-1
     irho(t,:)=1./rho(t,:);
    R=squeeze(rho(t,:)'*irho(t,:));
    %R=squeeze(irho(t,:)*rho(t,:)'); 
    A(t,:,:)=sparse(Q)'.*R;  
   % A(t,:,:)=sparse(Q).*R;  
end
for t=1:tau-1
    %B=squeeze(A(t,:,:)).*mu(t,:)';
    B=squeeze(A(t,:,:)).*mu(t,:);
    gam(t,:,:)=B;
end

end

