function [mu,rho,gam] = exact_CTBN_evidence_POST(Q,dt,X0,t0,TZ,Z,D,L)
%exact calculation global markov process
T=TZ(end)+dt*t0;
tau=ceil(T/dt);
rho=zeros(tau,D^L);
mu=zeros(tau,D^L);


rhoT=ones(1,D^L);
TSPAN_Y=linspace(T,TZ(end),t0);
%%%%%%%%%%ACTUAL TIME EVOLUTION OF LAGR. MULT. (BACKWARDS IN TIME)%

drho=@(y) -Q*y;
[TY,Y]=ode15s(@(t,y) drho(y), TSPAN_Y, rhoT);
R{length(TZ)+1}=Y;

for k=length(TZ):-1:1
    
    rhoT=Y(end,:).*Z(k,:);
    if k>1
        TSPAN_Y=linspace(TZ(k),TZ(k-1),ceil((TZ(k)-TZ(k-1))/dt));
    else
        TSPAN_Y=linspace(TZ(1),0,ceil(TZ(1)/dt));
    end
    
    [TY,Y]=ode15s(@(t,y) drho(y), TSPAN_Y, rhoT');
    R{k}=Y;
    
end



for i=1:D^L
    A=[];
    for k=length(TZ)+1:-1:1
        A=[A ;R{k}(:,i)];
    end
    rho(tau:-1:1,i)=A(1:tau);
end


%TSPAN=linspace(T,0,ceil(T/dt));
%[TY,Y]=ode15s(@(t,y) drho(y), TSPAN, sparse(eT));
%rho(tau:-1:1,:)=Y;

Qt=zeros(tau-1,D^L,D^L);
R=zeros(tau-1,D^L,D^L);

irho(1:tau-1,:)=1./rho(1:tau-1,:);
for i=1:D^L
    for j=1:D^L
        R(1:tau-1,i,j)=squeeze(rho(1:tau-1,i).*irho(1:tau-1,j));
        Qt(1:tau-1,i,j)=Q(j,i)*R(1:tau-1,i,j);
    end
end
for i=1:D^L
    Qt(1:tau-1,i,i)=0;
    Qt(1:tau-1,i,i)=-sum(Qt(1:tau-1,:,i),2);
end

Y(1,:)=X0;
for t=2:tau-1
    Y(t,:)=expm(squeeze(dt*Qt(t,:,:)))*Y(t-1,:)';
end
mu(1:tau-1,:)=Y(1:tau-1,:);

gam=zeros(tau-1,D^L,D^L);
for i=1:D^L
    for j=1:D^L
        B(1:tau-1,i,j)=squeeze(Qt(1:tau-1,i,j)).*mu(1:tau-1,j);
    end
end
gam=B;

end

