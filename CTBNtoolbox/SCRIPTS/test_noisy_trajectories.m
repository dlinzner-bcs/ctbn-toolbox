%Generate synthetic data
L=3;       %number of nodes
steps=10; %number of sampled transitions
dt=0.005; %simulation time step
SIGMA=0.8;%variance of gaussian observationmodel
%define ground-truth graph A (preset tree)
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
A(L,1)=1;
A=(A+1)'/2;

%parameters of ground-truth (kinetic ising-model)
ta=1; %rate scale
beta=0.6; % inverse temperature

%draw samples from ground-truth graph
node0=createLibOfNodes(L,A, -ones(1,L), beta,ta);
steps0=10; %number of sampled states
dt0=0.1; %granularity of sample trajectories

%draw sample trajectories
[TZ,DATA0] = CTBN_TRAJ(node0,1,steps,steps0,dt0,dt);
%add noise
[Z,DATA] = corrupted_observation_gaussian(DATA0,SIGMA);

%%
%calculate approximate posterior
M=20; %number of ODE-iterations
%%set initial state (see paper)
X0=zeros(2,L);
X0(1,:)=1/2;
X0(2,:)=1/2;
%%set initial time
t0=21;

[MU,RHO,F,~] = P_CVMCTBN_STAR_POSTV3(node0,2,dt,M,X0,t0,Z{1},TZ{1},0.01);

%%
%calculate exact posterior
[EZ,X0] = CTBN2CTMC_data(Z{1}); %transform CTBN into CTMC data

[Q] = jointIM_amalgation(node0); %transorm CTBM CIM into CTMC intensity matrices

%calculate exact posterior
[mu,rho,gam] = exact_CTBN_evidence_POST(Q,dt,X0,t0,TZ{1},EZ,2,L);

%calculate marginals
T=TZ{1}(end)+dt*t0;
tau=ceil(T/dt);
[mu_m,rho_m,gam_m] = componentize_markov_V2(tau-1,L,2,rho(1:tau-1,:),mu(1:tau-1,:),gam);
%%
%plot approximate and exact posterior
k=1; %node index

m=squeeze(MU(:,k,1)-MU(:,k,2)); %approx expected post. state
m_true=squeeze(mu_m(:,k,1)-mu_m(:,k,2));  %exact expected post. state
var=2*(abs(squeeze(MU(:,k,1).*(1-MU(:,k,1)))));
x=linspace(0,TZ{1}(end),length(m))';

clf
figure(1)
hold on
plot(TZ{1},DATA{1}(:,k),'o','LineWidth', 1)
% Plot posterior mean
plot(x,m, 'r', 'LineWidth', 2)
plot(x,m_true(1:(length(m_true)-t0+1)),':', 'LineWidth', 2,'Color','b')
ylim([-2,2])
xlabel('time')
ylabel('E[x]')
legend('data','approx. posterior','exact posterior')

curve1 = m+var;
curve2 = m-var;
% Plot posterior variance
plot(x, curve1, 'r', 'LineWidth', 0.1);
hold on;
plot(x, curve2, 'r', 'LineWidth', 0.1);
x2 = [x; flipud(x)]; % Use ; instead of ,
inBetween = [curve1; flipud(curve2)]; % Use ; instead of ,
fill(x2, inBetween, [0.95 0.95 0.95]);
grid on;
plot(x,m, 'r', 'LineWidth', 2)
plot(x,m_true(1:(length(m_true)-t0+1)),':', 'LineWidth', 2,'Color','b')

