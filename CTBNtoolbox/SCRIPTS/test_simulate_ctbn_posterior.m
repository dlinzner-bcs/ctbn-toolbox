%Generate synthetic data
L=40;     %number of nodes (posterior simulation of hundreds of nodes is possible in terms of minutes)
steps=10; %number of transitions in trajectorie
dt=0.005; %simulation time step
SIGMA=0.8;%variance of gaussian observation-model

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
node0=createLibOfNodes(L,A, -ones(1,L), beta,ta); %generate CTBN with glauber dynamics
steps0=10; %number of sampled states
dt0=0.1; %granularity of sample trajectories

%draw sample trajectories using gillespie sampling
[TZ,DATA0] = CTBN_TRAJ(node0,1,steps,steps0,dt0,dt);
%add noise
[Z,DATA] = corrupted_observation_gaussian(DATA0,SIGMA);

%%
%calculate approximate posterior given data
M=20; %number of ODE-iterations
%%set initial state
X0=zeros(2,L);
X0(1,:)=1/2;
X0(2,:)=1/2;
%%set initial time
t0=21;

[MU,RHO,F,~] = P_CVMCTBN_STAR_POSTV3(node0,2,dt,M,X0,t0,Z{1},TZ{1},0.01);
%%
%plot approximate posterior
k=1; %node index

m=squeeze(MU(:,k,1)-MU(:,k,2)); %approx expected post. state
var=2*(abs(squeeze(MU(:,k,1).*(1-MU(:,k,1)))));
x=linspace(0,TZ{1}(end),length(m))';

figure(2)
hold on
plot(TZ{1},DATA{1}(:,k),'o','LineWidth', 1)
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

% Plot posterior mean
plot(x,m, 'r', 'LineWidth', 2)

ylim([-2,2])
