function [DATAC,D] = corrupted_observation_gaussian(DATA,SIGMA)

N_TRAJ=length(DATA);
sD=size(DATA{1});
L=sD(2);
    
for j=1:N_TRAJ
   d=size( DATA{j});
   tau=d(1);
   Y=zeros(2,tau,L);
   
   D{j}=DATA{j}+normrnd(0,ones(d)*SIGMA);
   Y(1,:,:)=normpdf( D{j},1,SIGMA);
   Y(2,:,:)=normpdf( D{j},-1,SIGMA);
   
   Z=round(Y,5);
   Y(Z==0)=10^(-4); %add small probability for numerical stability
 
   DATAC{j}=Y;
   
end


end

