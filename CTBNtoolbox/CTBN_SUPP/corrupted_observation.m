function [DATAC] = corrupted_observation(DATA,SIGMA)

N_TRAJ=length(DATA);

for j=1:N_TRAJ
   d=size( DATA{j});
   D{j}=DATA{j}+normrnd(0,ones(d)*SIGMA);
   DATAC{j}(1,j,:)=normpdf( D{j},1,SIGMA);
   DATAC{j}(2,j,:)=normpdf( D{j},-1,SIGMA);
end


end

