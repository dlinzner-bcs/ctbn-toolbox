function [ mu_m ] = componentize_mu(mu,L)
%UNTITLED7 Summary of this function goes here
%   Detailed explanation goes here

for i=1:L
    for d=1:2
        mu_m(:,i,d)= sum(mu(:,loc_state(i,d-1,2,L)),2);
    end
end

end

