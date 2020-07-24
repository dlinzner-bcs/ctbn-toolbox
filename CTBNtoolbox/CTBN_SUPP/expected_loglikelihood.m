function [F] = expected_loglikelihood(mu,dt,Z,TZ)
%UNTITLED18 Summary of this function goes here
%   Detailed explanation goes here
sz=size(Z);
D=sz(1);
tau=sz(2);
L=sz(3);

F=0;

for j=1:tau-1
    for i=1:sz(3)
        for d=1:D
            
            F=F-Z(d,j,i)*mu(floor(TZ(j)/dt),i,d);
            
        end
    end
end

end

