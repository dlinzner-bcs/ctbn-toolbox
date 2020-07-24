function [node,mu,rho,NO,F] = CTBNCVM_STAR_EM_MULTIPLE_TRAJ_PRIOR_POSTV2(node,DATA0,time0,dt,M,H,tau,M_T,thresh)
%UNTITLED3 Summary of this function goes here
%   Detailed explanation goes here

F=0;
mu=0;
rho=0;
PENALTY=0;
for h=1:H

    NO{h}=node; 
    [mu,rho,TK,F(h)] = CVM_STAR_EXPECTATION_STEP_MULTIPLE_TRAJ_POSTV2(node,DATA0,time0,dt,M); 
    F
    %check convergence 
    if h>5
        tail=F(end:-1:end-4);
        stdF=(mean(tail)^2-std(tail))/mean(tail)^2;
        (1-stdF)
        if (1-stdF)<thresh
            break
        end
        
    end
    
    [node] = CVM_STAR_MAXIMIZATION_STEP_MULTIPLE_TRAJ_PRIOR(node,mu,rho,TK,dt,tau,M_T);
    
end

end

