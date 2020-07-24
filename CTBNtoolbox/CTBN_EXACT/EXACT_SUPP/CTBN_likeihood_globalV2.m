function [F,E,H] = CTBN_likeihood_globalV2(mu,gam,Q,dt)
%UNTITLED9 Summary of this function goes here
%   Detailed explanation goes here

sm=size(mu);
D=sm(2);
T=sm(1);
E=0;
H=0;
F=0;

gam=permute(gam,[1 3 2]);
for d=1:D
    E=E+trapz(linspace(0,T*dt,T),mu(:,d)).*Q(d,d);
    for d_=1:D
        if (d~=d_)
            
            E=E+trapz(linspace(0,T*dt,T),gam(:,d,d_)).*log(Q(d,d_));
            H=H-trapz(linspace(0,T*dt,T),gam(:,d,d_).*log(gam(:,d,d_)));
            H=H+trapz(linspace(0,T*dt,T),gam(:,d,d_).*log(mu(:,d)));
            
            H=H+trapz(linspace(0,T*dt,T),gam(:,d,d_));
            
        end
        
    end
    
end

F=H+E;

end

