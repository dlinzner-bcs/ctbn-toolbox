function [F,E,H] = CTBN_likeihood_global(mu,gam,Q,dt)
%UNTITLED9 Summary of this function goes here
%   Detailed explanation goes here

sm=size(mu);
D=sm(2);

E=zeros(sm(1)-2,1);
H=zeros(sm(1)-2,1);

gam=permute(gam,[1 2 3]);

for t=2:sm(1)-1
    em=0;
    eg=0;
    hm=0;
    hg=0;
    hg0=0;
    
    for d=1:D
        em=em+mu(t,d)*Q(d,d);
        for d_=1:D
            if (d~=d_)
                if gam(t,d,d_)~=0
                    eg=eg+gam(t,d,d_)*log(Q(d,d_));
                    hg0=hg0+gam(t,d,d_)*log(gam(t,d,d_));
                    
                    hm=hm+gam(t,d,d_)*log(mu(t,d));
                    
                    hg=hg+gam(t,d,d_);
                end
            end
            
        end
        
    end

    E(t-1)=em+eg;
    H(t-1)=hg+hm-hg0;
end
f=H+E;
xt=linspace(dt,dt*(sm(1)-1),sm(1)-2);

F=trapz(xt,f');
