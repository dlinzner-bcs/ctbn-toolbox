function [F,E,H] = CTBN_likeihood_local_star(node,mu,rho,dt)
%UNTITLED9 Summary of this function goes here
%   Detailed explanation goes here

sm=size(mu);
L=sm(2);
D=sm(3);

E=zeros(sm(1)-2,1);
H=zeros(sm(1)-2,1);
for t=2:sm(1)-1
    
    em=0;
    eg=0;
    hm=0;
    hg=0;
    hg0=0;
    
    for i=1:L
        
        %%%%%% calculate log(q)q
        ps=node(i).allPosStatesOfParents;
        ps(ps==-1)=2;
        ps(ps==1)=1;
        Q_avg=zeros(2,2);
        Q_avglog=zeros(2,2);
        
        Q_i=node(i).cellOfCondRM;
        if (isempty(node(i).parents)==0)
            
            sp=size(ps);
            for mn=1:sp(1)
                pa=1;
                ki=1;
                
                for k=node(i).parents
                    pa=pa*mu(t,k,ps(mn,ki));
                    ki=ki+1;
                end
                
                Q_avglog=Q_avglog+log(abs(Q_i{mn})).*Q_i{mn}*pa;
                Q_avg=Q_avg+Q_i{mn}*pa;
            end
            
        else
            
            Q_avglog=log(abs(Q_i{1})).*Q_i{1};
            Q_avg=Q_i{1};
            
        end
        
        
        %%%%%
        
        for d=1:D
            
            em=em+mu(t,i,d)*Q_avg(d,d);
            
            for d_=1:D
                if (d~=d_)
                    eg=eg+mu(t,i,d)*Q_avglog(d,d_)*rho(t,i,d_)/rho(t,i,d);
                    hg0=hg0+mu(t,i,d)*rho(t,i,d_)/rho(t,i,d)*(Q_avglog(d,d_)+Q_avg(d,d_)*log(rho(t,i,d_)/rho(t,i,d)));
                    hg=hg+mu(t,i,d)*rho(t,i,d_)*Q_avg(d,d_)/rho(t,i,d);
                end
            end
        end
    end
    
    E(t-1)=em+eg;
    H(t-1)=hg-hg0;
    
    
end
f=H+E;
xt=linspace(dt,dt*(sm(1)-1),sm(1)-2);

F=trapz(xt,f');
