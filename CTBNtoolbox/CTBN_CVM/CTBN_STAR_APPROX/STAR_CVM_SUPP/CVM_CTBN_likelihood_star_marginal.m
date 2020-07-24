function [F,E,H] = CVM_CTBN_likelihood_star_marginal(node,mu,rho,dt,tau,M)
%UNTITLED Summary of this function goes here
%   Detailed explanation goes here
sm=size(mu);
L=length(node);
T=sm(1);
E=0;
H=0;


%%%calculate Energy
for i=1:L
    
    [T1,M1] = CTBN_cond_stat_star(node,i,T-1,dt,mu(1:end-1,:,:),rho(1:end-1,:,:),2);
    Q_i=node(i).cellOfCondRM;
    
    if isempty(node(i).parents)
        M0=0;
        T0=0;
        M0=M1{1};
        T0=T1{1};
        
        Q=Q_i{1};
        
        for d=1:2
            for d_=1:2
                
                if d~=d_
                    E=E+M(1,i,d,d_)*log(tau(1,i,d))-(M(1,i,d,d_)+M0(d,d_))*log(tau(1,i,d)+T0(d))+(1/2)*log(M(1,i,d,d_))-(1/2)*log(M(1,i,d,d_)+M0(d,d_))+(M(1,i,d,d_)+M0(d,d_))*log(M(1,i,d,d_)+M0(d,d_))-(M(1,i,d,d_))*log(M(1,i,d,d_))-M0(d,d_);
                    E=E+M0(d,d_)*log(Q(d,d_));
                    H=H+M0(d,d_);
                else
                    E=E+T0(d)*Q(d,d);
                end
                
            end
        end
        
    else
        for k=1:length(node(i).allPosStatesOfParents)
            
            M0=0;
            T0=0;
            M0=M1{k};
            T0=T1{k};
            
            Q=Q_i{k};
            
            for d=1:2
                for d_=1:2
                    
                    if d~=d_
                        E=E+M(k,i,d,d_)*log(tau(k,i,d))-(M(k,i,d,d_)+M0(d,d_))*log(tau(k,i,d)+T0(d))+(1/2)*log(M(k,i,d,d_))-(1/2)*log(M(k,i,d,d_)+M0(d,d_))+(M(k,i,d,d_)+M0(d,d_))*log(M(k,i,d,d_)+M0(d,d_))-(M(k,i,d,d_))*log(M(k,i,d,d_))-M0(d,d_);
                        E=E+M0(d,d_)*log(Q(d,d_));
                        H=H+M0(d,d_);
                    else
                        E=E+T0(d)*Q(d,d);
                    end
                    
                end
            end
            
        end
    end
    
end

%%%calculate Entropy


for i=1:L
    
    ps=node(i).allPosStatesOfParents;
    ps(ps==-1)=2;
    ps(ps==1)=1;
    
    Q_i=node(i).cellOfCondRM;
    if (isempty(node(i).parents)==0)
        
        sp=size(ps);
        for mn=1:sp(1)
            pa=1;
            ki=1;
            
            for k=node(i).parents
                pa=pa.*mu(2:end-1,k,ps(mn,ki));
                ki=ki+1;
            end
            
            
            for d=1:2
                for d_=1:2
                    if d_~=d
                        H=H-trapz(linspace(dt,dt*T-dt,T-2),mu(2:end-1,i,d).*pa.*(rho(2:end-1,i,d_)./rho(2:end-1,i,d)).*Q_i{mn}(d,d_).*log(abs(Q_i{mn}(d,d_).*rho(2:end-1,i,d_)./rho(2:end-1,i,d))));
                    end
                end
            end
            
        end
        
    else
        
        for d=1:2
            for d_=1:2
                if d_~=d
                    H=H-trapz(linspace(dt,dt*T-dt,T-2),mu(2:end-1,i,d).*rho(2:end-1,i,d_)./rho(2:end-1,i,d).*Q_i{1}(d,d_).*(log(Q_i{1}(d,d_).*rho(2:end-1,i,d_)./rho(2:end-1,i,d))));
                end
            end
        end
    end
    
    
end

F=E+H;
