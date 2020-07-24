function [mu_m,rho_m,gam_m,gam_naive,gam_star] = componentize_markov(node,tau,L,D,rho,mu,gam)
%global

s=de2bi(0:(D^L)-1);

for t=1:tau
    for i=1:L
        for d=1:D

            mu_m(t,i,3-d)= sum(mu(t,loc_state(i,d-1,D,L)));
            rho_m(t,i,3-d)= sum(rho(t,loc_state(i,d-1,D,L)));       
            
            for d_=1:D
                
                A=0;
                if (d_~=d)
                    
                    x=loc_state(i,d-1,D,L);
                    y=loc_state(i,d_-1,D,L);
                    
                    for k=1:length(x)
                        for h=1:length(y)
                            if sum(abs(s(x(k),:)-s(y(h),:)))==1
                          
                                A=A+ gam(t,x(k),y(h));
                                
                            end
                        end
                    end
                    
                    gam_m(t,i,3-d,3-d_)= A;
                    
                end
            end
            
        end
    end
end

gam_star=zeros(L,D,D);
gam_naive=gam_star;
q=zeros(L,D);
logq=q;


for t=1:tau
    for i=1:L
        
        [q(i,:),logq(i,:)] = P_CVMCTBN_EFF_RATES(node,i,t,mu_m,D);
        
        for d=1:D
            for d_=1:D
                if (d~=d_)
                    
                    gam_star(t,i,d,d_)=mu_m(t,i,d)*q(i,d)*rho_m(t,i,d_)/rho_m(t,i,d_);
                    gam_naive(t,i,d,d_)=mu_m(t,i,d)*exp(logq(i,d))*rho_m(t,i,d_)/rho_m(t,i,d_);
                    
                end
                
            end
        end
        
    end
end


end

