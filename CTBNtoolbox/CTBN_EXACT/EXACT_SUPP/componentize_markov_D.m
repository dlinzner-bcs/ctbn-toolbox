function [mu_m,rho_m,gam_m] = componentize_markov_D(tau,L,D,rho,mu,gam)
%global

s=de2bi(0:(D^L)-1);


for i=1:L
    for d=1:D
        loc_state(i,d-1,D,L)
        mu_m(1:tau,i,d)= sum(mu(1:tau,loc_state(i,d-1,D,L)),D);
        rho_m(1:tau,i,d)= sum(rho(1:tau,loc_state(i,d-1,D,L)),D);
        
        for d_=1:D
            
            A=0;
            if (d_~=d)
                
                x=loc_state(i,d-1,D,L);
                y=loc_state(i,d_-1,D,L);
                
                for k=1:length(x)
                    for h=1:length(y)
                        if sum(abs(s(x(k),:)-s(y(h),:)))==1
                            
                            A=A+ gam(1:tau,x(k),y(h));
                            
                        end
                    end
                end
                
                gam_m(:,i,d_,d_)= A;
                
            end
        end
        
    end
end



end

