function [T,M] = CTBN_cond_stat_star(node,i,t,dt,mu,rho,D)
%CALCULATE COND STATISTICS OF CTBN

%%%get effective rate
ps=node(i).allPosStatesOfParents;
ps(ps==-1)=2;
ps(ps==1)=1;
Q_i=node(i).cellOfCondRM;

sp=size(ps);

T=cell(sp(1),1);
M=cell(sp(1),1);

if (isempty(node(i).parents)==0)
    %%%%averaging
    
    sp=size(ps);
    for mn=1:sp(1)
        T_i=0;
        M_i=0;
        
        pa=ones(t,1);
        ki=1;
       
        for k=node(i).parents
            pa=pa.*mu(1:t,k,ps(mn,ki));  
              
            ki=ki+1;
       
        end
      
        for d=1:D
            
           T_i(d)=trapz(linspace(0,dt*t,t),mu(1:t,i,d).*pa);
           
           for d_=1:D
               if (d~=d_)
                 M_i(d,d_)=trapz(linspace(0,dt*t,t),pa.*mu(1:t,i,d).*rho(1:t,i,d_)./rho(1:t,i,d)).*Q_i{mn}(d,d_);
               end
           end
        end
        
        T{mn}=T_i;
        M{mn}=M_i;
     
    end
    
else
    
    for d=1:D
         T_i(d)=trapz(linspace(0,dt*t,t),mu(1:t,i,d));
         
        for d_=1:D
             if (d~=d_)
                M_i(d,d_)=trapz(linspace(0,dt*t,t),mu(1:t,i,d).*(rho(1:t,i,d_)./rho(1:t,i,d))).*Q_i{1}(d,d_);
                
             end
        end
         
    end
    
     T{1}=T_i; 
     M{1}=M_i; 
end

end

