function [psi0] = P_CVMCTBN_COUP_V2(node,i,mu,rho,D)
%CALCULATE CHILD RESPONSE FOR NODE i
%AT TIME t IN ITERATION m

sm=size(mu);
tau=sm(1);

psi0=zeros(tau,D); % time dimension
%%%%% diamond term
for d=1:D
    
    %%%%AVERAGING RATES OF CHILDREN%%%%%
    
    for j=node(i).children
        %%%%%REMOVE CONSIDERED PARENT FROM PARENT SET NODE i AS STATE IS FIXED
        n=node(j).parents;
        n(n==i)=[];
        %%%%%STATE TRANSFORM ALWAYS NECESSARY
        ps=node(j).allPosStatesOfParents;
        ps(ps==-1)=2;
        ps(ps==1)=1;
        
        %%%%%CONDITION ON STATE OF NODE i BY REMOVING ALL NON CONSITENT STATES
        sp=size(ps);
        pk=[1:sp(1)];
        pk(ps(:,i==node(j).parents)~=d)=[];
        ps(ps(:,i==node(j).parents)~=d,:)=[];
        
        ps(:,i==node(j).parents)=[];
        
        %%%%%LOAD CHILDRENS RATE MATRIXES
        Q_i=node(j).cellOfCondRM;
        Q_avg=zeros(tau,2);
        sp=size(ps);
        
        %%%%%CALC. EXPECTED RATE MATRIXES
        for d_=1:D
            for mn=1:sp(1)
                pa=ones(tau,1);
                ki=1;
                
                for k=n
                    pa=pa.*mu(:,k,ps(mn,ki));
                    ki=ki+1;
                end
                
                Q_avg(:,d_)=Q_avg(:,d_)+Q_i{pk(mn)}(d_,3-d_).*pa;
                
            end
            
            if (isempty(node(j).parents))
                
                Q_avg(:,d_)=Q_i{1}(d_,3-d_).*ones(tau,1);
                
            end
        end
        %%%%%CALC. CHILD RESPONSE
        C=zeros(tau,1);
        
        for d_=1:D
            
            C=C+mu(:,j,d_).*(Q_avg(:,d_).*rho(:,j,3-d_)./rho(:,j,d_)-Q_avg(:,d_));
            
        end
        psi0(1:tau-1,d)=psi0(1:tau-1,d)+C(1:tau-1);%C(tau:-1:1);
    end
end

end

