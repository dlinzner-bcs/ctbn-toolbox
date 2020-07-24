function [q,logq] = P_CVMCTBN_EFF_RATES_V2(node,i,mu,D)
%CALCULATE EFFECTIVE RATES FOR NODE i GIVEN THE NEIGHBOURS
%AT TIME t IN ITERATION m

sm=size(mu);
tau=sm(1);

q=zeros(tau,D);
logq=zeros(tau,D);

%%%get effective rate
ps=node(i).allPosStatesOfParents;
ps(ps==-1)=2;
ps(ps==1)=1;

Q_avg=zeros(tau,D);
logQ_avg=zeros(tau,D);
Q_i=node(i).cellOfCondRM;
%%%%averaging
for d=1:D
    if (isempty(node(i).parents)==0)
        sp=size(ps);
        for mn=1:sp(1)
            pa=ones(tau,1);
            ki=1;
            
            for k=node(i).parents
                pa=pa.*mu(:,k,ps(mn,ki));
                ki=ki+1;
            end

            Q_avg(:,d)=Q_avg(:,d)+Q_i{mn}(d,3-d).*pa;
            logQ_avg(:,d)=logQ_avg(:,d)+log((Q_i{mn}(d,3-d))).*pa;
        end
        
        %%%%%%%
        
    else
        Q_avg(:,d)=Q_i{1}(d,3-d).*ones(tau,1);
        logQ_avg(:,d)=log((Q_i{1}(d,3-d)).*ones(tau,1));
    end
end


for d=1:D
    
    q=Q_avg;
    logq=logQ_avg;
    
end

end

