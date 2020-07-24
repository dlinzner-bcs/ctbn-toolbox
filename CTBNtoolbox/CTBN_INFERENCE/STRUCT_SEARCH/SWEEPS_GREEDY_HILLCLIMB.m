function [NET,T_SCORES,node] = SWEEPS_GREEDY_HILLCLIMB(MAX_PAR,L,node,DATAC,time0,dt,M,H,tau,M_T,thresh)
%PERFORM GREEDY HILL CLIMBING

for h=1:L
    for Kp=1:MAX_PAR+1
        p_comb{Kp}=K_Parents(h,L,Kp-1);
    end
    j=0;
    for Kp=1:MAX_PAR+1
        sp=size(p_comb{Kp});
        for i=1:sp(1)
            j=j+1;
            COMB{j}=p_comb{Kp}(i,:);
        end
    end
    COMB{j+1}=[]; %COMB is possible set of families of node h
    SCORES=zeros(1,length(COMB)); 
    parfor i=1:length(COMB)
        [node1]=shift_parent(node,h,COMB{i});
        node2net(node1)
        try
            [node1]=createLibOfNodesGammaPrior(L,node2net(node1), -ones(1,L),tau,M_T);
            [node1,mu,rho,~,~] =  CTBNCVM_STAR_EM_MULTIPLE_TRAJ_PRIOR_POSTV2(node1,DATAC,time0,dt,M,H,tau,M_T,thesh);
        catch
            [node1]=createLibOfNodesGammaPrior(L,node2net(node1), -ones(1,L),tau,M_T);
            [node1,mu,rho,~,~] =  CTBNCVM_STAR_EM_MULTIPLE_TRAJ_PRIOR_POSTV2(node1,DATAC,time0,dt,M,H,tau,M_T,thresh);
        end
        Fmax=0;
        for k=1:length(mu)
            [F_k,~,~] = CVM_CTBN_likelihood_star_marginal(node1,mu{k},rho{k},dt,tau,M_T);
            [lnP]=expected_loglikelihood(mu{k},dt,DATAC{k},time0{k});
            Fmax=Fmax+F_k+lnP;
        end
        SCORES(i)=Fmax;
        NET{h,i}=node2net(node1); 
        T_SCORES{h,i}=SCORES(i);
    end
    SCORES
    A=node2net(node);
    [maxi, max_ind]=max(SCORES);
    A=NET{h,max_ind}; %extract best scoring net
    [node]=createLibOfNodesGammaPrior(L,A, -ones(1,L),tau,M_T);
    A
end

end

