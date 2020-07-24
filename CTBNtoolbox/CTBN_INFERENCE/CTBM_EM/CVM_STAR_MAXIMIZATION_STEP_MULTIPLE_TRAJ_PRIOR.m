function [node] = CVM_STAR_MAXIMIZATION_STEP_MULTIPLE_TRAJ_PRIOR(node,mu,rho,TK,dt,tau,M_T)
%UNTITLED Summary of this function goes here
%   Detailed explanation goes here

L=length(node);
NO_TRAJ=length(TK);

for i=1:L
    QK=cell(1,length(node(i).allPosStatesOfParents));
    if isempty(node(i).allPosStatesOfParents)
        M0=0;
        T0=0;
        for j=1:NO_TRAJ
            [T1,M1] = CTBN_cond_stat_star(node,i,length(TK{j})-2,dt,squeeze(mu{j}(2:end-1,:,:)),squeeze(rho{j}(2:end-1,:,:)),2);
            M0=M0+M1{1};
            T0=T0+T1{1};
        end
        QK{1}=[-(M0(1,2)+M_T(1,i,1,2))/(T0(1)+tau(1,i,1)) (M0(1,2)+M_T(1,i,1,2))/(T0(1)+tau(1,i,1));(M0(2,1)+M_T(1,i,2,1))/(T0(2)+tau(1,i,2)) -(M0(2,1)+M_T(1,i,2,1))/(T0(2)+tau(1,i,2))];
    else
        for k=1:length(node(i).allPosStatesOfParents)
            M0=0;
            T0=0;
            for j=1:NO_TRAJ
                [T1,M1] = CTBN_cond_stat_star(node,i,length(TK{j})-2,dt,squeeze(mu{j}(2:end-1,:,:)),squeeze(rho{j}(2:end-1,:,:)),2);
                
                M0=M0+M1{k};
                T0=T0+T1{k};
            end
            QK{k}=[-(M0(1,2)+M_T(k,i,1,2))/(T0(1)+tau(k,i,1)) (M0(1,2)+M_T(k,i,1,2))/(T0(1)+tau(k,i,1));(M0(2,1)+M_T(k,i,2,1))/(T0(2)+tau(k,i,2)) -(M0(2,1)+M_T(k,i,2,1))/(T0(2)+tau(k,i,2))];
        end
    end
    node(i).cellOfCondRM=QK;
end

end

