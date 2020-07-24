function [time0,DATA0] = CTBN_TRAJ(node,N_TRAJ,steps,steps0,dt0,dt)

L=length(node);

for i=1:N_TRAJ
    tr=0;
    t0=0;
    
    st=real(2*(rand(L,1)>0.5)-ones(L,1));
    
    while (tr<=2*dt)||(t0==0)
        for j=1:L
            node(j).state=st(j);
        end
        [time0{i},DATA0{i}] = CTBN_Trajectories_Incomplete(node,steps,steps0,dt0);
        tr=min(diff(time0{i}));
        t0=(min(time0{i}));
        
    end
   
end

end

