function dydt = CVM_CTBN_rho_fast(D,t, y, ft, f, gt, g)

%%% AS ODE IS TIME DEPENDENT ON EFF. RATES WE HAVE TO INTERPOLATE

ind=find(ft>=t,1); %take current time-point fast - but only valid for small timesteps
for d=1:D
    f(d)=f(ind,d);
    g(d)=g(ind,d);
end

M=[-f(1)-g(2) f(1); f(2) -f(2)-g(1)];
dydt = -M*y; % Evalute ODE at times t
