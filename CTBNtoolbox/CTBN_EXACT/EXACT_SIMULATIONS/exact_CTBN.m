function [P] = exact_CTBN(L,Q,e0,tau,dt)
%UNTITLED6 Summary of this function goes here
%   Detailed explanation goes here

Pi=zeros(2^L,1);
%Pi(statePos(double(e0(1,:)),2,L))=1;
Pi=e0;
P=zeros(tau,2^L);
P(1,:)=Pi;
A=expm(dt*Q);
for t=1:tau-1
    P(t+1,:)=P(t,:)*A;
end

end

