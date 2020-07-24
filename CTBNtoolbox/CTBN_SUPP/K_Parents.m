function [m] = K_Parents(l,L,K)
%Returns all possible parent positions of K parents in L nodes for node l
%needs function uniqueperms


if K>0
    
    s=zeros(L-1,1);
    s(1:K)=1;
    
    n=[1:L];
    n(l)=[];
    
    states=uniqueperms(s);
    ss=size(states);
    m=zeros(ss(1),K);
    for i=1:ss(1)
        m(i,:)=n(states(i,:)==1);
    end
    
else
    
    m=[];
    
end

% if K>0
%     
%     s=zeros(L-1,1);
%     s(1:K)=1;
%     
%     n=[1:L];
%     n(l)=[];
%     
%     states=uniqueperms(s);
%     for i=1:length(states)
%         m(i,:)=n(states(i,:)==1);
%     end
%     
% else
%     
%     m=[];
%     
% end

%end

% s=zeros(L-1,1);
% s(1:K)=1;
%
% n=[1:L];
% n(l)=[];
%
% states=uniqueperms(s);
% for i=1:length(states)
%     m(i,:)=n(states(i,:)==1);
% end
%
% end

