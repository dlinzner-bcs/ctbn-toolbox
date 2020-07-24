function [Q] = jointIM_amalgation(node)
%UNTITLED9 Summary of this function goes here
%   Detailed explanation goes here
L=length(node);
D=2;

Q=zeros(D^L,D^L);
s=de2bi(0:(D^L)-1);

sz=size(s);
for i=1:sz(1)
    for j=1:sz(1)
        
        a=s(i,:);
        b=s(j,:);
        
        if sum(abs(a-b))~=1
            
            Q(i,j)=0;
            
        else
            
            ind=find(abs(a-b)==1);
            
            pa=node(ind).parents;
            if isempty(pa)==0
                mn=bi2de(a(pa))+1;
                
            else
                mn=1;
            end
            d=a(ind)+1;
            d_=b(ind)+1;
            
            
            w=node(ind).cellOfCondRM{mn}(d,d_);
            
            Q(i,j)=w;
            
        end
    end
end

for i=1:sz(1)
    Q(i,i)=-sum(Q(i,:));
end

end

