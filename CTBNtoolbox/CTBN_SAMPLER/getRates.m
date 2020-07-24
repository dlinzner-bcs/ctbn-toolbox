function [rates, generalRate]=getRates(numbOfV,node)
%for each node find its current conditional rate matrix
%depending on its parent`s states
for w=1:numbOfV
    if isempty(node(w).parents)   
        node(w).currentCondRM=node(w).cellOfCondRM{1};
    else       
        parentsState = checkParentsState(node, w);
        choosedRateM= chooseConditionalRM(parentsState, node, w);
        node(w).currentCondRM=choosedRateM;
    end
end

%vector and sum of rates depending on state
rates=[];
for w=1:numbOfV
    if node(w).state==+1
        rates=[ rates, node(w).currentCondRM(1,2)];
    end
    if node(w).state==-1
        rates=[ rates, node(w).currentCondRM(2,1)];
    end
end
generalRate=sum(rates);
end