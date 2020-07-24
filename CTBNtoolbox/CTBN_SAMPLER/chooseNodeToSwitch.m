function choosedNodeNumb=chooseNodeToSwitch(rates, generalRate)
%returns number of node that should change it`s state depending on vector
%of rates and general rate
for w=1:length(rates)
    normRates(w)=rates(w)/generalRate;
end
interval=[normRates(1)];
for w=2:length(rates)
    interval=[interval, interval(w-1)+normRates(w)];
end
randToChoose=rand(1);
for w=1:length(rates)
    if randToChoose<interval(w)
        choosedNodeNumb=w;
        break;
    end
end
end