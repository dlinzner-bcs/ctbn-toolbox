function edges = generateEdges(numbOfV, probEdgeExist)
% Generation of random adjency matrix
edges = sign(probEdgeExist - rand(numbOfV, numbOfV));
for i=1:numbOfV
    edges(i,i)=-1;
end
end