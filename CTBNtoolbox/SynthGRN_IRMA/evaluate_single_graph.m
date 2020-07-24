% mN is inferred network

%%define ground truth
Z=zeros(L);
Z(1,3)=1;
Z(2,1)=1;Z(2,5)=1;
Z(3,2)=1;Z(3,4)=1;
Z(4,1)=1;Z(4,3)=1;
Z(5,1)=1;

adjacency2test(mN)
adjacency2gold(Z)

%run evaluation script
%creates file with evaluation results
Dream_Evaluation_Script('BCS_DIRECTED-UNSIGNED_FiveGene_qPCR.txt', 'GoldStandard_UNDIRECTED-UNSIGNED_FiveGene_qPCR.txt');
