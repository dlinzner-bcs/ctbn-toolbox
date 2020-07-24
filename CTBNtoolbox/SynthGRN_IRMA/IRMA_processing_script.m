%run experiment on IRMA data
%IRMA dataset needs to be downloaded and expression data
%of individual genes need to be extracted and saved as
%'so_*genename*_*traj_number*.mat' for switch on and as
%'soff_*genename*_*traj_number*.mat' for switch off 

%%
%switch on

Z=cell(5,1);
time=cell(5,1);
genes=cell(5,1);
genes{1}='SWI5';
genes{2}='CBF1';
genes{3}='GAL4';
genes{4}='GAL80';
genes{5}='ASH1';

for n=1:length(genes)
    curr_gene=genes{n};
    for d=1:5
        name=sprintf('so_%s_%d.mat',curr_gene,d);
        A=load(name);
        H{d}=cell2mat(struct2cell(A));
    end
    DAT=[];
    for d=1:5
        DAT=[DAT; H{d}(:,6)];
    end
    m=mean(DAT);
    sm=var(DAT);
    for d=1:5
        x=H{d}(:,6);
        p_over=((x>m)-(erf((x-m)/sm)));
        p_under=((x<m)-(erf((m-x)/sm)));
        Z{d}(:,:,n)=[p_over,p_under]';
        time{d}= H{d}(:,1)/10;
    end
end

save('switch-on.mat','Z','time');
%%
%%
Z=cell(4,1);
time=cell(4,1);
genes=cell(5,1);
genes{1}='SWI5';
genes{2}='CBF1';
genes{3}='GAL4';
genes{4}='GAL80';
genes{5}='ASH1';

for n=1:length(genes)
    curr_gene=genes{n};
    for d=1:4
        name=sprintf('soff_%s_%d.mat',curr_gene,d);
        A=load(name);
        H{d}=cell2mat(struct2cell(A));
    end
    DAT=[];
    for d=1:4
        DAT=[DAT; H{d}(:,6)];
    end
    m=mean(DAT);
    sm=var(DAT);
    for d=1:4
        x=H{d}(:,6);
        p_over=((x>m)-(erf((x-m)/sm)));
        p_under=((x<m)-(erf((m-x)/sm)));
        Z{d}(:,:,n)=[p_over,p_under]';
        time{d}= H{d}(:,1)/10;
    end
end

save('switch-off.mat','Z','time');
