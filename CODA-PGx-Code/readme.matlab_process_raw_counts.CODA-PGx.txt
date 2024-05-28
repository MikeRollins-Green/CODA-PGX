
%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%
%
% This is a suite of matlab functions that expects as input a matrix of guides x samples, where the entries in the matrix are NGS sequencing read counts
% see sample inputs in folder "sample_inputs"
% 
%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%


%% 

%
%
%  _                 _ 
% | |               | |
% | | ___   __ _  __| |
% | |/ _ \ / _` |/ _` |
% | | (_) | (_| | (_| |
% |_|\___/ \__,_|\__,_|
%                      
%                      
%
% load
addpath('matlab_code');
cd 'my_test_dir'
[pos.original_sort,pos.gwz_id,pos.collabels,pos.operator,pos.cell_line,pos.pt_pool,pos.drug,pos.conc,pos.group]=textread('sample_inputs/collabels.txt','%n%s%s%s%s%n%s%n%n','delimiter','\t');
[pos.rowlabels,pos.gene,pos.pool]=textread('sample_inputs/rowlabels.txt','%s%s%n','delimiter','\t');
pos.counts=load('sample_inputs/matrix.tab');
pos.counts(isnan(pos.counts))=0;
pos.cpm=to_cpm(pos.counts);
pos.pool(isnan(pos.pool))=6;
pos.data=quantilenorm(pos.cpm);
save pos pos

%
%
%   ____   _____                             _   _                 _       _     _                       _   _   _                  ___  
%  / __ \ / ____|                           | | | |               (_)     | |   | |                     | | | | | |                |__ \ 
% | |  | | |       ______    __ _ _ __ ___  | |_| |__   ___   _ __ _  __ _| |__ | |_   _ __   ___   ___ | | | |_| |__   ___ _ __ ___  ) |
% | |  | | |      |______|  / _` | '__/ _ \ | __| '_ \ / _ \ | '__| |/ _` | '_ \| __| | '_ \ / _ \ / _ \| | | __| '_ \ / _ \ '__/ _ \/ / 
% | |__| | |____           | (_| | | |  __/ | |_| | | |  __/ | |  | | (_| | | | | |_  | |_) | (_) | (_) | | | |_| | | |  __/ | |  __/_|  
%  \___\_\\_____|           \__,_|_|  \___|  \__|_| |_|\___| |_|  |_|\__, |_| |_|\__| | .__/ \___/ \___/|_|  \__|_| |_|\___|_|  \___(_)  
%                                                                     __/ |           | |                                                
%                                                                    |___/            |_|                                                
%
% check how things cluster
cd 'my_test_dir'
load pos

% visual qc
tt=pos;
tt.rowlabels=tt.pool;
%tt.collabels=tt.column_pur;
tt.data=asinh(tt.data);
close all; [mc2,xind,yind]=myCluster(tt,0,0,1,1,'euclidean');
% pools look good

% dump low-coverage guides and see if impacts clustering
x=downSize(x,min(x.cpm,[],2)>200);
x=downSize(x,std(x.data')>0.2);
close all; [mc2,xind,yind]=myCluster(x,0,0,1,1,'correlation'); colorbar; colormap by

% keep all DMSO?
tt=pos;
tt.collabels=my_num2cell(tt.group)
my_pca(tt)
%% corr-all
corr_all(tt)

% distribution of guides
plot(sort(pos.counts),'linewidth',2)
set(gca,'yscale','log')
grid on
xlabel('810 expected guides sorted on coverage')
ylabel('# times sequenced')


%
%
%            _   _           
%           | | (_)          
%  _ __ __ _| |_ _  ___  ___ 
% | '__/ _` | __| |/ _ \/ __|
% | | | (_| | |_| | (_) \__ \
% |_|  \__,_|\__|_|\___/|___/
%                            
%                            
%
%
cd 'my_test_dir'
load pos

% min reads in DMSO control to be counted == 100
pos.data=quantilenorm(pos.cpm);
%pos2.data=pos2.cpm;   % sometimes see better results
posr=norm_em2(pos,50); % normalize within the same group
posdmso=norm_em_dmsoXdmso2(pos,50);
save posr posr
save posdmso posdmso

% calc_intra_sig(posr,filter_on_expr,ignore_bad_guides)
cd 'my_test_dir'
load posdmso
load posr

signi=calc_intra_sig2(posr,1,0)
signi_low_conf=calc_intra_sig2(posr,1,-1)
signi_dmso=calc_intra_sig2(posdmso,1,0)
save signi signi
save signi_dmso signi_dmso
plot_volcanoes(signi_dmso,'volcanoes_dmso')
plot_volcanoes(signi,'volcanoes')
plot_volcanoes(signi_low_conf,'volcanoes_low_conf')


%% view raw data for individual samples (QC/sanity check)
cd 'my_test_dir'
addpath('C:\Users\TomasBabak\Dropbox (personal)\co\science\pipe\matlab_code');
load pos
load posr
view_raw44('STAG2','C93',posr,pos,0)



             _       _         _       _        
            (_)     | |       | |     | |       
  _ __  _ __ _ _ __ | |_    __| | __ _| |_ __ _ 
 | '_ \| '__| | '_ \| __|  / _` |/ _` | __/ _` |
 | |_) | |  | | | | | |_  | (_| | (_| | || (_| |
 | .__/|_|  |_|_| |_|\__|  \__,_|\__,_|\__\__,_|
 | |                                            
 |_|                                            

cd 'my_test_dir'
addpath('C:\Users\TomasBabak\Dropbox (personal)\co\science\pipe\matlab_code');
load signi
mkdir outputs
x.data=signi.fc;
x.rowlabels=signi.gene;
x.collabels=signi.collabels;
print_excel_sheet(x,'outputs/fc.tab');
x.data=signi.p;
print_excel_sheet(x,'outputs/p.tab');

