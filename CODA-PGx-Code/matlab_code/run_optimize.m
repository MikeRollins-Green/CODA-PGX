function run_optimize(signi,dir_out)

fclose all
delete [dir_out './*'];
rng shuffle

addpath('C:\Users\Tomas\Dropbox\co\science\pipe\matlab_code')
load C:\Users\Tomas\Dropbox\co\science\pipe\survival_package % lof = cnv<0, gg_taste>0; gg = gg_taste>0 (mutations only); gof = cnv>0 & gg_taste>0;
[look.mydrug,look.tcga_drug]=textread('C:\Users\Tomas\Dropbox\co\science\pipe\mydrug2tcga_drug.txt','%s%s','delimiter','\t');
[dd.gene,dd.driver]=textread('C:\Users\Tomas\Dropbox\celgene_backup\H\tcga\drivers.txt','%s%n','delimiter','\t');

nn=1;

% one
trt_line=1;
max_hit_p=0.1;
min_fc=0.5;
min_pts=30;
max_km_p=0.2;
min_J=1;
fid_str=signi2survival_pipe(signi,clin,drugs,gof,lof,gg,min_fc,max_hit_p,min_pts,trt_line,max_km_p,look,min_J);
optimize_em(dd,[dir_out '/out.' num2str(nn) '.tab'],fid_str);
nn=nn+1;

% two
trt_line=1;
max_hit_p=0.1;
min_fc=1;
min_pts=30;
max_km_p=0.2;
min_J=1;
fid_str=signi2survival_pipe(signi,clin,drugs,gof,lof,gg,min_fc,max_hit_p,min_pts,trt_line,max_km_p,look,min_J);
optimize_em(dd,[dir_out '/out.' num2str(nn) '.tab'],fid_str);
nn=nn+1;

% three
trt_line=1;
max_hit_p=0.1;
min_fc=0.5;
min_pts=30;
max_km_p=0.1;
min_J=1;
fid_str=signi2survival_pipe(signi,clin,drugs,gof,lof,gg,min_fc,max_hit_p,min_pts,trt_line,max_km_p,look,min_J);
optimize_em(dd,[dir_out '/out.' num2str(nn) '.tab'],fid_str);
nn=nn+1;

% four
trt_line=1;
max_hit_p=0.05;
min_fc=0.5;
min_pts=30;
max_km_p=0.2;
min_J=1;
fid_str=signi2survival_pipe(signi,clin,drugs,gof,lof,gg,min_fc,max_hit_p,min_pts,trt_line,max_km_p,look,min_J);
optimize_em(dd,[dir_out '/out.' num2str(nn) '.tab'],fid_str);
nn=nn+1;

% five
trt_line=1;
max_hit_p=0.05;
min_fc=0.5;
min_pts=30;
max_km_p=0.1;
min_J=1;
fid_str=signi2survival_pipe(signi,clin,drugs,gof,lof,gg,min_fc,max_hit_p,min_pts,trt_line,max_km_p,look,min_J);
optimize_em(dd,[dir_out '/out.' num2str(nn) '.tab'],fid_str);
nn=nn+1;

% six
trt_line=1;
max_hit_p=0.1;
min_fc=0.5;
min_pts=30;
max_km_p=0.1;
min_J=2;
fid_str=signi2survival_pipe(signi,clin,drugs,gof,lof,gg,min_fc,max_hit_p,min_pts,trt_line,max_km_p,look,min_J);
optimize_em(dd,[dir_out '/out.' num2str(nn) '.tab'],fid_str);
nn=nn+1;

% seven
trt_line=1;
max_hit_p=0.1;
min_fc=1;
min_pts=30;
max_km_p=0.2;
min_J=2;
fid_str=signi2survival_pipe(signi,clin,drugs,gof,lof,gg,min_fc,max_hit_p,min_pts,trt_line,max_km_p,look,min_J);
optimize_em(dd,[dir_out '/out.' num2str(nn) '.tab'],fid_str);
nn=nn+1;

% eight
trt_line=1;
max_hit_p=0.05;
min_fc=0.5;
min_pts=30;
max_km_p=0.2;
min_J=2;
fid_str=signi2survival_pipe(signi,clin,drugs,gof,lof,gg,min_fc,max_hit_p,min_pts,trt_line,max_km_p,look,min_J);
optimize_em(dd,[dir_out '/out.' num2str(nn) '.tab'],fid_str);
nn=nn+1;

% nine
trt_line=2;
max_hit_p=0.1;
min_fc=0.5;
min_pts=30;
max_km_p=0.2;
min_J=3;
fid_str=signi2survival_pipe(signi,clin,drugs,gof,lof,gg,min_fc,max_hit_p,min_pts,trt_line,max_km_p,look,min_J);
optimize_em(dd,[dir_out '/out.' num2str(nn) '.tab'],fid_str);
nn=nn+1;

% ten
trt_line=2;
max_hit_p=0.01;
min_fc=0.5;
min_pts=30;
max_km_p=0.2;
min_J=1;
fid_str=signi2survival_pipe(signi,clin,drugs,gof,lof,gg,min_fc,max_hit_p,min_pts,trt_line,max_km_p,look,min_J);
optimize_em(dd,[dir_out '/out.' num2str(nn) '.tab'],fid_str);
nn=nn+1;

% eleven
trt_line=2;
max_hit_p=0.1;
min_fc=0.5;
min_pts=30;
max_km_p=0.2;
min_J=1;
fid_str=signi2survival_pipe(signi,clin,drugs,gof,lof,gg,min_fc,max_hit_p,min_pts,trt_line,max_km_p,look,min_J);
optimize_em(dd,[dir_out '/out.' num2str(nn) '.tab'],fid_str);
nn=nn+1;
