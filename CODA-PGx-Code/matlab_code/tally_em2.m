function tally_em2(lof_gof_km,lof_gof_hits,hit_anti_hit,dd,fid,fid_str)

[aa.gene,aa.hit,~,aa.km_lof,~,~,aa.km_p,aa.J]=textread(fid_str,'%s%s%n%s%s%s%n%n','delimiter','\t');
aa.driver=vlookup(aa.gene,dd.gene,dd.driver);

% rows to keep
keep1=[];
keep2=[];
keep3=[];

% keep km_lof/gof/mut_only
for nn=1:length(lof_gof_km) % column 4
    keep1=[keep1,strcmp(lof_gof_km(nn),aa.km_lof)];
end
keep1=max(keep1,[],2);

% keep hits/antihits
for nn=1:length(hit_anti_hit) % column 2
    keep2=[keep2,strcmp(hit_anti_hit(nn),aa.hit)];
end
keep2=max(keep2,[],2);

% lof/gof drivers
for nn=1:length(lof_gof_hits) % from aa.driver
    keep3=[keep3,aa.driver==lof_gof_hits(nn)];
end
keep3=max(keep3,[],2);

keep=min([keep1,keep2,keep3],[],2);

% only keep hits where km_lof agrees with hit_lof
% aa.driver == -1/0/1
% aa.km_lof == gof/mut_only/lof
agree = ( strcmp(aa.km_lof,'LOF')|strcmp(aa.km_lof,'MUT_ONLY') )&( aa.driver==-1| aa.driver==0 )| strcmp(aa.km_lof,'GOF') & aa.driver==1;
keep(~agree)=0;

aa=downSize(aa,logical(keep));

% tally-em-up
total=sum(keep);
in=sum( (strcmp(aa.hit,'HIT')&aa.J<0) | (strcmp(aa.hit,'ANTIHIT')&aa.J>0) );


for ii=1:length(lof_gof_km)-1
    fprintf(fid,[lof_gof_km{ii} ',']);
end
fprintf(fid,[lof_gof_km{end} '\t']);

for ii=1:length(lof_gof_hits)-1
    fprintf(fid,[num2str(lof_gof_hits(ii)) ',']);
end
fprintf(fid,[num2str(lof_gof_hits(end)) '\t']);

for ii=1:length(hit_anti_hit)-1
    fprintf(fid,[hit_anti_hit{ii} ',']);
end
fprintf(fid,[hit_anti_hit{end} '\t']);

if in==0
    total=1;
end

fprintf(fid,[num2str(in) '\t' num2str(total) '\t' num2str(in/total) '\t' fid_str '\n']);
    












