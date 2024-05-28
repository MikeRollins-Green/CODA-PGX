function optimize_em2(dd,my_file,fid_str)

fid=fopen(my_file,'w');

km_lof_gof=[{'LOF'};{'GOF'};{'MUT_ONLY'}];
hit_antihit=[{'HIT'};{'ANTIHIT'}];
drivers=[-1;0;1];

tr=[1,1,1;1,1,0;1,0,1;0,1,1;1,0,0;0,1,0;0,0,1];
[ss3,~]=size(tr);
dbl=[1,1;1,0];
[ss2,~]=size(dbl);

tr2=[1,1,1;1,1,0;1,0,1;1,0,0];
[ss4,~]=size(tr2);


for ii=1:ss3
    for jj=1:ss4
        for kk=1:ss2
            
            tally_em2(km_lof_gof(find(tr(ii,:))'),drivers(find(tr2(jj,:)')),hit_antihit(find(dbl(kk,:)')),dd,fid,fid_str);
        end
    end
end

fclose(fid);