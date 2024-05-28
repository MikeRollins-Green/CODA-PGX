function pos=min_rr(pos,min_rr)

cl1=unique(pos.cell_line);
dr=unique(pos.drug);
pl=unique(pos.pt_pool);

keep=zeros(length(pos.collabels),1);

for ii=1:length(cl1)
    for jj=1:length(dr)
        for kk=1:length(pl)
            
            % same line, drug, and pool
            oo1=find(strcmp(cl1(ii),pos.cell_line)&strcmp(dr(jj),pos.drug)&pos.pt_pool==pl(kk));
            
            if ~isempty(oo1)
                rr=corr(pos.data(:,oo1));
                
                
                [ss,~]=size(rr);
                
                % dump the self-hits
                for iii=1:ss
                    for jjj=iii
                        rr(iii,jjj)=0;
                    end
                end
                
                keep_rr=max(rr,[],2)>=min_rr;
                keep(oo1(keep_rr))=1;

            end
        end
    end
end

pos=downSize(pos,keep==1);





