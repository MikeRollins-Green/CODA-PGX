function signi=col_signi(msigni)

ss=length(msigni.rowlabels);
tt=length(msigni.collabels);

signi=msigni;
signi.p=nan(ss,tt);
signi.fc=nan(ss,tt);
signi.n=nan(ss,tt);

% collapse to 2-dim p
for ii=1:ss
    for jj=1:tt
        
        my_pow=msigni.powered(ii,jj,:);
        if max(my_pow)==1
            signi.p(ii,jj)=prod(msigni.p(ii,jj,my_pow));
            %signi.p(ii,jj)=min(msigni.p(ii,jj,my_pow));
            signi.n(ii,jj)=sum(my_pow);
        end
    end
end

% collapse to 2-dim fc
for ii=1:ss
    for jj=1:tt
        
        my_pow=msigni.powered(ii,jj,:);
        if max(my_pow)==1
            
            fc_temp=msigni.fc(ii,jj,my_pow);
            hw=abs(fc_temp)>=1;
            
            % make sure all the signs align, then take max
            if max(fc_temp(hw))<0 % all high weight vals negative
                signi.fc(ii,jj)=min(fc_temp);
            elseif min(fc_temp(hw))>0 % all high weight vals positive
                signi.fc(ii,jj)=max(fc_temp);
            else % all values are >-1 and <1
                
                if sum(fc_temp<0)/length(fc_temp)>0.5 % mostly negative
                    signi.fc(ii,jj)=min(fc_temp);
                else
                    signi.fc(ii,jj)=max(fc_temp); % mostly positive
                end
                
            end
        end
    end
end

signi.filt=max(msigni.filt,[],3);
