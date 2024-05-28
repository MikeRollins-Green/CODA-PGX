function msigni=col_powered(msigni)

msigni.powered=logical(msigni.filt)&isfinite(msigni.p);

% dump powered when have gof in filt (ignore data in lines where don't have gof mutation)
for ii=1:length(msigni.rowlabels)
	for jj=1:length(msigni.collabels)
		
		if max(msigni.filt(ii,jj,:))==2
			msigni.powered(ii,jj,:)=msigni.filt(ii,jj,:)==2&isfinite(msigni.p(ii,jj,:));
		end
	end
end
