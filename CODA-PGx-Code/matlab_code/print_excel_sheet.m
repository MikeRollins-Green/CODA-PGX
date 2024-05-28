function print_excel_sheet(all_data,filename)

% print_excel_sheet(all_data,filename)


pid=fopen(filename,'w');

fprintf(pid,['Rowlabels']);

for n=1:length(all_data.collabels)
    fprintf(pid,['\t' all_data.collabels{n}]);
end
fprintf(pid,'\n');

for n=1:length(all_data.rowlabels)
    
%     if n==1
%         fprintf(pid,['']);
%     else
        fprintf(pid,[all_data.rowlabels{n}]);
%     end
    
    for j=1:length(all_data.collabels)
        
        fprintf(pid,['\t' num2str(all_data.data(n,j)) ]);
        
    end
    fprintf(pid,'\n');
end

fclose(pid);
