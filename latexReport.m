function [  ] = latexReport( out_name, errL2, errH1, h, nt, n, t)

fid = fopen(out_name,'wt');

doc_class = '\documentclass{standalone}';
fprintf(fid, '%s\n\n',doc_class);
begin_doc = '\begin{document}';
clearpage = '\clearpage';
fprintf(fid, '%s\n\n', begin_doc);
fprintf(fid, '%s\n', clearpage);

disp('[*] Writing Latex Table Report:');
fprintf('');

sh = 'h	&	$';

for i=1:n
   sh = [sh, num2str(h(i))];
   
   if (i==n)
      sh = [sh, '$ \\ \hline'];
   else
      sh = [sh, '$ & $'];
   end
end % end for i

for j=1:nt
fprintf(' - t = %s', num2str(t(j)));
    
sL2 = 'ErrL2	&	$';
sH1 = 'ErrH1	&	$';
    
   
for i=1:n
   sL2 = [sL2, num2str(errL2(i,j))];
   sH1 = [sH1, num2str(errH1(i,j))];
        
   if (i==n)
      sL2 = [sL2, '$ \\ \hline'];
      sH1 = [sH1, '$ \\ \hline'];
   else
      sL2 = [sL2, '$ & $'];
      sH1 = [sH1, '$ & $'];
    end
end % end for i
    
% Used tmp variables to avoid problems with the escape character \,
% which is pervasive in Latex commands
hfill = '\hfill \\';
fprintf(fid, '%s \n',hfill);
begin_table = '\begin{table}[!h]';
centering = '\centering';
fprintf(fid, '%s\n%s\n', begin_table, centering);
begin_tabular = '\begin{tabular}{ | c | c | c | c | c | c | c | }';
fprintf(fid, '%s\n', begin_tabular);
h_line = '\hline';
fprintf(fid, '%s\n', h_line);
fprintf(fid, '%s\n', sh);
fprintf(fid, '%s\n', sL2);
fprintf(fid, '%s\n', sH1);
end_tabular = '\end{tabular}';
caption = '\caption{t = ';
end_table = '\end{table}';
 fprintf(fid, '%s\n%s%s }\n%s\n\n', end_tabular, caption, num2str(t(j)), end_table);
    
fprintf(' ... OK \n');

end

fprintf(fid, '\n%s\n', clearpage);
end_doc = '\end{document}';
fprintf(fid, '\n\n%s', end_doc);

fclose(fid);


end

