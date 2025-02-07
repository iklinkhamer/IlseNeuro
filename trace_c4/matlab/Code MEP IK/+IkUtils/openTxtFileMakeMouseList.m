ff = fileread('Shank2 KO MUT');
f2 = regexprep(ff,'\s+',"', '");
fid = fopen('/home/mick/Desktop/Ilse/Stuff/Shank2 KO MUT','w');
fprintf(fid, '%s', f2);
fclose(fid);





% % fid  = fopen('Shank2 KO MUT','r');
% %f=fread(fid, 'char*1');
% ff = fileread('Shank2 KO MUT');
% %f2 = strrep(f, "-08\s+", "'MI2");
% f2 = regexprep(ff,'\s+',"', '");
% %f3 = regexprep(ff,'\s+',"'\s+'");
% fid = fopen('/home/mick/Desktop/Ilse/Stuff/Shank2 KO MUT','w');
% %fid2 = fopen('/home/mick/Desktop/Ilse/Stuff/MouseNames3','w');
% fprintf(fid, '%s', f2);
% %fprintf(fid2, '%s', f3);
% fclose(fid);


