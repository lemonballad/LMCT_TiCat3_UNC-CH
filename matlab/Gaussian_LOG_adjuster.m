%% Define
path_source='C:\Users\Thomas\Desktop'; % An example would be 'C:\Users\Thomas\Desktop\GaussView 3.09';
path_old=path_source; % Same as source directory
path_new='C:\Users\Thomas\Desktop'; % An example would be 'C:\Users\Thomas\Desktop\GaussView 3.09';

% The name of the file to be corrected
filename_source='TiCat3';
% The name of the file to conatain the old data
filename_old=[filename_source '_old'];
% The name of the file to contain the new data
filename_new=[filename_source '_new'];
% Suffix of file usually .log or .out
file_ext='.log';

%% Corrects if path names were not entered
if strcmp(path_source,'')    
    file_source=[filename_source file_ext];
else
    file_source=[path_source '\' filename_source file_ext];
end

if strcmp(path_old,'')    
    file_old=[filename_old file_ext];
else
    file_old=[path_old '\' filename_old file_ext];
end

if strcmp(path_new,'')    
    file_new=[filename_new file_ext];
else
    file_new=[path_new '\' filename_new file_ext];
end

%% Corrects if filenames weren't entered
if strcmp(filename_source,'')
    skip_flag=1; % Will skip correction
else
    skip_flag=0; % Continue with correction
end;
if strcmp(filename_old,''),filename_old=[filename_source '_old'];end;
if strcmp(filename_new,''),filename_new=[filename_source '_new'];end;
if strcmp(file_ext,'')&&~skip_flag
    skip_flag=1; % Will skip correction
else
    skip_flag=0; % Will continue with correction
end;

%% Copy contents of old file and correct for new file
if ~skip_flag
    file_contents_old=fileread(file_source);
    file_contents_new=regexprep(file_contents_old,'Atom  AN','Atom AN');
    
    fout_old=fopen(file_old,'w+');
    fout_new=fopen(file_new,'w+');
    
    fprintf(fout_old,'%s',file_contents_old);
    fprintf(fout_new,'%s',file_contents_new);
    
    fclose('all');
end
