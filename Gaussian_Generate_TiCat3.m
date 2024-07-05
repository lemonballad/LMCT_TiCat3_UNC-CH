path_source='C:\Users\Thomas\Desktop\tic';
path_dest='C:\Users\Thomas\Desktop\tic';
filename_source='TiCat3_1';
filename_prefix_dest='TiCat3_';
file_ext_gjf='.gjf';
file_source_gjf=[path_source '\' filename_source file_ext_gjf];

file_contents=fileread(file_source_gjf);

for state=2:20
    file_contents_current=regexprep(file_contents,'%chk=TiCat3_1.chk',['%chk=TiCat3_' int2str(state) '.chk']);
    file_contents_current=regexprep(file_contents_current,'nstates=20,root=1)',...
        ['nstates=20,root=' int2str(state) ')']);

    file_dest_gjf=...
        [path_dest '\' filename_prefix_dest int2str(state) file_ext_gjf];
    
    fout=fopen(file_dest_gjf,'w+');

    fprintf(fout,'%s',file_contents_current);

    fclose(fout);
end

fclose('all');