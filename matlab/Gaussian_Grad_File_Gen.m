%% Define files and parameters

norm_modes=[27 39];%[6 10 27 39];

path_read='C:\Users\Thomas\Desktop\TiCat3\Gaussian\Input\6-21-16';
path_write='C:\Users\Thomas\Desktop\TiCat3\Gaussian\Input\6-21-16';
molecule='TiCat3_C3';
file_ext_gjf='.gjf';
file_ext_log='.log';
file_ext_chk='.chk';
input_file_0=[path_read '\' molecule '_0' file_ext_gjf];
output_file_0=[path_read '\' molecule '_0' file_ext_log];
chk_file_0=[path_read '\' molecule '_0' file_ext_chk];

%% Read normal mode displacements
start_line=' Full mass-weighted force constant matrix:';
end_line=' - Thermochemistry -';
line_prefix='\s+(?<Atom_num>\d+)\s+\d+\s+';
coord_collect_1='(?<x1>[-]*\d+[.]\d+)\s+(?<y1>[-]*\d+[.]\d+)\s+(?<z1>[-]*\d+[.]\d+)';
coord_collect_2='(?<x2>[-]*\d+[.]\d+)\s+(?<y2>[-]*\d+[.]\d+)\s+(?<z2>[-]*\d+[.]\d+)';
coord_collect_3='(?<x3>[-]*\d+[.]\d+)\s+(?<y3>[-]*\d+[.]\d+)\s+(?<z3>[-]*\d+[.]\d+)';
match_line=[line_prefix coord_collect_1 '\s+' coord_collect_2 '\s+' coord_collect_3];

start_flag=false;
end_flag=false;
last_index=-1;
NMode=0;

f_in=fopen(output_file_0);
line=fgetl(f_in);
while ischar(line)
    if regexp(line,start_line),start_flag=true;end
    if start_flag&&~end_flag
        NMCoord_DATA=regexp(line,match_line,'names');
        if size(NMCoord_DATA,1)~=0
            Atom_Index=str2double(NMCoord_DATA.Atom_num);
        else
            Atom_Index=-1;
        end
        if Atom_Index~=last_index&&Atom_Index~=-1
            NMCoords(Atom_Index,1:3,NMode+1)=[str2double(NMCoord_DATA.x1) str2double(NMCoord_DATA.y1) str2double(NMCoord_DATA.z1)];
            NMCoords(Atom_Index,1:3,NMode+2)=[str2double(NMCoord_DATA.x2) str2double(NMCoord_DATA.y2) str2double(NMCoord_DATA.z2)];
            NMCoords(Atom_Index,1:3,NMode+3)=[str2double(NMCoord_DATA.x3) str2double(NMCoord_DATA.y3) str2double(NMCoord_DATA.z3)];
            if Atom_Index<last_index,NMode=NMode+3;end
            last_index=Atom_Index;
            NMCoord_DATA=[];
        end
        if regexp(line,end_line),end_flag=true;end
    end
    line=fgetl(f_in);
end
fclose('all');

NAtoms=size(NMCoords,1);

%% Read ground state cartesian coordinates

f_in=fopen(input_file_0);
line=fgetl(f_in);
Atom(1:NAtoms)=cellstr(' ');
x(1:NAtoms)=0;
y(1:NAtoms)=0;
z(1:NAtoms)=0;
atom_index=0;
cart_line='\s+(?<Atom>\S+)\s+(?<x>[-]*\d+[.]\d+)\s+(?<y>[-]*\d+[.]\d+)\s+(?<z>[-]*\d+[.]\d+)';
header=[];
tail=[];
prefix_flag=true;
while ischar(line)
    if length(line)==58
        atom_index=atom_index+1;
        Coord_DATA=regexp(line,cart_line,'names');
        Atom(atom_index)=cellstr(Coord_DATA.Atom);
        x(atom_index)=str2double(Coord_DATA.x);
        y(atom_index)=str2double(Coord_DATA.y);
        z(atom_index)=str2double(Coord_DATA.z);
        prefix_flag=false;
    elseif prefix_flag
        header=[header line char(10)];
    else
        tail=[tail line char(10)];
    end
    line=fgetl(f_in);
end
fclose('all');

%% Write to gjf file

% Loop over electronic excitations
for root=1%1:6
    % Loop over modes of interest
    for iMode=1:length(norm_modes)
        % Account for which mode will be used to calculate 2nd derivatives
%         if iMode==1
%             disp=-0.02:0.01:0.02;
%         else
%             disp=-0.01:0.01:0.02;
            disp=-0.04:0.02:0.04;

%         end
        % Loop over steps away from equilibrium for derivatives
        for idisp=1:length(disp)
            % Write header of gjf file
            header_current=...
                regexprep(header,'%chk=TiCat3_C3_0.chk',['%chk=TiCat3_C3_'...
                int2str(norm_modes(iMode)) '_' num2str(disp(idisp)) '.chk']);
%             header_current=regexprep(header_current,'nstates=6,root=1)',...
%                 ['nstates=6,root=' int2str(root) ')']);
            
            % Write current body of gjf file
            body_current=[];
            for iatom=1:length(Atom)
                line_current=sprintf(' %-16s%14.8f%14.8f%14.8f\n'...
                    ,cell2mat(Atom(iatom))...
                    ,x(iatom)+disp(idisp)*NMCoords(iatom,1,iMode)...
                    ,y(iatom)+disp(idisp)*NMCoords(iatom,2,iMode)...
                    ,z(iatom)+disp(idisp)*NMCoords(iatom,3,iMode));
                body_current=[body_current line_current];
            end
            
            % Write contents of gjf file
            file_current=[header_current body_current tail];
            
            % Define path and name for current gjf file
            file_dest_gjf=...
                [path_write '\' molecule '_' int2str(norm_modes(iMode)) '_' num2str(disp(idisp)) file_ext_gjf];
            chk_file_new=[path_write '\' molecule '_' int2str(norm_modes(iMode)) '_' num2str(disp(idisp)) file_ext_chk];

                        % Open current gjf file with write permission
            fout=fopen(file_dest_gjf,'w+');
            % Write to gjf file
            fprintf(fout,'%s',file_current);
            % Close current gjf file
            fclose(fout);
            copyfile(chk_file_0,chk_file_new);
        end % End loop over displacement steps
    end % End loop over modes
end % End loop over electronic excitations

% Close all files
fclose('all');
