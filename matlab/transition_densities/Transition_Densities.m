tic
read_flag=false;
write_flag=false;
path_read='C:\Users\Thomas\Desktop\TiCat3\Gaussian\Good Result\El_C3';
path_write='C:\Users\Thomas\Desktop\TiCat3\Gaussian\Good Result\El_C3';
molecule='TiCat3_TD';
file_ext_cub='.cub';
file_ext_log='.log';
cube_1=[path_read '\' molecule '_1' file_ext_cub];
log_1=[path_read '\' molecule '_1' file_ext_log];

[NAtoms,r0,N,e,AtNum,q,r,rho_1]=read_cube(cube_1);
toc

f_in=fopen(cube_1);
line=fgets(f_in);
line_num=1;
header='';
while ischar(line)&&line_num<=6+NAtoms
    header=[header line];
    line=fgets(f_in);
    line_num=line_num+1;
end
fclose('all');

start_line=' Excitation energies and oscillator strengths:';
exp_match=' Excited State\s+(?<lvl>\d+):\s+\S+\s+(?<element_1>\d+[.]\d+) eV\s+(?<element_2>\d+[.]\d+) nm\s+f=(?<element_3>\d+[.]\d+)\s+[<]S[*]+2[>]=\d+[.]\d+';
start_flag=false;
end_flag=false;
f_in=fopen(log_1);
line=fgets(f_in);
while ischar(line)
    if regexp(line,start_line),start_flag=true;end
    if start_flag&&~end_flag
        len_line=length(line);
        if len_line==85
            UV_DATA=regexp(line,exp_match,'names');
            lvl=str2double(UV_DATA.lvl);
            eV(lvl)=str2double(UV_DATA.element_1); %#ok<*SAGROW>
            wvln(lvl)=str2double(UV_DATA.element_2);
            Osc_f(lvl)=str2double(UV_DATA.element_3);
            if lvl==20,end_flag=true;end
            UV_DATA=[];
            lvl=[];
        end
    end
    line=fgets(f_in);
end
fclose('all');

center=[r(1,:);sum(r([9 12],:),1)/2;sum(r([20 23],:),1)/2;sum(r([30 33],:),1)/2]';
radius=[4.0 3.5 3.5 3.5];

X=e(1,1)*(0:N(1)-1)+r0(1);
Y=e(2,2)*(0:N(2)-1)+r0(2);
Z=e(3,3)*(0:N(3)-1)+r0(3);

if read_flag,trans_rho=zeros(N(3),N(2),N(1),20,'double');end %#ok<*UNRCH>
cond_trans_rho=zeros(4,20,'double');
err_term(1:20)=0;
% Loop over cube files/electronic energy levels
parfor lvl=1:20
    tic
    if read_flag
        current_cube=[path_read '\' molecule '_' int2str(lvl) file_ext_cub];
        [~,~,~,~,~,~,~,rho]=read_cube(current_cube);
        trans_rho(:,:,:,lvl)=rho-rho_1;
        rho=[];
    end
    [cond_trans_rho(:,lvl),err_term(lvl)]=condense_transition_density(center,radius,X,Y,Z,trans_rho(:,:,:,lvl));
    if write_flag
        file_write=[path_write '\' molecule '_' int2str(lvl) '_1' file_ext_cub];
        write_cube(header,file_write,trans_rho(:,:,:,lvl));
    end
    toc
end
 fclose('all');