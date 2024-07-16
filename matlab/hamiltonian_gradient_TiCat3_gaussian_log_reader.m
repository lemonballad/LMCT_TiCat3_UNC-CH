% Define paths to read/write file from/to, file prefix, file suffix, etc
path_read='C:\Users\Thomas\Desktop\TiCat3\Gaussian\Output\3-17-16';
path_write='C:\Users\Thomas\Desktop\TiCat3\Gaussian\Output\3-17-16';
molecule='TiCat3_1_';
file_ext_log='.log';

% Define array of normal modes
NMode=[6 10 27 39];

% Initialize data matrices
Delta_Ex=zeros(length(NMode),5,9,'double');
E0=zeros(length(NMode),5,'double');
Ex=zeros(length(NMode),5,9,'double');
mu=zeros(length(NMode),5,9,'double');
Osc=zeros(length(NMode),5,9,'double');

% Define expressions to flag or match
E0_expr=' SCF Done:\s+\S+\s+[=]+\s+(?<E0>[-]*\d+[.]+\d+)\s+';
ExSt_End_expr=' Leave Link  914 ';
ExSt_flag=false;
ExSt_flag_expr=' Excitation energies and oscillator strengths:';
Ex_pre_str='\s+Excited State\s+\d+[:]+\s+\S+\s+\d+[.]+\d+\s+eV\s+';
wln_str='(?<wln>\d+[.]+\d+)';
osc_str='(?<Osc>\d+[.]+\d+)';
ExSt_expr=[Ex_pre_str wln_str '\s+nm\s+f[=]+' osc_str '\s+'];
mu_End_expr='Ground to excited state transition velocity dipole moments';
mu_flag=false;
mu_flag_expr='Ground to excited state transition electric dipole moments';
mu_expr='\s+(?<state>\d+)\s+(?<mu_x>[-]*\d+[.]+\d+)\s+(?<mu_y>[-]*\d+[.]+\d+)\s+(?<mu_z>[-]*\d+[.]+\d+)\s+\d+[.]+\d+\s+\d+[.]+\d+';

% Loop over normal modes
for inmode=1:length(NMode)
    % Define array of displacements in normal modes
    if NMode(inmode)==6
        Disp=-0.02:0.01:0.02;
    else
        Disp=-0.01:0.01:0.01;
    end
    
    % Loop over displacement in normal modes
    for idisp=1:length(Disp)
        % Initialize excitation level
        lvl=0;
        state=0;
        % Define file name to read
        log_file=[path_read '\' molecule int2str(NMode(inmode)) '_' num2str(Disp(idisp)) file_ext_log];
        
        % Open file to read
        f_in=fopen(log_file);
        % Read in first line of file
        line=fgetl(f_in);
        
        % Loop over lines in file. If file is empty or end of file is hit
        % loop will not execute, (ischar==false).
        while ischar(line)
            % Read ground state energy in a.u
            Read_DATA=regexp(line,E0_expr,'names');
            if size(Read_DATA,1)~=0,E0(inmode,idisp)=2.2*10^5*str2double(Read_DATA.E0);end;
            
            % Read transition dipole
            if regexp(line,mu_flag_expr),mu_flag=true;end;
            if mu_flag
                Read_DATA=regexp(line,mu_expr,'names');
                if size(Read_DATA,1)~=0
                    state=str2double(Read_DATA.state);
                    mu(inmode,idisp,state)=...
                        sqrt(str2double(Read_DATA.mu_x).^2+...
                        str2double(Read_DATA.mu_y).^2+...
                        str2double(Read_DATA.mu_z).^2);
                end
                if regexp(line,mu_End_expr),mu_flag=false;end;
                if state==9,mu_flag=false;end;
            end
            
            % Read excited state energies and oscillator strengths
            if regexp(line,ExSt_flag_expr),ExSt_flag=true;end;
            if ExSt_flag
                Read_DATA=regexp(line,ExSt_expr,'names');
                if size(Read_DATA,1)~=0
                    lvl=lvl+1;
                    Ex(inmode,idisp,lvl)=10^7/str2double(Read_DATA.wln);
                    Delta_Ex(inmode,idisp,lvl)=Ex(inmode,idisp,lvl)+E0(inmode,idisp);
                    Osc(inmode,idisp,lvl)=str2double(Read_DATA.Osc);
                end
                if regexp(line,ExSt_End_expr),ExSt_flag=false;end;
                if lvl==9,ExSt_flag=false;end;
            end
            
            % Initialize and read in next line from file
            Read_DATA=[];
            line=fgetl(f_in);
        end % End loop over lines in file
        % Close files
        fclose('all');
    end % End loop over displacements in normal mode
end % End loop over normal modes

% Adjust for variable length data
Delta_Ex(2:end,2:4,:)=Delta_Ex(2:end,1:3,:);
Delta_Ex(2:end,[1 5],:)=NaN;
E0(2:end,2:4)=E0(2:end,1:3);
E0(2:end,[1 5],:)=NaN;
Ex(2:end,2:4,:)=Ex(2:end,1:3,:);
Ex(2:end,[1 5],:)=NaN;
mu(2:end,2:4,:)=mu(2:end,1:3,:);
mu(2:end,[1 5],:)=NaN;
Osc(2:end,2:4,:)=Osc(2:end,1:3,:);
Osc(2:end,[1 5],:)=NaN;

% Compute extinction coefficient
omega_max=25000;
lnwidth=3000;
fwhm=2.355*lnwidth;
epsilon=mu.^2*omega_max/0.009186/sqrt(2*pi)/lnwidth;

% Compute electronic coupling
dist=3.5848;
H_AD=0.0206/dist*sqrt(epsilon*omega_max*fwhm);
H_ad_rev=(Ex(:,:,6)-Ex(:,:,4))/3
diag_el=(Ex(:,:,6)+2*Ex(:,:,4))/3;

% Compute first derivative electron coupling
% Compute gradient 
for inmode=1:4
    for idisp=2:4
        for ilvl=1:9
            dH_AD(inmode,ilvl)=(abs(diff(H_AD(inmode,2:3,ilvl)))+abs(diff(H_AD(inmode,3:4,ilvl))))/0.02;
            grad_E(inmode,ilvl)=(abs(diff(Ex(inmode,2:3,ilvl)))+abs(diff(Ex(inmode,3:4,ilvl))))/0.02;
        end
        dH_ad_rev(inmode)=(abs(diff(H_ad_rev(inmode,2:3)))+abs(diff(H_ad_rev(inmode,3:4))))/0.02;
        dH_ad_rev2(inmode)=((diff(H_ad_rev(inmode,2:3)))+(diff(H_ad_rev(inmode,3:4))))/0.02;
    end
end

for inmode=1:4
    for idisp=2:4
        for ilvl=1:9
            dhad(inmode,ilvl)=(H_AD(inmode,2,ilvl)-H_AD(inmode,4,ilvl))/-0.02;
            grade(inmode,ilvl)=(Ex(inmode,2,ilvl)-Ex(inmode,4,ilvl))/-0.02;
        end
        dhadrev(inmode)=(H_ad_rev(inmode,2)-H_ad_rev(inmode,4))/-0.02;
    end
end


