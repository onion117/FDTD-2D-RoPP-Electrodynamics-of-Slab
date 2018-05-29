clear all
clc
addpath EM_functions        % folder with sub-functions
addpath material_data       % folder with material data

% Define Models the simulation will run through each model, all must be 3
% charecters long i.e 'AB '

    models=cellstr([ 'AB ';'MN ';'AMP';'EL ';'Chu']); % Run all models
    %models=cellstr([ 'Chu';]);

   
% Define which kind of slab the pulse is going through
% All strings must be the same length
    simcases=cellstr([  'Lossless        ';...
                        'matched_lossy   ';...
                        'Ag_reflector    ';...
                        'negative_matched';]);                 
 
    sim_case=char(simcases(4));                       % This will simulate Lossless
                                                % a negative index
                                                % Impedance matched slab
                                                                

    save_mode=1;                                % Turn on '1' save all workspace data
                                                % 0 dont save, this is
                                                % workspace not .dat text
                                                % files
for sim_j=1:length(models)
    
    close all  
    keep models sim_j save_mode sim_case       % clear everything except which model to use.
    model=models(sim_j);
% FIGURE INITIALIZATION

    fig_count=200;              % time steps per plot
    plot_on=1;                  % Turn plotting on or, off, simulation much faster if off
    n_count=0;                  % Plot Counter
    s_d=get(0,'ScreenSize');    % Screen [0 0 width height]
    sw=s_d(3);                  % Screen width
    sh=s_d(4);                  % Screen height

% Figure Positions, relatibe to screen highe 'sh' and screen width 'sw'
% convention is [xo,yo,width,height] in pixels

    p1 = [5 45 sw/2 sh/1.8];
    p2=[5+p1(3) 45 sw/2 sh/1.1];
    p3=[5 sh/2+15 sw/3 (sh-35)/2];
    p4=[p3(1)+p3(3) sh/2+15 sw/3 (sh-35)/2];
    p5=[p4(1)+p4(3) sh/2+15 sw/3 (sh-35)/2];
    p6=[5 45 sw/3 sh/2.1];
    p7=[p6(3) 45 sw/3 sh/2.1];
    p8=[p6(3)+p7(3) 45 sw/3 sh/2.1];

    fig_1=figure(1);                
    set(fig_1,'name','Source Properties')   
    set(fig_1,'DefaulttextFontSize',14)
    set(fig_1,'OuterPosition',p1)

    fig_2=figure(2);
    set(fig_2,'name','Material Properties')
    set(fig_2,'DefaulttextFontSize',14)
    set(fig_2,'OuterPosition',p2)

    fig_3=figure(3);
    set(fig_3,'name','Hz field')
    set(fig_3,'DefaulttextFontSize',14)
    set(fig_3,'OuterPosition',p3)

    fig_4=figure(4);
    set(fig_4,'name','X-Momentum')
    set(fig_4,'DefaulttextFontSize',14)
    set(fig_4,'OuterPosition',p4)

    fig_5=figure(5);
    set(fig_5,'name','Y-Momentum')
    set(fig_5,'DefaulttextFontSize',14)
    set(fig_5,'OuterPosition',p5)

    fig_6=figure(6);
    set(fig_6,'name','X Center of Mass')
    set(fig_6,'DefaulttextFontSize',14)
    set(fig_6,'OuterPosition',p6)

    fig_7=figure(7);
    set(fig_7,'name','Y Center of Mass')
    set(fig_7,'DefaulttextFontSize',14)
    set(fig_7,'OuterPosition',p7)


    fig_8=figure(8);
    set(fig_8,'name','Total_Energy')
    set(fig_8,'DefaulttextFontSize',14)
    set(fig_8,'OuterPosition',p8)

    fig_9=figure(9);
    set(fig_9,'name','Displacement')
    set(fig_9,'DefaulttextFontSize',14)
    set(fig_9,'OuterPosition',p8)


% DEFINE SI UNITS

    meters=1;
    nm=meters*1e-9;
    fs=1e-15;
    mu_o=4*pi*10^-7;
    c=299792458;
    eps_o=(1/(c*c*mu_o));
    eta_o=sqrt(mu_o/eps_o);
    Kg=1;


%% Simulation Parameters

    LAMBDA=500*nm;                  % Center wavelength of pulse
    f=c./LAMBDA;                    % Center frequency of pulse
    T=1./f;                         % Center period

    dx=28*nm;
    dy=28*nm;
    
% TEMPORAL GRID
    
    dt=dx/(2*c);
    sim_cycles=38;
    sim_time=(sim_cycles*(1/f));
    Nt=round(sim_time/dt);

% GEOMETRY

    M_slab=.00005*Kg;                   % Mass of slab
    M_system=M_slab;                    % Assume mass of system is same slab
    
    spc_x1=6500*nm;                     % Space on the left of slab
    spc_x2=500*nm;                      % Space on the right of slab  
    spc_y1=200*nm;                      % Bottom space  (y=0 and grater)  below slab            
    spc_y2=200*nm;                      % Space between slab and Y_PML  above slab
    
    source_type='x';                    % source propagation direction
    gauss_y_width=3*LAMBDA;             % Width of gaussian beam
    gauss_x_width=1*LAMBDA;             % Width of gaussian beam
    
    NPML=1.0*LAMBDA;                    % PML size
    NPML_x=round(NPML/dx);              % size of PML along x    
    NPML_y=round(NPML/dy);              % size of PML along y    

    object_type='box';                  % set tp 'box' or 'circle'
    tilt_in_degrees=45;
    
    L_x=3.5*gauss_x_width;              % x Length of object
    L_y=4.5*gauss_y_width;              % y length of object    
    
    x_size=spc_x1+spc_x2+max([ L_x gauss_x_width]);% x[m] size of simulation
    y_size=max([L_y gauss_y_width])+spc_y1+spc_y2; % y[m] size of simulation
    
    Nx=round(x_size/dx)+2*NPML_x;           % Grid x length
    Ny=round(y_size/dy)+2*NPML_y;           % Grid y length

    x=[0:1:Nx-1]*dx;                        % x axis 
    y=[0:1:Ny-1]*dy;                        % y axis 
    t=[0:dt:sim_time];                      % time axis in [s]    
    
    c_x=NPML_x*dx+spc_x1+L_x/2-L_y/4;        % center x location of slab
    c_y=mean(y);                             % center y of grid used before rotation
    c_y_object=c_y;                             % y location of slab
   
    gauss_x_avg=c_x;                    
    gauss_y_avg=c_y_object-600*nm;

    source1_x=NPML_x+round(5*dx/dx);            % x-location of source
    source1_y=Ny-NPML_y-round(spc_y2/(5*dy));   % y-location of source


%% Sim Cases
switch sim_case
    
case ('Ag_reflector')
    
    er_1=1;                             % Background \eps_r
    mr_1=1;                             % Background \mu_r      

    er_slab=4;                              
    mr_slab=1;                                

    gamma_e_material=6E13;                
    omega_e_material=1.4E16;    

    gamma_m_material=0;                 % damping frequency of object
    omega_m_material=0;                 % electric plasma frequency of

    sigma_m_material=0;
    sigma_e_material=0;    
    
case ('Lossless') 

    er_1=1;                             % Background \eps_r
    mr_1=1;                             % Background \mu_r        

    er_slab=3;                              
    mr_slab=1;                                

    gamma_e_material=0E13;                
    omega_e_material=.25E16;    

    gamma_m_material=0*5E14;                % damping frequency of object
    omega_m_material=0*6E15;                % electric plasma frequency of object

    sigma_e_material=0;                     % Electrical Conductivity
    sigma_m_material=0;                     % Magnetic Conductivity    

    
case ('matched_lossy') 
    
    er_1=1;                             % Background \eps_r
    mr_1=1;                             % Background \mu_r      

    er_slab=1.4;                             % e_inf for the slab
    mr_slab=1.4;                             % m_inf for the slab 

    gamma_e_material=4.5*5E14;                
    omega_e_material=.2*6E15;    

    gamma_m_material=4.5*5E14;               % damping frequency of object
    omega_m_material=.2*6E15;                % electric plasma frequency of

    sigma_m_material=0;
    sigma_e_material=0;


case ('negative_matched') 
    
    er_1=1;                             % Background \eps_r
    mr_1=1;                             % Background \mu_r      

    er_slab=1.02;                             % e_inf for the slab
    mr_slab=1.02;                             % m_inf for the slab 

    gamma_e_material=0*5E14;                
    omega_e_material=6.0000000E+15;    

    gamma_m_material=0*5E14;               % damping frequency of object
    omega_m_material=6.0000000E+15;                % electric plasma frequency of

    sigma_m_material=0;
    sigma_e_material=0;

case ('electric_lossy') 
    
    er_1=1;                             % Background \eps_r
    mr_1=1;                             % Background \mu_r      

    er_slab=3;                             % e_inf for the slab
    mr_slab=1;                             % m_inf for the slab 

    gamma_e_material=1.5*5E14;                
    omega_e_material=.3*6E15;    

    gamma_m_material=0;               % damping frequency of object
    omega_m_material=0;                % electric plasma frequency of

    sigma_m_material=0;
    sigma_e_material=0;

 end




%% Material Create

    er=er_1.*ones(Nx,Ny);
    mr=mr_1.*ones(Nx,Ny);

  switch object_type
      
      case('box')
    
    [er,cx1,cx2,cy1,cy2]=create_box(L_x,L_y,c_x,c_y,er_slab,er,dx,dy);
    [mr                ]=create_box(L_x,L_y,c_x,c_y,mr_slab,mr,dx,dy);
    
    
    case('circle')  
        [er,cx1,cx2,cy1,cy2]=create_circle(max([L_x/2,L_y/2]),c_x,c_y,er_slab,er,dx,dy);
        [mr                ]=create_circle(max([L_x/2,L_y/2]),c_x,c_y,mr_slab,mr,dx,dy);
%     
  end
  material_location_matrix=1.*(er==er_slab);
  if ((tilt_in_degrees~=0)&&(strcmp('box',object_type)))
    
    material_location_matrix=imrotate(material_location_matrix,45,'bilinear','crop');
    
    material_location_matrix=1.*(material_location_matrix~=0);
    material_location_matrix=center_y_at(y,material_location_matrix,c_y_object);
    mat_er1=er_slab.*material_location_matrix;
    mat_er2=er_1.*(mat_er1~=er_slab);
    mat_mr1=mr_slab.*material_location_matrix;
    mat_mr2=mr_1.*(mat_mr1~=mr_slab);
    er=mat_er1+mat_er2;
    mr=mat_mr1+mat_mr2;
    
    x_proj=(L_y/2)*cos((2*pi*tilt_in_degrees)/(360));
    y_proj=(L_x/2)*cos((2*pi*tilt_in_degrees)/(360));
%     
     cx1=cx1-round(x_proj/dx);
     cx2=cx2+round(x_proj/dx);
     
     cy1=cy1-round(y_proj/dy);
     cy2=cy2+round(y_proj/dy);
   
  end
    % Initialize drude matricies
    gamma_e=gamma_e_material.*material_location_matrix;
    omega_e=omega_e_material.*material_location_matrix;
    gamma_m=gamma_m_material.*material_location_matrix;
    omega_m=omega_m_material.*material_location_matrix;

    sigma_m=sigma_m_material.*material_location_matrix;
    sigma_e=sigma_e_material.*material_location_matrix;
    
    
% add two cells around the slab to include in slabe calculations
% This will catch derivative terms that extend +2 dx at the boundary

    [Fx,Fy]=gradient(material_location_matrix);
    mat_1=1.*(Fx~=0);
    mat_2=1.*(Fy~=0);
    slab_mat=material_location_matrix+mat_1+mat_2;
    slab_mat=1.*(~(slab_mat==0));

    [Fx,Fy]=gradient(slab_mat);
    mat_1=1.*(Fx~=0);
    mat_2=1.*(Fy~=0);
    slab_mat=slab_mat+mat_1+mat_2;
    slab_mat=1.*(~(slab_mat==0));
    system_mat=1.*(~(slab_mat==1));

    clear mat_1 mat_2 mat_er1 mat_er2 mat_mr1 mat_mr2


%% PML Implement

    sigma_PML_min=0;                          % Taper PML conductivity between min and max
    sigma_PML_max=90000;                      % Found by trial and error           

    eps_1=er(NPML_x+2,NPML_y+2);        	% Sample Surrounding Medium, call this medium 1     
    mu_1=mr(NPML_x+2,NPML_y+2);
     
    eps_2=(eps_1);                          % Solve for Medium 2 in PML, assuming eps2=eps1
    mu_2=(eps_2/eps_1)*mu_1;                % This will only work if mr=1 for surrounding materials

    m_fac=(mu_2/(eps_2))*(mu_o/eps_o);      % sigma_m=m_fac*sigma_e
    
    % X PML
    er(1:NPML_x,:)=eps_2;
    er((Nx-NPML_x):Nx,:)=eps_2;
    mr(1:NPML_x,:)=mu_2;
    mr((Nx-NPML_x):Nx,:)=mu_2;
    
    % Y PML
    er(:,1:NPML_y)=eps_2;
    er(:,(Ny-NPML_y):Ny)=eps_2;
    mr(:,1:NPML_y)=mu_2;
    mr(:,(Ny-NPML_y):Ny)=mu_2;  

    slope_x=(sigma_PML_min-sigma_PML_max)/(NPML_x-1);
    slope_y=(sigma_PML_min-sigma_PML_max)/(NPML_y-1);
    b_x=sigma_PML_min-slope_x*NPML_x;
    b_y=sigma_PML_min-slope_y*NPML_y;
    
    
    % set x PML
for i=1:NPML_x
    
        sigma_e(i,:)=slope_x*i+b_x;        
        sigma_e(Nx-i+1,:)=slope_x*i+b_x;
        sigma_m(i,:)=m_fac*sigma_e(i,:);
        sigma_m(Nx-i+1,:)=m_fac*sigma_e(Nx-i+1,:);
    
end

    % set y PML
for j=1:NPML_y
    

        sigma_e(:,j)=slope_y*j+b_y+sigma_e(:,j);        
        sigma_e(:,Ny-j+1)=slope_y*j+b_y+sigma_e(:,Ny-j+1);
        sigma_m(:,j)=m_fac*sigma_e(:,j);
        sigma_m(:,Ny-j+1)=m_fac*sigma_e(:,Ny-j+1);
    
end
    
 %% Source Properties 

    a=5000000000000000;                          % amplitude scaleing
    f_start=.2E15;                              % lowest frequency component   
    f_end=1.2E15;                               % highest frequency component
    FWHM=.2*1/2.9e-15;                             % Full Width Half Maximum
    N=50;                                       % Number of carrier signals
    df=(f_end-f_start)/N;                       % Frequency spaceing
    tau=5*T;                                    % shift pulse in time domain
    G=gauss_pulse(tau,f,FWHM,N,f_start,f_end);  % Create Spectrum
    norm_G=max((G(:,2)));                       % Normalize sourcespc/nm

    G(:,2)=G(:,2)./norm_G;                      % Normalize..                            

% Multiple by pulse in time for control over when pulse ends.
    t_avg=(tau);            
    sig_t=1*(1/f);   
    gauss_t=1*exp(-1.*((t-t_avg).^2)/(2*sig_t^2));  

% Y-Source GAUSSIAN PROFILE

    sig_y=gauss_y_width/1;
    sig_y=sig_y/2;


    N_y_gauss=round(gauss_y_width/dy);
    y_gauss=[1:Ny]*dy;

    gauss_y=gauss_create(gauss_y_avg,sig_y,y_gauss);
    gauss_y=gauss_y./max(gauss_y);

    gauss_y_start=gauss_y_avg-4.0*sig_y;
    gauss_y_end=gauss_y_avg+4.0*sig_y;

    nyg1=round(gauss_y_start/dy);
    nyg2=round(gauss_y_end/dy);

% Y-Source GAUSSIAN PROFILE

    sig_x=gauss_x_width/1;
    sig_x=sig_x/2;


    N_x_gauss=round(gauss_x_width/dx);
    x_gauss=[1:Nx]*dx;

    gauss_x=gauss_create(gauss_x_avg,sig_x,x_gauss);
    gauss_x=gauss_x./max(gauss_x);

    gauss_x_start=gauss_x_avg-1.0*sig_x;
    gauss_x_end=gauss_x_avg+1.0*sig_x;

    nxg1=round(gauss_x_start/dx);
    nxg2=round(gauss_x_end/dx);




% Store drude model in space and check its locatons

    er_d=zeros(Nx,Ny);  % Drude permitivity
    mr_d=zeros(Nx,Ny);  % Drude permeability

    parfor i=1:Nx
       for j=1:Ny


        lam0=LAMBDA;
        er_d(i,j)=drude_calc(gamma_e(i,j),omega_e(i,j),er(i,j),lam0);
        mr_d(i,j)=drude_calc(gamma_m(i,j),omega_m(i,j),mr(i,j),lam0);

        end
    end

    eta_d=sqrt(mr_d./er_d);


% check frequency drude
    er_w=zeros(1,length(G(:,1)));
    mr_w=zeros(1,length(G(:,1)));
    LAMBDA_drude=zeros(size(G(:,1)));
       
for j=1:length(G(:,1))
    lam0=c./G(j,1);
    i_test=cx1+round((cx2-cx1)/2);
    j_test=cy1+round((cy2-cy1)/2);
    er_w(1,j)=drude_calc(gamma_e(i_test,j_test),omega_e(i_test,j_test),er(i_test,j_test),lam0);
    mr_w(1,j)=drude_calc(gamma_m(i_test,j_test),omega_m(i_test,j_test),mr(i_test,j_test),lam0);
    LAMBDA_drude(j)=lam0;
end

% Drude dispersion
eta_w=sqrt(mr_w./er_w);


% Figure 1: X-Source Properties

figure(1)
subplot(1,2,1)
    stem(c./(G(:,1)*nm),abs(G(:,2)),'color','black')
    xlabel('\lambda')
    ylabel(' Signal Power ')
subplot(1,2,2)
    plot(LAMBDA_drude*1E9,real(er_w),'color','r')
    hold on
    plot(LAMBDA_drude*1E9,imag(er_w),'--','color','r')
    hold on
    plot(LAMBDA_drude*1E9,real(mr_w),'--','color','b')
    hold on
    plot(LAMBDA_drude*1E9,imag(mr_w),'color','b')
    xlabel(' Wavelength (nm) ')
    xlim([ 200 1100])
    ylabel(' \epsilon_r(\lambda) ')
    legend('Re(\epsilon_r)','Im(\epsilon_r)','Re(\mu_r)','Im(\mu_r)')
    
%% Ouput Fields
% Initialize Fields and Current matricies
    Hz=zeros(Nx,Ny);
    Bz=zeros(Nx,Ny);
    Hz_at_x=zeros(Nx,Ny);
    Hz_at_y=zeros(Nx,Ny);    
    Bz_at_x=zeros(Nx,Ny);
    Bz_at_y=zeros(Nx,Ny);   
    
    Ex=zeros(Nx,Ny);
    Ey=zeros(Nx,Ny);
    Dx=zeros(Nx,Ny);
    Dy=zeros(Nx,Ny);

    % Dispersive current densities
    Jxd=zeros(Nx,Ny);
    Jyd=zeros(Nx,Ny);
    Jmzd=zeros(Nx,Ny);
        
    % Static current densities
    Jxs=zeros(Nx,Ny);
    Jys=zeros(Nx,Ny);
    Jmzs=zeros(Nx,Ny);
    
    % Dispersive polarization and magnetization.
    Pxd=zeros(Nx,Ny);
    Pyd=zeros(Nx,Ny);
    Mzd=zeros(Nx,Ny);
    
    % Static Polarization and magnetizaton
    Pxs=zeros(Nx,Ny);
    Pys=zeros(Nx,Ny);
    Mzs=zeros(Nx,Ny);
    
    % Total polarization and currents
    
    Px=zeros(Nx,Ny);
    Py=zeros(Nx,Ny);
    Mz=zeros(Nx,Ny);
    Jmz=zeros(Nx,Ny);
    Jex=zeros(Nx,Ny);
    Jey=zeros(Nx,Ny);
    
%% Output Variables

    Hz_along_center=zeros(1,Nt);       % Store the maximum H field along center, used for plotting purposes

    % Matrix to store image of Hz
    Hz_field_fig=zeros(Nx,Ny,round(Nt/fig_count)+1);
    Hz_field_counter=1;

    % Enclosure/"system" variables
    system_momentum_x=zeros(1,Nt);
    system_momentum_y=zeros(1,Nt);
    % COM
    x_bar_system=zeros(1,Nt);
    y_bar_system=zeros(1,Nt);
    % Force
    F_system_x=zeros(1,Nt);
    F_system_y=zeros(1,Nt); 

    system_acceleration_x=zeros(1,Nt);    %[m/s^2]        Total acceleration
    system_velocity_x=zeros(1,Nt);        %[m/s]          Total Velocity        
    system_displacement_x=zeros(1,Nt);    %[m]            Displacement

    system_acceleration_y=zeros(1,Nt);    %[m/s^2]        Total acceleration      
    system_velocity_y=zeros(1,Nt);        %[m/s]          Total Velocity        
    system_displacement_y=zeros(1,Nt);    %[m]            Displacement

    x_bar_system_o=mean(x);
    x_bar_system_contribution=zeros(1,Nt);
    y_bar_system_o=mean(y);
    y_bar_system_contribution=zeros(1,Nt);

        


% Total variables

    Tx=zeros(Nx,Ny);                                % Tx component
    Ty=zeros(Nx,Ny);                                % Ty component

    g_mech_x=zeros(Nx,Ny);                         %Momentum density
    g_mech_y=zeros(Nx,Ny);                         %Momentum density

    x_bar_total=zeros(1,Nt);    
    y_bar_total=zeros(1,Nt);  
                
% Pulse_Output_Variables
    W=zeros(Nx,Ny);                     %[J/m^3]        Energy Density
    M_pulse=zeros(1,Nt);                %[Kg]           Total pulse mass
    pulse_energy=zeros(1,Nt);           %[J]            Total Energy

    % pulse x-variables
    pulse_momentum_x=zeros(1,Nt);       %[kg(m/s)]      Total Momentum
    G_x=zeros(Nx,Ny);                   %[Kg/(m^2*s)]   Momentum Density
    G_x_n_prev=0;
    x_bar_pulse=zeros(1,Nt);            %[m]            x- center of mass
    x_bar_pulse_contribution=zeros(1,Nt);%[m]per kg      % contribution of ceneter of mass

    % y-variables
    pulse_momentum_y=zeros(1,Nt);       %[kg(m/s)]      Total Momentum
    G_y=zeros(Nx,Ny);                   %[Kg/(m^2*s)]   Momentum Density
    G_y_n_prev=0;
    y_bar_pulse=zeros(1,Nt);            %[m]            x- center of mass
    y_bar_pulse_contribution=zeros(1,Nt);           %[m]per kg      % contribution of ceneter of mass

       
% Slab Output variables


% x-variables

    F_slab_x=zeros(1,Nt);
    F_slab_y=zeros(1,Nt);

    slab_momentum_x=zeros(1,Nt);        %[kg(m/s)]      Total Momentum          
    slab_acceleration_x=zeros(1,Nt);    %[m/s^2]        Total acceleration        
    slab_velocity_x=zeros(1,Nt);        %[m/s]          Total Velocity        
    slab_displacement_x=zeros(1,Nt);    %[m]            Displacement

    x_bar_slab_o=(1/(sum(sum(slab_mat))*dx*dy)).*sum(x(1:Nx).*sum((slab_mat(1:Nx,1:Ny))'.*dy)).*dx;

    x_bar_slab=zeros(1,Nt);             %[m]            x-center of mass
    x_bar_slab_contribution=zeros(1,Nt);%[m]per kg      % contribution of ceneter of mas


% y-variables        
    slab_momentum_y=zeros(1,Nt);        %[kg(m/s)]      Total Momentum          
    slab_acceleration_y=zeros(1,Nt);    %[m/s^2]        Total acceleration
    slab_velocity_y=zeros(1,Nt);        %[m/s]          Total Velocity
    slab_displacement_y=zeros(1,Nt);    %[m]            Displacement
    y_bar_slab=zeros(1,Nt);             %[m]            x-center of mass
    y_bar_slab_contribution=zeros(1,Nt);%[m]per kg      % contribution of ceneter of mas

    y_bar_slab_o=(1./(sum(sum(slab_mat))*dx*dy)).*sum(y(1:Ny).*sum((slab_mat(1:Nx,1:Ny)).*dy)).*dx;

 figure(2)
 subplot(2,1,1)
     title(' Permittivity ')
     surf(x*1e6,y*1e6,real(er_d)')
     shading flat
     view([0 90])
     colorbar east
     xlabel('x-axis ({\mu|m)')
     xlabel('y-axis ({\mu|m)')
     h_l=line([ x_bar_slab_o x_bar_slab_o]*1e6,[ y_bar_slab_o-L_x y_bar_slab_o+L_x]*1e6,[max(real(er_w)) max(real(er_w))],'color','black');
     line([ x_bar_slab_o-L_x x_bar_slab_o+L_x]*1e6,[ y_bar_slab_o y_bar_slab_o]*1e6,[max(real(er_w)) max(real(er_w))],'color','black')
     legend(h_l,'x_o slab,y_o slab')
     xlabel('x (mu m)')
     ylabel( 'y (mu m)')
 subplot(2,1,2) 
     title(' Drude Impedance ')
     surf(x*1e6,y*1e6,real(eta_d)')
     shading flat
     view([0 90])
    %  axis equal
     colorbar east
     h_l=line([ x_bar_slab_o x_bar_slab_o]*1e6,[ y_bar_slab_o-L_x y_bar_slab_o+L_x]*1e6,[max(real(er_w)) max(real(er_w))],'color','black');
     line([ x_bar_slab_o-L_x x_bar_slab_o+L_x]*1e6,[ y_bar_slab_o y_bar_slab_o]*1e6,[max(real(er_w)) max(real(er_w))],'color','black')
     legend(h_l,'x_o slab,y_o slab')
     xlabel('x (mu m)')
     ylabel( 'y (mu m)')
%% Update Coefficients
% Ex UPDATE COEFFCIENTS

C1=eps_o*(er(2:end,2:end)-1);

P1=(1/dt-gamma_e/2)./(1/dt+gamma_e/2);      
P2=eps_o*omega_e.^2./(2.*(1/dt+gamma_e/2));
den_x=(sigma_e/2+eps_o*er/dt+1/2.*P2);
X1=(eps_o*er/dt+sigma_e/2+P2/2);
X2=(eps_o*er/dt-sigma_e/2-P2/2);

X3=X2./X1;
X4=0.5.*((P1+1)./X1);
X5=(1./(dy.*X1));

P1=P1(2:end,2:end);
P2=P2(2:end,2:end);

X3=X3(2:end,2:end);
X4=X4(2:end,2:end);
X5=X5(2:end,2:end);

% Hz Update Coefficients

C2=mu_o*(mr-1);

M1=(1/dt-gamma_m./2)./(gamma_m./2+1/dt);                                    
M2=mu_o.*omega_m.^2./(2.*(gamma_m./2+1/dt));                                   

A1=(mu_o*mr./dt+sigma_m./2+1/2.*M2);
A2=(mu_o*mr./dt-sigma_m./2-1/2.*M2);

A3=A2./A1;

A4=(1./(2.*A1));
M3=M1+1;
A5=1./A1;

M1=M1(1:end-1,1:end-1);
M2=M2(1:end-1,1:end-1);
M3=M3(1:end-1,1:end-1);
A3=A3(1:end-1,1:end-1);
A4=A4(1:end-1,1:end-1);
A5=A5(1:end-1,1:end-1);


% Ey UPDATE COEFFICIENTS

P1Y=(1/dt-gamma_e/2)./(1/dt+gamma_e/2);
P2Y=eps_o*omega_e.^2./(2*(1/dt+gamma_e/2));
den_y=(sigma_e/2+eps_o*er/dt+1/2.*P2Y);

Y1=(eps_o*er/dt-sigma_e/2-P2Y/2)./den_y;
Y2=(P1Y+1)./(2*(den_y));
Y3=((1/dy)*(1./(den_y)));

P1Y=P1Y(2:end,2:end);
P2Y=P2Y(2:end,2:end);

Y1=Y1(2:end,2:end);
Y2=Y2(2:end,2:end);
Y3=Y3(2:end,2:end);

%% Indexes

i=4:Nx-4;
j=4:Ny-4;
switch source_type
case ('x')
i_W=[4:1:source1_x-4, source1_x+4:1:Nx-4];
j_W=j;
case ('y')
i_W=i;
j_W=[4:1:source1_y-4, source1_y+4:1:Nx-4];         
case ('xy') 
i_W=[4:1:source1_x-4, source1_x+4:1:Nx-4];
j_W=[4:1:source1_y-4, source1_y+4:1:Nx-4];
            
end




 clear gamma_e gamma_m  sigma_m sigma_e A1 A2 X1 X2 E omega_e omega_m
 clear FWHM E_source H_source E_y E_tot 
 clear N NPML N_y_gauss W_source b_x b_y den_x den_y
 clear eps_1 eps_2 er_1 er_Im_Si er_Re_Si  
 clear f_end f_start gauss_y_end gauss_y_start gauss_x_avg 
 clear gamma_e_material gamma_m_material mu_2 n_Si_exp
 clear k_Si_exp k_mat lam0 lambda_Si norm_G omega_e 
 clear m_fac omega_m 
 clear q s_d sh sig_t sig_y sigma_PML_e sigma_PML_max sim_cycles sim_time
 clear slope_x slope_y  sw t_avg
clear mr_d  er_d
    for n= 1:Nt
      

% Store Previous values

    Ex_n_prev=Ex(2:Nx,2:Ny);
    Dx_n_prev=Dx(2:Nx,2:Ny);   
    Jxd_n_prev=Jxd(2:Nx,2:Ny);      
    Jxs_n_prev=Jxs(2:Nx,2:Ny);      

    Ey_n_prev=Ey(2:Nx,2:Ny);
    Dy_n_prev=Ey(2:Nx,2:Ny);
    Jyd_n_prev=Jyd(2:end,2:end);   
    Jys_n_prev=Jys(2:end,2:end);   


    Hz_n_prev=Hz;
    Bz_n_prev=Bz;
    Jmz_n_prev_fr=Jmz;
    Jmzd_n_prev=Jmzd;    
              

%% Update Hz (n+1/2)

    E_term=-1*(1/dx).*(Ey(2:Nx,1:Ny-1)-Ey(1:Nx-1,1:Ny-1))...
            +(1/dy).*(Ex(1:Nx-1,2:Ny)-Ex(1:Nx-1,1:Ny-1));

    Hz(1:Nx-1,1:Ny-1)=  ((A3).*Hz(1:Nx-1,1:Ny-1)...    
                        -(A4).*((M3).*Jmzd(1:Nx-1,1:Ny-1))...
                        +(A5).*E_term);

    Jmzd(1:Nx-1,1:Ny-1)=M1.*Jmzd(1:Nx-1,1:Ny-1)+M2.*(Hz(1:Nx-1,1:Ny-1)...
                        +Hz_n_prev(1:Nx-1,1:Ny-1));            

    Jmzs(1:Nx-1,1:Ny-1)=-1*Jmzs(1:Nx-1,1:Ny-1)+(2/dt).*C2(1:Nx-1,1:Ny-1).*(Hz(1:Nx-1,1:Ny-1)-Hz_n_prev(1:Nx-1,1:Ny-1));
    Jmz(1:Nx-1,1:Ny-1)=Jmzs(1:Nx-1,1:Ny-1)+Jmzd(1:Nx-1,1:Ny-1);

        
    Mzd(1:Nx-1,1:Ny-1)=Mzd(1:Nx-1,1:Ny-1)+(dt/2).*(Jmzd(1:Nx-1,1:Ny-1)+Jmzd_n_prev(1:Nx-1,1:Ny-1));
    Mzs(1:Nx-1,1:Ny-1)=C2(1:Nx-1,1:Ny-1).*Hz(1:Nx-1,1:Ny-1);
    Mz=Mzd+Mzs;
    Bz=mu_o.*Hz+Mz;
             
    % Hz Source    
%     for q=1:length(G(:,1))   
%         % Y-GUASSIAN %
%     %Hz(source1_x,nyg1:nyg2)=gauss_y(nyg1:nyg2).*real(G(q,2)).*gauss_t(n).*a/eta_o.*cos(2*pi*G(q,1)*(t(n)))+ Hz(source1_x,nyg1:nyg2);
%         % X-GUASSIAN %
%      Hz(nxg1:nxg2,source1_y)=gauss_x(nxg1:nxg2)'.*real(G(q,2)).*gauss_t(n).*a/eta_o.*cos(2*pi*G(q,1)*(t(n)))+ Hz(nxg1:nxg2,source1_y);    
%     end    

    for q=1:length(G(:,1))
    switch source_type
        
    case ('x')
%Hz(source1_x,nyg1:nyg2)=gauss_y(nyg1:nyg2).*real(G(q,2)).*gauss_t(n).*a/eta_o.*cos(2*pi*G(q,1)*(t(n)))+ Hz(source1_x,nyg1:nyg2);            % X- Guassian %         
%Hz(source1_x,nyg1:nyg2)=gauss_y(nyg1:nyg2).*real(G(q,2)).*a/eta_o.*cos(2*pi*G(q,1)*(t(n)))+ Hz(source1_x,nyg1:nyg2);            % X- Guassian %         
Hz(source1_x,nyg1:nyg2)=gauss_t(n).*gauss_y(nyg1:nyg2).*a/eta_o.*cos(2*pi*f*(t(n)))+ Hz(source1_x,nyg1:nyg2);            % X- Guassian %         

    
    case ('y')
    
Hz(nxg1:nxg2,source1_y)=gauss_x(nxg1:nxg2)'.*real(G(q,2)).*gauss_t(n).*a/eta_o.*cos(2*pi*G(q,1)*(t(n)))+ Hz(nxg1:nxg2,source1_y);   
 
    case ('xy') 
                  
    Hz(nxg1:nxg2,source1_y)=gauss_x(nxg1:nxg2)'.*real(G(q,2)).*gauss_t(n).*a/eta_o.*cos(2*pi*G(q,1)*(t(n)))+ Hz(nxg1:nxg2,source1_y);    
    Hz(source1_x,nyg1:nyg2)=gauss_y(nyg1:nyg2).*real(G(q,2)).*gauss_t(n).*a/eta_o.*cos(2*pi*G(q,1)*(t(n)))+ Hz(source1_x,nyg1:nyg2);            % X- Guassian %
end

        
        
     
    end 



if (strcmp(model,'MN'))
    
    [W]=calculate_W_MN(c^2,i_W,j_W,Dx,Dx_n_prev,Dy,Dy_n_prev,Bz,Bz_n_prev,dx,dy,dt,W);


    Bz_at_y(i,j)=(1/4)*(Bz(i,j)+Bz(i-1,j)+Bz_n_prev(i,j)+Bz_n_prev(i-1,j));
%     Bz_at_x(i,j)=(1/4)*(Bz(i,j)+Bz(i,j-1)+Bz_n_prev(i,j)+Bz_n_prev(i,j-1));
    G_x_n_prev=G_x;
    G_x(i,j)=Dy(i,j).*Bz_at_y(i,j);                  % G_MN=DXB

    G_y_n_prev=G_y;
    Dx_av=(0.5).*(Dx(i,j)+Dx(i,j+1));
    G_y(i,j)=-1*Dx_av.*Bz(i,j);                  

    [Tx ] = Calculate_Tx_MN( i,j,Ex,Ey,Dx,Dy,Hz,Hz_n_prev,Bz,Bz_n_prev,dx,dy );
    [Ty,t1,t2,t3,t4 ] = Calculate_Ty_MN_2018( i,j,Ex,Ex_n_prev,Ey,Ey_n_prev,Dx,Dx_n_prev,...
                            Dy,Dy_n_prev,Hz,Hz_n_prev,Bz,Bz_n_prev,dx,dy );
end


if (strcmp(model,'Chu'))
    
    [W]=calculate_W_AB(i_W,j_W,Ex,Ex_n_prev,Ey,Ey_n_prev,Hz,Hz_n_prev,dx,dy,dt,W);  
    
    Hz_at_y(i,j)=(1/4)*(Hz(i,j)+Hz(i-1,j)+Hz_n_prev(i,j)+Hz_n_prev(i-1,j)); 
    Hz_at_x(i,j)=(1/4)*(Hz(i,j)+Hz(i,j-1)+Hz_n_prev(i,j)+Hz_n_prev(i,j-1));
    
    G_x_n_prev=G_x;        
    G_x(i,j)=eps_o*mu_o.*(Ey(i,j).*Hz_at_y(i,j)); % Chu, EL and AB  
    [ Tx] = Calculate_Tx_Chu(i,j,Ex,Ey,Dx,Dy,Hz,Hz_n_prev,Bz,Bz_n_prev,dx,dy );   
    
    G_y_n_prev=G_y;
    G_y(i,j)=-1*eps_o*mu_o.*(Ex(i,j).*Hz_at_x(i,j));                  % 
    [Ty,t1,t2,t3,t4 ] = Calculate_Ty_Chu( i,j,Ex,Ey,Dx,Dy,Hz,Hz_n_prev,Bz,Bz_n_prev,dx,dy );
   
end


if (strcmp(model,'EL'))
    
    [W]=calculate_W_AB(i_W,j_W,Ex,Ex_n_prev,Ey,Ey_n_prev,Hz,Hz_n_prev,dx,dy,dt,W);
    
    Hz_at_y(i,j)=(1/4)*(Hz(i,j)+Hz(i-1,j)+Hz_n_prev(i,j)+Hz_n_prev(i-1,j)); 
    Hz_at_x(i,j)=(1/4)*(Hz(i,j)+Hz(i,j-1)+Hz_n_prev(i,j)+Hz_n_prev(i,j-1));
    G_x_n_prev=G_x;        
    G_x(i,j)=eps_o*mu_o.*(Ey(i,j).*Hz_at_y(i,j)); % Chu, EL and AB    
    [Tx] = Calculate_Tx_EL(i,j,Ex,Ey,Dx,Dy,Hz,Hz_n_prev,Bz,Bz_n_prev,dx,dy );   
    
    
    G_y_n_prev=G_y;
    G_y(i,j)=-1*eps_o*mu_o.*(Ex(i,j).*Hz_at_x(i,j));                  % 
    [Ty,t1,t2,t3,t4 ] = Calculate_Ty_EL( i,j,Ex,Ey,Dx,Dy,Hz,Hz_n_prev,Bz,Bz_n_prev,dx,dy );           % 
    
    
    
    
end

if (strcmp(model,'AB'))
    
    [W]=calculate_W_AB(i_W,j_W,Ex,Ex_n_prev,Ey,Ey_n_prev,Hz,Hz_n_prev,dx,dy,dt,W);      
    Hz_at_y(i,j)=(1/4)*(Hz(i,j)+Hz(i-1,j)+Hz_n_prev(i,j)+Hz_n_prev(i-1,j)); 
    Hz_at_x(i,j)=(1/4)*(Hz(i,j)+Hz(i,j-1)+Hz_n_prev(i,j)+Hz_n_prev(i,j-1));
    G_x_n_prev=G_x;        
    G_x(i,j)=eps_o*mu_o.*(Ey(i,j).*Hz_at_y(i,j)); % Chu, EL and AB    
    [ Tx] = Calculate_Tx_AB(i,j,Ex,Ey,Dx,Dy,Hz,Hz_n_prev,Bz,Bz_n_prev,dx,dy );   
 
    G_y_n_prev=G_y;
    G_y(i,j)=-1*eps_o*mu_o.*(Ex(i,j).*Hz_at_x(i,j));                  % 
    [Ty,t1,t2,t3,t4 ] = Calculate_Ty_AB( i,j,Ex,Ey,Dx,Dy,Hz,Hz_n_prev,Bz,Bz_n_prev,dx,dy );
    
end

if (strcmp(model,'AMP'))
    
        [W]=calculate_W_MN(eps_o*c^2,i_W,j_W,Ex,Ex_n_prev,Ey,Ey_n_prev,Bz,Bz_n_prev,dx,dy,dt,W);    
        Bz_at_y(i,j)=(1/4)*(Bz(i,j)+Bz(i-1,j)+Bz_n_prev(i,j)+Bz_n_prev(i-1,j));   
        Bz_at_x(i,j)=(1/4)*(Bz(i,j)+Bz(i,j-1)+Bz_n_prev(i,j)+Bz_n_prev(i,j-1));
        G_x_n_prev=G_x;        
        G_x(i,j)=eps_o*(Ey(i,j).*Bz_at_y(i,j));
        [Tx] = Calculate_Tx_AMP(i,j,Ex,Ey,Dx,Dy,Hz,Hz_n_prev,Bz,Bz_n_prev,dx,dy );  
        G_y_n_prev=G_y;
        G_y(i,j)=-1*eps_o*(Ex(i,j).*Bz_at_x(i,j));                 % 
        [Ty,t1,t2,t3,t4 ] = Calculate_Ty_AMP( i,j,Ex,Ey,Dx,Dy,Hz,Hz_n_prev,Bz,Bz_n_prev,dx,dy );

 
end
% 
if n>round(7*T/dt)
     
    i_W=i;
    j_W=j;
    
end%
%         
         
%         % Calculate G_mech from Tx and G_EM
%         g_mech(i_f,j_f)=g_mech(i_f,j_f)-dt*Tx(i_f,j_f)-1*(G_x(i_f,j_f)-G_x_n_prev(i_f,j_f));
%         G_mech(1,n)=sum(sum(g_mech))*dx*dy;
 
%% Update Ex,Ey (n+1)

% Update Ex
    Ex(2:Nx,2:Ny)=(X3).*Ex(2:Nx,2:Ny)-X4.*Jxd(2:Nx,2:Ny)+X5.*(Hz(2:Nx,2:Ny)-Hz(2:Nx,1:Ny-1));
    
    Jxd(2:Nx,2:Ny)=P1.*Jxd(2:Nx,2:Ny)+P2.*(Ex(2:Nx,2:Ny)+Ex_n_prev);
    Pxd(2:Nx,2:Ny)=Pxd(2:Nx,2:Ny)+(dt/2).*(Jxd(2:Nx,2:Ny)+Jxd_n_prev);

    Jxs(2:Nx,2:Ny)=-Jxs(2:Nx,2:Ny)+(2/dt)*C1.*(Ex(2:Nx,2:Ny)-Ex_n_prev);
    Pxs(2:Nx,2:Ny)=Pxs(2:Nx,2:Ny)+(dt/2)*(Jxs(2:Nx,2:Ny)+Jxs_n_prev);

    Jex(2:Nx,2:Ny)=Jxs(2:Nx,2:Ny)+Jxd(2:Nx,2:Ny);    
    Px(2:Nx,2:Ny)=Pxd(2:Nx,2:Ny)+Pxs(2:Nx,2:Ny);
    Dx=eps_o.*Ex+Px;

% UPDATE Ey


    Ey(2:Nx,2:Ny)=Y1.*Ey(2:Nx,2:Ny)-Y2.*Jyd(2:Nx,2:Ny)-Y3.*(Hz(2:Nx,2:Ny)-Hz(1:Nx-1,2:Ny));

    Jyd(2:Nx,2:Ny)=P1Y.*Jyd(2:Nx,2:Ny)+P2Y.*(Ey(2:Nx,2:Ny)+Ey_n_prev);
    Pyd(2:end,2:end)=Pyd(2:end,2:end)+(dt/2)*(Jyd(2:end,2:end)+Jyd_n_prev);

    Jys(2:Nx,2:Ny)=-Jys(2:Nx,2:Ny)+(2/dt)*C1.*(Ey(2:Nx,2:Ny)-Ey_n_prev);
    Pys(2:end,2:end)=Pys(2:end,2:end)+(dt/2)*(Jys(2:end,2:end)+Jys_n_prev);

    Jey(2:Nx,2:Ny)=Jys(2:Nx,2:Ny)+Jyd(2:Nx,2:Ny);
    Py(2:end,2:end)=Pyd(2:end,2:end)+Pys(2:end,2:end);
    Dy=eps_o.*Ey+Py;
 
    % Source
%      for q=1:length(G(:,1))    
%          
%          % Y-Guassian %
%         % Ey(source1_x,nyg1:nyg2)=gauss_y(nyg1:nyg2).*real(G(q,2)).*gauss_t(n).*a.*cos(2*pi*G(q,1)*(t(n)))+ Ey(source1_x,nyg1:nyg2);            
%             % X- Guassian %
%             
%              Ex(nxg1:nxg2,source1_y)=gauss_x(nxg1:nxg2)'.*real(G(q,2)).*gauss_t(n).*a.*cos(2*pi*G(q,1)*(t(n)))+ Ex(nxg1:nxg2,source1_y); 
%      
%      end
    % Source
     for q=1:length(G(:,1))    
          switch source_type
%          
              case ('x')
            % Ey(source1_x,nyg1:nyg2)=gauss_y(nyg1:nyg2).*real(G(q,2)).*gauss_t(n).*a.*cos(2*pi*G(q,1)*(t(n)))+ Ey(source1_x,nyg1:nyg2);            
%            
             Ey(source1_x,nyg1:nyg2)=gauss_t(n).*gauss_y(nyg1:nyg2).*a.*cos(2*pi*f*(t(n)))+ Ey(source1_x,nyg1:nyg2);            

% X- Guassian %
%             
             case ('y')
             Ex(nxg1:nxg2,source1_y)=gauss_x(nxg1:nxg2)'.*real(G(q,2)).*gauss_t(n).*a.*cos(2*pi*G(q,1)*(t(n)))+ Ex(nxg1:nxg2,source1_y); 
%             
              case ('xy') 
                 
            Ey(source1_x,nyg1:nyg2)=gauss_y(nyg1:nyg2).*real(G(q,2)).*gauss_t(n).*a.*cos(2*pi*G(q,1)*(t(n)))+ Ey(source1_x,nyg1:nyg2);            
            Ex(nxg1:nxg2,source1_y)=gauss_x(nxg1:nxg2)'.*real(G(q,2)).*gauss_t(n).*a.*cos(2*pi*G(q,1)*(t(n)))+ Ex(nxg1:nxg2,source1_y); 
          end
     end

% if n<round(5*T/dt)
%     W=0.*W;
%     
% end


   
%     y_bar_calc(1,j)=y(j).*sum((W(:,j)./c^2).*dx);
%     


     
%% Pulse Calculations
    M_pulse(1,n)=sum(sum(W./c^2))*dx*dy;
    pulse_energy(1,n)=sum(sum(W))*dx*dy;   
           
    pulse_momentum_x(1,n)=sum(sum(G_x(i,j)))*dx*dy;
    pulse_momentum_y(1,n)=sum(sum(G_y(i,j)))*dx*dy;
    x_bar_pulse(1,n)=(1/M_pulse(1,n)).*sum(x(i).*sum((W(i,j)./c^2)'.*dy)).*dx;
    y_bar_pulse(1,n)=(1/M_pulse(1,n)).*sum(y(j).*sum((W(i,j)./c^2).*dx)).*dy;
%% Slab Calculations

%     fx_T(i,j)=-Tx(i,j)-1*(1/dt)*(G_x(i,j)-G_x_n_prev(i,j));
%     fy_T(i,j)=-Ty(i,j)-1*(1/dt)*(G_y(i,j)-G_y_n_prev(i,j));
    
    g_mech_x(i,j)=g_mech_x(i,j)-dt.*Tx(i,j)-1*(G_x(i,j)-G_x_n_prev(i,j));
    g_mech_y(i,j)=g_mech_y(i,j)-dt.*Ty(i,j)-1*(G_y(i,j)-G_y_n_prev(i,j));
  
 
    
%  
%         slab_momentum_x(1,n)=sum(sum(g_mech_x(i_f,j_f)))*dx*dy;
%     system_momentum_x(1,n)=sum(sum(g_mech_x(i,j)))*dx*dy-slab_momentum_x(1,n);
%     
%     slab_momentum_y(1,n)=sum(sum(g_mech_y(i_f,j_f)))*dx*dy;
%     system_momentum_y(1,n)=sum(sum(g_mech_y(i,j)))*dx*dy-slab_momentum_y(1,n);
% %     
        slab_momentum_x(1,n)=sum(sum(g_mech_x(i,j).*slab_mat(i,j)))*dx*dy;
   % system_momentum_x(1,n)=sum(sum(g_mech_x(i,j)))*dx*dy-slab_momentum_x(1,n);
       system_momentum_x(1,n)=sum(sum(g_mech_x(i,j).*system_mat(i,j)))*dx*dy;
    
    
slab_momentum_y(1,n)=sum(sum(g_mech_y(i,j).*slab_mat(i,j)))*dx*dy;
% system_momentum_y(1,n)=sum(sum(g_mech_y(i,j)))*dx*dy-slab_momentum_y(1,n);
system_momentum_y(1,n)=sum(sum(g_mech_y(i,j).*system_mat(i,j)))*dx*dy;
   
    
    
    if n>1
       
F_slab_x(1,n)=(slab_momentum_x(1,n)-slab_momentum_x(1,n-1))/dt;
slab_acceleration_x(1,n)=1/M_slab*F_slab_x(1,n);     
slab_velocity_x(1,n)=sum(slab_acceleration_x(1,1:n))*dt;
slab_displacement_x(1,n)=sum(slab_velocity_x(1,1:n))*dt;
x_bar_slab(1,n)=x_bar_slab_o+slab_displacement_x(1,n);
 
F_slab_y(1,n)=(slab_momentum_y(1,n)-slab_momentum_y(1,n-1))/dt;
slab_acceleration_y(1,n)=1/M_slab*F_slab_y(1,n);
slab_velocity_y(1,n)=sum(slab_acceleration_y(1,1:n))*dt;
slab_displacement_y(1,n)=sum(slab_velocity_y(1,1:n))*dt;

y_bar_slab(1,n)=y_bar_slab_o+slab_displacement_y(1,n);

F_system_x(1,n)=(system_momentum_x(1,n)-system_momentum_x(1,n-1))/dt;
system_acceleration_x(1,n)=1/M_system*F_system_x(1,n);     
system_velocity_x(1,n)=sum(system_acceleration_x(1,1:n))*dt;
system_displacement_x(1,n)=sum(system_velocity_x(1,1:n))*dt;
x_bar_system(1,n)=x_bar_system_o+system_displacement_x(1,n);
 
F_system_y(1,n)=(system_momentum_y(1,n)-system_momentum_y(1,n-1))/dt;
system_acceleration_y(1,n)=1/M_system*F_system_y(1,n);
system_velocity_y(1,n)=sum(system_acceleration_y(1,1:n))*dt;
system_displacement_y(1,n)=sum(system_velocity_y(1,1:n))*dt;

y_bar_system(1,n)=y_bar_system_o+system_displacement_y(1,n);

    
    
    
    end
%     slab_momentum_x(1,n)=slab_velocity_x(1,n).*M_slab;

% slab_energy(1,n)=1/2*M_slab.*slab_velocity_x(1,n);
% work_done_x(1,n)=sum((slab_acceleration_x(1:n)*M_slab).*slab_displacement_x(1:n))*dt;
% x-components
	
% y-components    

 
    
%% System Calculations
M_total=M_slab+M_pulse(1,n)+M_system;
x_bar_total(1,n)= (1/M_total)*(M_slab.*x_bar_slab(1,n)+M_pulse(1,n).*x_bar_pulse(1,n)+M_system*x_bar_system(1,n));
x_bar_pulse_contribution(1,n)=((1/M_total)*(M_pulse(1,n).*x_bar_pulse(1,n)));
x_bar_slab_contribution(1,n)=(1/M_total)*(M_slab.*x_bar_slab(1,n));
 x_bar_system_contribution(1,n)=(1/M_total)*(M_system.*x_bar_system(1,n));
% Y Center of Mass

y_bar_total(1,n)= (1/M_total)*(M_slab.*y_bar_slab(1,n)+M_pulse(1,n).*y_bar_pulse(1,n)+M_system*y_bar_system(1,n));
y_bar_pulse_contribution(1,n)=(1/M_total)*(M_pulse(1,n).*y_bar_pulse(1,n));
y_bar_slab_contribution(1,n)=(1/M_total)*(M_slab.*y_bar_slab(1,n));
 y_bar_system_contribution(1,n)=(1/M_total)*(M_system.*y_bar_system(1,n)); 
 
% 
%         PML_momentum_x_flux(1,n)=   sum(sum(G_x(i(1):i(2),j)))*dy*dx...         % Left PML
%                                     +sum(sum(G_x(i(end-1):i(end),j)))*dy*dy...         % Right PML
%                                     +sum(sum(G_x(i(2:end-1),j(1):j(2))))*dx...         % Bottom PML
%                                     +sum(sum(G_x(i(2:end-1),j(end-1):j(end))))*dx*dy;         % Bottom PML
%                             
%         PML_momentum_x(1,n)=sum(PML_momentum_x_flux)*dt;
 
%% Plots




if n_count==fig_count
    
    
    Hz_field_fig(:,:,Hz_field_counter)=Hz;
Hz_field_counter= Hz_field_counter+1; 
% cmax=max(max(max(Hz_field_fig)));

Hz_along_center(1,n)=max(Hz(NPML_x:round(Nx/2),round(Ny/2)));
cmax=max(Hz_along_center);

if (plot_on==1)


% Plot pule interacting sith slab
figure(3)
    title(model)
    surf(x*1e6,y*1e6,Hz')

       hold on
       contour(x*1e6,y*1e6,er'-er_slab,1,'color','black')
       hold off
       xlabel('x axis [{\mu}m]')
       ylabel('y axis [{\mu}m]')
     c_lim=    [-cmax cmax];  
       caxis(c_lim)
       shading flat
% Draw PML limes
        line([NPML_x*dx  NPML_x*dx]*1e6, [(NPML_y)*dy (Ny-NPML_y)*dy]*1e6,'color','black');  
        line([(Nx-NPML_x)*dx (Nx-NPML_x)*dx]*1e6, [(NPML_y)*dy (Ny-NPML_y)*dy]*1e6,'color','black');  
        line([NPML_x*dx (Nx-NPML_x)*dx]*1e6, [(NPML_y)*dy (NPML_y)*dy]*1e6,'color','black');  
        line([NPML_x*dx (Nx-NPML_x)*dx]*1e6, [(Ny-NPML_y)*dy (Ny-NPML_y)*dy]*1e6,'color','black'); %         

 % Draw force integration locations
%  axis equal
view([0 90])


% X momentum figure
figure(4)
if(t(n)>1.5*tau)
P_o=(max(pulse_energy)./c);
clf
else
P_o=1;
end
plot(t(1:n),pulse_momentum_x(1:n)./P_o,'color','r')
hold on
plot(t(1:n),slab_momentum_x(1:n)./P_o,'color','b')
hold on
plot(t(1:n),system_momentum_x(1:n)./P_o,'--','color','black')
hold on
plot(t(1:n),(slab_momentum_x(1:n)+pulse_momentum_x(1:n)+system_momentum_x(1:n))./P_o,'color','g')


xlabel('time(s)')
ylabel(' X-Linear Momentum [kgm/s]')
legend('Pulse','Object','System','Total','Location','northwest')


% Y momentum figure
figure(5)
plot(t(1:n),pulse_momentum_y(1:n)./P_o,'color','r')
hold on
plot(t(1:n),slab_momentum_y(1:n)./P_o,'color','b')
hold on
plot(t(1:n),system_momentum_y(1:n)./P_o,'--','color','black')
hold on
plot(t(1:n),(slab_momentum_y(1:n)+pulse_momentum_y(1:n)+system_momentum_y(1:n))./P_o,'color','green')

xlabel('time(s)')
ylabel(' Y-Linear Momentum [kgm/s]')
legend('Pulse','Object','System','Total','Location','northwest')

% X Center of Mass figure
figure(6)
M_tot=M_pulse+M_slab+M_system;
xb_pulse_displacement(5:n)=(M_pulse(5:n)./M_tot(5:n)).*(x_bar_pulse(5:n)-x_bar_pulse(1,5));
xb_slab_displacement(5:n)=(M_slab./M_tot(5:n)).*(x_bar_slab(5:n)-x_bar_slab(1,5));
xb_system_displacement(5:n)=((M_system)./M_tot(5:n)).*(x_bar_system(5:n)-x_bar_system(1,5));
xb_total_displacement=xb_pulse_displacement...
                            +xb_slab_displacement...
                            +xb_system_displacement;

    h_1=plot(t(5:n)/fs,xb_pulse_displacement(5:n),'color','r');
hold on
    h_2=plot(t(5:n)/fs,xb_slab_displacement(5:n),'color','b');
hold on
    plot(t(5:n)/fs,xb_system_displacement(5:n),'--','color','black')
hold on
    h_3=plot(t(5:n)/fs,(xb_total_displacement(5:n)),'color','g');
h_vec=[h_1,h_2,h_3];
%     set(h_vec,'linewidth',Plw)
h_4= legend('Pulse','Slab','Enclosure','Total');
set(h_4,'location','northwest')
% set(h_4,'fontsize',legend_fontsize)

xlabel(' Time (fs) ')
ylabel(' %M \cdot \Delta{x} ({\mu}m) ')

% % plot(t(2:n),x_bar_pulse_contribution(2:n)+x_bar_slab_contribution(1,2),'color','r')
% % hold on
% % plot(t(2:n),x_bar_slab_contribution(2:n),'color','blue')
% % hold on
% % plot(t(2:n),x_bar_system_contribution(2:n),'color','black')
% % legend(' pulse contribution + z_o ','slab contribution',' Total' ,'Location','northwest')
% % xlabel('Time (s)')
% % ylabel('X-Center of Mass location (m)')
% % % hold on
% % % 
% % % plot(t(1:n),x_bar_total(1:n),'color','black')
% % % legend(' pulse contribution','slab contribution','System Contribution' ,'Location','southwest')
% % xlabel('Time (s)')
% % ylabel('X-Center of Mass location (m)')
% % % subplot(1,2,1)
% % %     plot(x,W(:,round(Ny/2)))
% % %     xlabel('x-axis (m) ')
% % %     ylabel(' W along cneter [J/m] ')
% % % subplot(1,2,2)
% % %     plot(t(1:n),pulse_energy(1:n),'color','r')
% % %     xlabel(' time(s) ')
% % %     ylabel(' Total Pulse Energy [J] ')

% Y center of mass figure
figure(7)
plot(t(1:n),y_bar_pulse_contribution(1:n)+y_bar_slab_o,'color','r')
hold on
plot(t(1:n),y_bar_slab_contribution(1:n),'color','blue')
hold on
plot(t(1:n),y_bar_system_contribution(1:n),'--','color','black')
hold on
plot(t(1:n),y_bar_total(1:n),'color','green')
legend('pulse contribution + z_o ','slab contribution','Enclosusure','Total','Location','southwest')
xlabel(' Time (s)')
ylabel('Y-Center of Mass location (m)')

figure(8)
plot(t(1:n),pulse_energy(1:n),'color','r')
xlabel(' Time (s) ')
ylabel(' Total Pulse Energy (J) ')
% % 
figure(9)
subplot(1,2,1)
plot(t(1:n),slab_displacement_x(1:n).*1e9)
xlabel(' Time (s) ')
ylabel (' {\Delta}{x} (nm) ')
subplot(1,2,2)
plot(t(1:n),slab_displacement_y(1:n).*1e9,'color','g')
xlabel(' Time (s) ')
ylabel (' {\Delta}{y} (nm) ')
% % subplot(1,2,1)
% % plot(t(1:n),slab_displacement_x(1:n))
% % hold on
% % plot(t(1:n),slab_displacement_y(1:n),'--','color','g')
% % xlabel(' time(s) ')
% % ylabel(' Displacement (m) ')
% % 
% % subplot(1,2,2)
% % plot(t(1:n),slab_velocity_x(1:n),'color','r')
% % xlabel(' time(s) ')
% % ylabel(' x-velocity (m/s) ')
% %   
% end
drawnow

end
    
end

if n_count==fig_count
    n_count=0;
    
end
% 
% t(n);
n_count=n_count+1;
% x_test=x_bar_pulse(1,n)*1e6
    end
    

if save_mode==1
 model_type=char(models(sim_j));
%  save(strcat(model_type,'_Workspace_data','-',month,'-',day,'-',hour,'-',minutes));
 save(strcat(pwd,'\output_data\','Hz_field_',sim_case),'Hz_field_fig')
 clear Hz_field_fig
save(strcat(pwd,'\output_data\',model_type,'_Workspace_data','_',sim_case));
 
end
close all
end


