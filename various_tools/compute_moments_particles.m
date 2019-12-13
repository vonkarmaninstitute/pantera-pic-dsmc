close all
clear
clc

page_screen_output(0);

% This script loads particle files in the dumps folder and computes moments

files_list = dir('../dumps/proc_00000_time_*');
% files_list = dir('/media/starlight/Maxtor/PANTERA_data/proc_00000_time_*');

M  = 6.63E-26; % [kg] Molecular mass
nb = 1.4e20;     % Number density
dt = 5e-8;    % Timestep

t_vect = zeros(numel(files_list), 1);

for ii = 1:numel(files_list)

    filename = ['../dumps/', files_list(ii).name];
    % filename = ['/media/starlight/Maxtor/PANTERA_data/', files_list(ii).name];

    fprintf('Loading file: %s...\n', filename)
    dd = load(filename);    

    fprintf('  Now computing moments...\n')
   
    % Time vector
    timestep   = dd(1,1);
    t_vect(ii) = dt*timestep;

    % Extract particles 
   
    vx_p = dd(:,5);
    vy_p = dd(:,6);
    vz_p = dd(:,7);
    
    % % Set variable names for compatibility (I copied the script from another project)
    % % VR_range = vx_particles;
    % % Vz_range = vy_particles;
    % % Vt_range = vz_particles;
    % % 
    % % mass_spec = M;
    
    % ++++++  NUMBER OF PARTICLES  ++++++
    num_part = numel(vx_p);
    
    % ++++++  NUMBER DENSITY  ++++++
    % Average density is imposed in this problem
    mom_n = nb; % [1/m^3] number density
    
    % ++++++  VELOCITY  ++++++
    mom_ux = sum(vx_p)/num_part
    mom_uy = sum(vy_p)/num_part
    mom_uz = sum(vz_p)/num_part
    
    mom_u2 = sum(vx_p.^2 + vy_p.^2 + vy_p.^2)/num_part;
    
    % ###### some commodity vectors #####
    Cx = vx_p - mom_ux;
    Cy = vy_p - mom_uy;
    Cz = vz_p - mom_uz;
    
    C2 = Cx.^2 + Cy.^2 + Cz.^2;
    
    % % Load ionization cross section
    % V2_range = Vz_range.^2 + Vt_range.^2 + VR_range.^2;
    % V_range  = sqrt(V2_range);
    % EE_range = 1/2*9.109e-31*V2_range;
    % 
    % sig_data_eV = load('sigma_ionization.dat');
    % E_eV_data   = sig_data_eV(:,1);
    % sig_m2_data = sig_data_eV(:,2);
    % E_J_data    = E_eV_data*1.602e-19;
    % 
    % sig = interp1(E_J_data, sig_m2_data, EE_range);
    % 
    % mom_kf = mean(sig.*V_range);
    
    % ++++++  TEMPERATURES  ++++++
    mom_Tx = 1/(num_part*1.38e-23)*M*sum(Cx.^2);
    mom_Ty = 1/(num_part*1.38e-23)*M*sum(Cy.^2);
    mom_Tz = 1/(num_part*1.38e-23)*M*sum(Cz.^2);
    mom_T  = 1/(3*num_part*1.38e-23)*M*sum(C2);
    
    % ++++++  PRESSURE TENSOR  ++++++
    % Note that the tensor is symmetric by construction 
    % (by kinetic definition), therefore it only has 6
    % independent components
    mom_Pxx = 1/num_part*M*mom_n*sum(Cx.^2)
    mom_Pxy = 1/num_part*M*mom_n*sum(Cx.*Cy)
    mom_Pxz = 1/num_part*M*mom_n*sum(Cx.*Cz)
    mom_Pyy = 1/num_part*M*mom_n*sum(Cy.^2)
    mom_Pyz = 1/num_part*M*mom_n*sum(Cy.*Cz)
    mom_Pzz = 1/num_part*M*mom_n*sum(Cz.^2)
    
    % ++++++  HEAT FLUX TENSOR  ++++++
    % Qijk is a third order tensor, symmetric with respect to permutations 
    % of the indices, therefore it has only 10 independent components
    mom_Qxxx = 1/num_part*M*mom_n*sum(Cx.^3)
    mom_Qxxy = 1/num_part*M*mom_n*sum(Cx.*Cx.*Cy) 
    mom_Qxyy = 1/num_part*M*mom_n*sum(Cx.*Cy.*Cy) 
    mom_Qyyy = 1/num_part*M*mom_n*sum(Cy.^3) 
    mom_Qyyz = 1/num_part*M*mom_n*sum(Cy.*Cy.*Cz) 
    mom_Qyzz = 1/num_part*M*mom_n*sum(Cy.*Cz.*Cz) 
    mom_Qzzz = 1/num_part*M*mom_n*sum(Cz.^3) 
    mom_Qxxz = 1/num_part*M*mom_n*sum(Cx.*Cx.*Cz)
    mom_Qxzz = 1/num_part*M*mom_n*sum(Cx.*Cz.*Cz)
    mom_Qxyz = 1/num_part*M*mom_n*sum(Cx.*Cy.*Cz) 
    
    % ++++++  HEAT FLUX VECTOR  ++++++
    mom_qx = 1/num_part*M*mom_n*sum(C2.*Cx)
    mom_qy = 1/num_part*M*mom_n*sum(C2.*Cy)
    mom_qz = 1/num_part*M*mom_n*sum(C2.*Cz)
    
    % ++++++  Rijkl  ++++++
    % I have some contractions: Riijj = \sum_i \sum_j
    mom_Riijj = 1/num_part*M*mom_n*sum(C2.*C2)
    
    mom_Rxxjj = 1/num_part*M*mom_n*sum(Cx.*Cx.*C2)
    mom_Rxyjj = 1/num_part*M*mom_n*sum(Cx.*Cy.*C2)
    mom_Rxzjj = 1/num_part*M*mom_n*sum(Cx.*Cz.*C2)
    mom_Ryyjj = 1/num_part*M*mom_n*sum(Cy.*Cy.*C2)
    mom_Ryzjj = 1/num_part*M*mom_n*sum(Cy.*Cz.*C2)
    mom_Rzzjj = 1/num_part*M*mom_n*sum(Cz.*Cz.*C2)
    
    mom_Ryyyy = 1/num_part*M*mom_n*sum(Cy.*Cy.*Cy.*Cy)
    mom_Ryyyz = 1/num_part*M*mom_n*sum(Cy.*Cy.*Cy.*Cz)
    mom_Ryyzz = 1/num_part*M*mom_n*sum(Cy.*Cy.*Cz.*Cz)
    mom_Rzzyz = 1/num_part*M*mom_n*sum(Cz.*Cz.*Cy.*Cz)
    mom_Ryyxy = 1/num_part*M*mom_n*sum(Cy.*Cy.*Cx.*Cy)
    mom_Rxxyy = 1/num_part*M*mom_n*sum(Cx.*Cx.*Cy.*Cy)
    mom_Rxxxy = 1/num_part*M*mom_n*sum(Cx.*Cx.*Cx.*Cy)
    mom_Rzzzz = 1/num_part*M*mom_n*sum(Cz.*Cz.*Cz.*Cz)
    mom_Rzzxz = 1/num_part*M*mom_n*sum(Cz.*Cz.*Cx.*Cz)
    mom_Rxxzz = 1/num_part*M*mom_n*sum(Cx.*Cx.*Cz.*Cz)
    mom_Rxxxz = 1/num_part*M*mom_n*sum(Cx.*Cx.*Cx.*Cz)
    mom_Rxxxx = 1/num_part*M*mom_n*sum(Cx.*Cx.*Cx.*Cx)
    mom_Ryyxz = 1/num_part*M*mom_n*sum(Cy.*Cy.*Cx.*Cz)
    mom_Rzzxy = 1/num_part*M*mom_n*sum(Cz.*Cz.*Cx.*Cy)
    mom_Rxxyz = 1/num_part*M*mom_n*sum(Cx.*Cx.*Cy.*Cz)
    
    % ++++++  Siijjk +++++++
    mom_Siijjx = 1/num_part*M*mom_n*sum(C2.*C2.*Cx)
    mom_Siijjy = 1/num_part*M*mom_n*sum(C2.*C2.*Cy)
    mom_Siijjz = 1/num_part*M*mom_n*sum(C2.*C2.*Cz)
    
    mom_Sjjyyy = 1/num_part*M*mom_n*sum(C2.*Cy.*Cy.*Cy)
    mom_Sjjyyz = 1/num_part*M*mom_n*sum(C2.*Cy.*Cy.*Cz)
    mom_Sjjyzz = 1/num_part*M*mom_n*sum(C2.*Cy.*Cz.*Cz)
    mom_Sjjxyy = 1/num_part*M*mom_n*sum(C2.*Cx.*Cy.*Cy)
    mom_Sjjxxy = 1/num_part*M*mom_n*sum(C2.*Cx.*Cx.*Cy)
    mom_Sjjxyz = 1/num_part*M*mom_n*sum(C2.*Cx.*Cy.*Cz)
    mom_Sjjzzz = 1/num_part*M*mom_n*sum(C2.*Cz.*Cz.*Cz)
    mom_Sjjxzz = 1/num_part*M*mom_n*sum(C2.*Cx.*Cz.*Cz)
    mom_Sjjxxz = 1/num_part*M*mom_n*sum(C2.*Cx.*Cx.*Cz)
    mom_Sjjxxx = 1/num_part*M*mom_n*sum(C2.*Cx.*Cx.*Cx)
    
    % Save moments into vector
    U_14(1,ii)  = mom_n*M;
    U_14(2,ii)  = mom_ux;
    U_14(3,ii)  = mom_uy;
    U_14(4,ii)  = mom_uz;
    U_14(5,ii)  = mom_Pxx;
    U_14(6,ii)  = mom_Pxy;
    U_14(7,ii)  = mom_Pxz;
    U_14(8,ii)  = mom_Pyy;
    U_14(9,ii)  = mom_Pyz;
    U_14(10,ii) = mom_Pzz;
    U_14(11,ii) = mom_qx;
    U_14(12,ii) = mom_qy;
    U_14(13,ii) = mom_qz;
    U_14(14,ii) = mom_Riijj;

end

% Save stuff
save sol_pantera_test U_14 t_vect

figure
plot(t_vect, U_14(1,:)', '-x', 'linewidth', 2)
xlabel('Time [s]')
ylabel('Density [kg/m3]')

figure
plot(t_vect, U_14(2,:)', '-xb', 'linewidth', 2)
hold on
plot(t_vect, U_14(3,:)', '-xr', 'linewidth', 2)
plot(t_vect, U_14(4,:)', '-xg', 'linewidth', 2)
xlabel('Time [s]')
ylabel('Velocities [m/s]')

figure
plot(t_vect, U_14(5:10,:)', '-x', 'linewidth', 2)
hold on
xlabel('Time [s]')
ylabel('Pressure tensor components [Pa]')


