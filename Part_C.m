clc;clear all;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% Author: Noah Sadaka
% Date: February 17 2019 - April 1 2019
% Course: Gas Turbine Design

% Purpose: Solve for the meanline design of a Turbine based on set ranges
% and criteria.

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


% Stage Nomenclature
% Station 1 = Vane Inlet
% Station 2 = Vane Exit and Blade Inlet
% Station 3 = Blade Exit

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%% Interpolation Data Import %%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

stagger_data = csvread('stagger_angle.csv',1); % Stagger angle data
tmaxc_data = csvread('tmaxc_v_betas.csv',1); % max thickness over chord
Yp_beta_0_data = csvread('Yp_beta1_0.csv',1); % Profile loss coeff beta=0
Yp_beta_alpha_data = csvread('Yp_beta1_alfa1.csv',1); % Profile loss coeff beta=alpha
delta_phi_data = csvread('fig_14_enecoef.csv',1); % trailing edge energy coefficient
kin_visc = csvread('kinematic_viscosity.csv',1); % kinematic viscosity with temperature

%%%%%%%%%%%%%%%%%%%%%%
%%% INITIALIZATION %%%
%%%%%%%%%%%%%%%%%%%%%%


% RAYMOND: values to put in the input data csv table:
% alpha_3, M_3, R, AN2, U_h (U hub), the inc at blade inlet
% zweif blade and vane
% Thanks!

eta_i = 0.85; % desired efficiency

% Cycle Analysis Values
To_1 = 1147.98; % Turbine Inlet Temperature [K]
Po_1 = 315940.86; % Turbine Inlet Pressure [Pa]
m_1 = 4.895; % Turbine Inlet Massflow [kg/s]

To_2 = 1126.67; % Total temperature after bleed air addition [K]

m_2 = 5; % Mass flow through blade [kg/s]

To_3 = 974.42; % Total temperature after station [K]
Po_3 = 151212.2; % Total pressure exiting blade, before ITD [Pa]
m_3 = m_2; % Mass flow exiting blade [kg/s]

W = 883830.57; % Stage work [W]

% Gas Constants and Coefficients
Cp = 1148; % Gas specific heat capacity [kg/JK]
gamma = 1.333; % Gas Gamma
Rg = 287; % Gas Constant
% 
% alpha_3_inp = -5:.1:5;
% M_3_imp = 0.3:0.01:0.45;
% R_imp = 0.3:0.01:0.4;
% AN2_mult = 0.9:0.1:1;
% U_mult = 0.9:0.01:1;
% 
% for alf_ind = 1:length(alpha_3_inp)
%     for m_ind = 1:length(M_3_imp)
%         for r_ind = 1:length(R_imp)
%             for U_ind = 1:length(U_mult)
%                 for AN2_ind = 1:length(AN2_mult)
error = 0; % flag for errors
error_p = 0; % flag for pressure errors


% % Initial Conditions
%                 alpha_3 = alpha_3_inp(alf_ind); % Blade exit swirl angle [deg] Range: -5 to 30
%                 M_3 = M_3_imp(m_ind); % Blade exit mach number. Range: 0.3-0.45
%                 R = R_imp(r_ind); % Reaction at the meanline.
alpha_3 = 24;
M_3 = .32;
R = .45;
inc_1_des = 0; % design incidence [deg]
inc_2_des = 3; % design incidence
% U_h = U_mult(U_ind)*1100*0.3048; % Max Blade speed at hub [m/s]
% U_h=    .85*1100*.3048;
U_h=295.05;
% AN2 = AN2_mult(AN2_ind)*4.5E10; % AN2 [in2 rpm^2]
AN2 = 3.6E10;
zweif_vane = 0.89; % Range: 0.7-0.8
zweif_blade = 0.91; % Range: 0.85-0.95
blade_tip_clearance = 0.009;

% Given Vane Variables
AR_v = 0.7; % Vane Aspect Ratio
TE_v = 0.045*0.0254; % Vane TE Thickness [m]

% Given Blade Variables
AR_b = 1.45; % Blade Aspect Ratio
TE_b = 0.025*0.0254; % Blade TE Thickness [m]
% for i=1:1:1;
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%% MEANLINE CALCULATIONS %%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% Station 1 Velocity Triangle

M_1 = 0.125; % Turbine Inlet Mach Number
alpha_1 = -10; % Inlet Swirl Angle [deg]
T_1 = To_1/(1+0.5*(gamma-1)*M_1^2); % Temperature at Vane [K]
P_1 = Po_1*(T_1/To_1)^(gamma/(gamma-1)); % Pressure at Vane [Pa]
rho_1 = P_1/(Rg*T_1); % Density at Vane Inlet [kg/m3]
V_1 = M_1 * sqrt(gamma * Rg * T_1); % Vane Inlet Velocity
Va_1 = V_1*cosd(alpha_1); % Vane Inlet Axial Velocity
Vu_1 = V_1 * sind(alpha_1); % Vane Inlet Tangential Velocity
A_1 = m_1/(rho_1 * Va_1); % Annulus Area at station 1 [m^2]


% Station 3 Absolute Velocity Triangle
T_3 = To_3 / (1 + 0.5 * (gamma - 1) * M_3^2); % Temperature at Blade Exit [K]
P_3 = Po_3 * (T_3 / To_3)^(gamma / (gamma - 1)); % Pressure at Blade Exit [Pa]
rho_3 = P_3 / (Rg * T_3); % Density at Blade Exit [kg/m3]
V_3 = M_3 * sqrt(gamma * Rg * T_3); % Blade Exit Inlet Velocity [m/s]
Va_3 = V_3 * cosd(alpha_3); % Axial velocity at station 3 [m/s]
A_3 = m_3 / (rho_3 * Va_3); % Annulus Area at station 3 [m^2]
Vu_3 =  sqrt(V_3^2 - Va_3^2); % Swirl velocity [m/s]

% Station 2 Absolute Velocity Triangle (part 1)
Ws = W/m_2; % Specific work [W/kg]
A_2 = A_3; % Assumption: rotor area is constant
T_2 = R * (T_1 - T_3) + T_3; % Temperature at station 2 [K]
V_2 = sqrt(2*Cp*(To_2 - T_2)); % Absolute Velocity [m/s]
M_2 = V_2/sqrt(T_2 * gamma * Rg); % Mach


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%% RADII AND STRUCTURAL LIMITATIONS %%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% AN2_max = 4.5E10;
A_rpm = 1550*A_2; % A_2 in in2
N_rpm = sqrt(AN2/A_rpm); % Rotation speed [rpm]
% N_rpm = 22524;
% N_rpm_max = sqrt(AN2_max/A_rpm);
N_rads = N_rpm * (1/60) * (2*pi); % Rotation speed [rad/s]

r_h_2 = U_h/N_rads; % hub radius [m]

r_t_2 = sqrt(A_2/pi +r_h_2^2);

r_t_3 = r_t_2; % Tip radius at 3 [m]

r_m_2 = 0.5*(r_t_2+r_h_2);

U = N_rads * r_m_2;

r_m_3 = r_m_2; % Mean radius at 3 [m]

r_h_3 = r_h_2; % Hub radius at 3 [m]
r_h_1 = r_h_2; % Hub radius at 1 [m] Assuming hub radius is cst


r_t_1 = sqrt(A_1/pi + r_h_1^2); % Tip radius at 1 [m]
r_m_1 = 0.5*(r_t_1 + r_h_1); % Mean radius at 1 [m]

if r_m_2 < r_h_2
    error = 1;
    fprintf('Mean radius is less than hub')
else
    fprintf('Hub radius is %5.4f inches\n',r_h_2*39.37)
    fprintf('Vane inlet mean radius is %5.4f inches\n',r_m_1*39.37)
    fprintf('Vane inlet tip radius is %5.4f inches\n',r_t_1*39.37)
    fprintf('Mean blade radius is %5.4f inches\n',r_m_2*39.37)
    fprintf('blade tip radius is %5.4f inches\n',r_t_2*39.37)
    fprintf('Rotational Speed of %7f RPM\n',N_rpm)
end

% Station 2 Absolute Velocity Triangle (Part 2)
Vu_2 = Ws/U - Vu_3;
%                 if Vu_2 > V_2
%                     continue
%                 end
Va_2 = sqrt(V_2^2-Vu_2^2);
alpha_2 = atand(Vu_2/Va_2);
rho_2 = m_2/(Va_2*A_2);
P_2 = rho_2*Rg*T_2;
Po_2 = P_2 * (To_2/T_2)^(gamma/(gamma-1));

% Station 2 Relative Velocity Triangle
Vru_2 = Vu_2 - U; % Relative swirl velocity [m/s]
Vr_2 = sqrt(Va_2^2 + Vru_2^2); % Relative velocity [m/s]
alpha_r_2 = atand(Vru_2 / Va_2); % Relative swirl angle [deg]
To_r_2 = T_2 + Vr_2^2/(2*Cp);
Po_r_2 = P_2 * (T_2 / To_r_2)^(-gamma / (gamma - 1)); % Relative Total Pressure[Pa]
M_r_2 = Vr_2/sqrt(gamma*Rg*T_2); % Relative Mach

% Station 3 Relative Velocity Triangle
Vru_3 = U + Vu_3; % Relative swirl velocity [m/s]
Vr_3 = sqrt(Vru_3^2 + Va_3^2); % Relative Velocity [m/s]
alpha_r_3 = atand(Vru_3/Va_3); % Relative swirl angle [deg]
To_r_3 = T_3 + Vr_3^2/(2*Cp);
Po_r_3 = P_3 * (T_3 / To_r_3)^(-gamma / (gamma - 1)); % Relative Total Pressure[Pa]
M_r_3 = Vr_3/sqrt(gamma*Rg*T_3); % Relative Mach

% double check the velocity reaction
R_check = (Vr_3^2 - Vr_2^2)/(Vr_3^2 - Vr_2^2 + V_2^2 - V_3^2);
R_check2 = (Vr_3^2 - Vr_2^2)/(2*U*(Vu_2-Vu_3));

if Po_1 > Po_2 && Po_2 > Po_3 && Po_r_2 > Po_r_3
else
    error = 1;
    error_p=1;
    fprintf('!!!!! Pressure physics is not respected !!!!!\n')
end

%%% Incidence
beta_1 = alpha_1 - inc_1_des;
beta_r_2 = alpha_r_2 - inc_2_des;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%% HUB VELOCITY TRIANGLES %%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

U_h = r_h_2 * U/r_m_2;

% Station 1
Vu_h_1 = (r_m_1/r_h_1) * Vu_1; % hub swirl vel [m/s]
Va_h_1 = Va_1; % hub axial vel [m/s]
V_h_1 = sqrt(Va_h_1^2 + Vu_h_1^2); % hub velocity [m/s]
alpha_h_1 = atand(Vu_h_1/Va_h_1); % Swirl angle [deg]
T_h_1 = To_1 - V_h_1^2/(2*Cp);
M_h_1 = V_h_1 / sqrt(gamma*Rg*T_h_1);

% Station 2
Va_h_2 = Va_2;
Vu_h_2 = (r_m_2/r_h_2) * Vu_2; % hub swirl vel [m/s]
V_h_2 = sqrt(Va_h_2^2 + Vu_h_2^2); % hub velocity [m/s]
alpha_h_2 = atand(Vu_h_2/Va_h_2); % Swirl angle [deg]
Vru_h_2 = Vu_h_2 - U_h; % Relative swirl velocity [m/s]
Vr_h_2 = sqrt(Vru_h_2^2 + Va_h_2^2); % Relative velocity [m/s]
alpha_r_h_2 = atand(Vru_h_2/Va_h_2); % Relative swirl angle [m/s]
T_h_2 = To_2 - V_h_2^2/(2*Cp);
M_h_2 = V_h_2 / sqrt(gamma*Rg*T_h_2);
Tor_h_2 = T_h_2 + Vr_h_2^2/(2*Cp);
Mr_h_2 = Vr_h_2/sqrt(gamma*Rg*T_h_2);

% Station 3
Va_h_3 = Va_3;
Vu_h_3 = (r_m_3/r_h_3) * Vu_3; % hub swirl vel [m/s]
V_h_3 = sqrt(Va_h_3^2 + Vu_h_3^2); % hub velocity [m/s]
alpha_h_3 = atand(Vu_h_3/Va_h_3); % Swirl angle [deg]
Vru_h_3 = U_h + Vu_h_3; % Relative swirl velocity [m/s]
Vr_h_3 = sqrt(Va_h_3^2 + Vru_h_3^2); % Relative velocity [m/s]
alpha_r_h_3 = atand(Vru_h_3/Va_h_3); % Relative swirl angle [m/s]
T_h_3 = To_3 - V_h_3^2/(2*Cp);
M_h_3 = V_h_3 / sqrt(gamma*Rg*T_h_3);
Tor_h_3 = T_h_3 + Vr_h_3^2/(2*Cp);
Mr_h_3 = Vr_h_3/sqrt(gamma*Rg*T_h_3);

R_hub = (Vr_h_3^2 - Vr_h_2^2)/(Vr_h_3^2 - Vr_h_2^2 + V_h_2^2 - V_h_3^2);

if R_hub < 0
    error = 1;
    fprintf('!!!!!! Hub reaction is negative !!!!!\n')
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%% TIP VELOCITY TRIANGLES %%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

U_t = U * (r_t_2/r_m_2);

% Station 1
Vu_t_1 = (r_m_1/r_t_1) * Vu_1; % tip swirl vel [m/s]
Va_t_1 = Va_1; % tip axial vel [m/s]
V_t_1 = sqrt(Va_t_1^2 + Vu_t_1^2); % tip velocity [m/s]
alpha_t_1 = atand(Vu_t_1/Va_t_1); % Swirl angle [deg]
T_t_1 = To_1 - V_t_1^2/(2*Cp);
M_t_1 = V_t_1 / sqrt(gamma*Rg*T_t_1);

% Station 2
Va_t_2 = Va_2;
Vu_t_2 = (r_m_2/r_t_2) * Vu_2; % tip swirl vel [m/s]
V_t_2 = sqrt(Va_t_2^2 + Vu_t_2^2); % tip velocity [m/s]
alpha_t_2 = atand(Vu_t_2/Va_t_2); % Swirl angle [deg]
Vru_t_2 = Vu_t_2 - U_t; % Relative swirl velocity [m/s]
Vr_t_2 = sqrt(Vru_t_2^2 + Va_t_2^2); % Relative velocity [m/s]
alpha_r_t_2 = atand(Vru_t_2/Va_t_2); % Relative swirl angle [m/s]
T_t_2 = To_2 - V_t_2^2/(2*Cp);
M_t_2 = V_t_2 / sqrt(gamma*Rg*T_t_2);
Tor_t_2 = T_t_2 + Vr_t_2^2/(2*Cp);
Mr_t_2 = Vr_t_2/sqrt(gamma*Rg*T_t_2);
Po_r_2_t = P_2 * (T_2 / Tor_t_2)^(-gamma / (gamma - 1)); % Relative Total Pressure[Pa]


% Station 3
Va_t_3 = Va_3;
Vu_t_3 = (r_m_3/r_t_3) * Vu_3; % hub swirl vel [m/s]
V_t_3 = sqrt(Va_t_3^2 + Vu_t_3^2); % hub velocity [m/s]
alpha_t_3 = atand(Vu_t_3/Va_t_3); % Swirl angle [deg]
Vru_t_3 = U_t + Vu_t_3; % Relative swirl velocity [m/s]
Vr_t_3 = sqrt(Va_t_3^2 + Vru_t_3^2); % Relative velocity [m/s]
alpha_r_t_3 = atand(Vru_t_3/Va_t_3); % Relative swirl angle [m/s]
T_t_3 = To_3 - V_t_3^2/(2*Cp);
M_t_3 = V_t_3 / sqrt(gamma*Rg*T_t_3);
Tor_t_3 = T_h_3 + Vr_t_3^2/(2*Cp);
Mr_t_3 = Vr_t_3/sqrt(gamma*Rg*T_t_3);
Po_r_3_t = P_3 * (T_3 / Tor_t_3)^(-gamma / (gamma - 1)); % Relative Total Pressure[Pa]

R_tip = (Vr_t_3^2 - Vr_t_2^2)/(Vr_t_3^2 - Vr_t_2^2 + V_t_2^2 - V_t_3^2); % Reaction at tip

if Po_r_2_t > Po_r_3_t
else
    error = 1;
    error_p=1;
    fprintf('!!!!! Pressure physics is not respected !!!!!\n')
end

%%%%%%%%%%%%%%%%
%%% GEOMETRY %%%
%%%%%%%%%%%%%%%%

vane_height = 0.5*((r_t_1-r_h_1)+(r_t_2-r_h_2)); % Vane height [m]
blade_height = r_t_3-r_h_3; % Blade height [m]

vane_actual_chord = vane_height/AR_v; % Vane actual chord [m]
blade_actual_chord = blade_height/AR_b; % Blade actual chord [m]

%%%%%%%%%%%%%%%%%%%%%
%%% STAGGER ANGLE %%%
%%%%%%%%%%%%%%%%%%%%%

beta_1_d = stagger_data(:,1); % beta before vane or blade
beta_2_d = stagger_data(:,2:end); % beta after vane or blade
beta_2_vector = [80,75,70,65,60,55,50]; % set of curves that define the plot

stagger_vane = interp2(beta_2_vector,beta_1_d,beta_2_d,alpha_2,beta_1); % Vane stagger angle [deg]
stagger_blade = interp2(beta_2_vector,beta_1_d,beta_2_d,alpha_r_3,beta_r_2); % Blade stagger angle [deg]

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%% NUMBER OF BLADES/VANES %%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

vane_axial_chord = vane_actual_chord * cosd(stagger_vane); % Axial chord of the vane [m]
blade_axial_chord = blade_actual_chord * cosd(stagger_blade); % Axial chord of the blad [m]

vane_pitch = (zweif_vane*vane_axial_chord*0.5)/((tand(alpha_1) + tand(alpha_2))* ((cosd(alpha_2))^2)); % Vane pitch [m]
vane_number = pi*(r_m_2+r_m_1)/vane_pitch; % Number of vanes

blade_pitch = (zweif_blade*blade_axial_chord*0.5)/((tand(alpha_r_2) + tand(alpha_r_3))* ((cosd(alpha_r_3))^2)); % Blade pitch [m]
blade_number = 2*pi*r_m_2/blade_pitch; % Number of blades

fprintf('%3.2f vanes \n',vane_number);
fprintf('%3.2f blades \n',blade_number);

vane_num_int = round(vane_number,0);
blade_num_int = round(blade_number,0);

if abs(vane_num_int-vane_number) > 0.1 || abs(blade_num_int-blade_number) > 0.1
%     error = 1;
%     fprintf('Non-integer number of blades or vanes\n');
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%% THICKNESS OVER CHORD %%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%

vane_max_thickness = vane_actual_chord*interp1(tmaxc_data(:,1),tmaxc_data(:,2),beta_1 + alpha_2);
blade_max_thickness = blade_actual_chord*interp1(tmaxc_data(:,1),tmaxc_data(:,2),beta_r_2 + alpha_r_3);

%%%%%%%%%%%%%%
%%% LOSSES %%%
%%%%%%%%%%%%%%

% k_Accel (component of profile losses)

if M_2 <= 0.2
    K1_vane = 1;
else
    K1_vane = 1-1.25*(M_2-0.2);
end
K2_vane = (M_1/M_2)^2;
k_accel_vane = 1-K2_vane *(1- K1_vane);

if M_r_3<=0.2
    K1_blade = 1;
else
    K1_blade = 1-1.25*(M_r_3-0.2);
end
K2_blade = (M_r_2/M_r_3)^2;
k_accel_blade = 1-K2_blade*(1-K1_blade);


% Ksh, will always be equal to zero for vane though since M1 is 0.125
% Shock Losses
if M_h_1 <= 0.4
    P_Q_sh_vane=0;
else
    P_Q_sh_vane = (r_h_1/(0.5*(r_t_1+r_t_2)))*(0.75*(M_h_1-0.4)^1.75);
end
Ksh_vane = P_Q_sh_vane*(P_1/P_2)*(1-(1+0.5*(gamma-1)*M_1^2)^(gamma/(gamma-1)))/(1-(1+0.5*(gamma-1)*M_2^2)^(gamma/(gamma-1)));

if Mr_h_2 <= 0.4
    P_Q_sh_blade=0;
else
    P_Q_sh_blade = (r_h_3/r_t_3)*(0.75*(Mr_h_2-0.4)^1.75);
end
Ksh_blade = P_Q_sh_blade*(P_2/P_3)*(1-(1+0.5*(gamma-1)*M_r_2^2)^(gamma/(gamma-1)))/(1-(1+0.5*(gamma-1)*M_r_3^2)^(gamma/(gamma-1)));


% Yp beta = 0
S_C_vane = vane_pitch / vane_actual_chord; % Span to chord ratio
S_C_blade = blade_pitch / blade_actual_chord;


S_over_C = Yp_beta_0_data(:,1);
alpha_two = [40, 50, 60, 65, 70, 75, 80];

Yp_beta_0_vane = interp2(alpha_two,S_over_C, Yp_beta_0_data(:,2:end),alpha_2,S_C_vane); % Yp beta=0 for vane
Yp_beta_0_blade = interp2(alpha_two , S_over_C , Yp_beta_0_data(:,2:end) , alpha_r_3 , S_C_blade); % Yp beta=0 for blade


% Yp beta = alpha

S_over_C_2 = Yp_beta_alpha_data(:,1);
alpha_two = [40 , 50 , 55 , 60 , 65 , 70];

if alpha_2 > 70
    Yp_beta_alpha_vane = interp1(40:0.01:70,interp2(alpha_two,S_over_C_2, Yp_beta_alpha_data(:,2:end),40:0.01:70,S_C_vane),alpha_2,'linear','extrap');
elseif alpha_2 > 80
    error=1;
    fprintf('!!!!!! alpha 2 is over 80 !!!!!!\n')
else
    Yp_beta_alpha_vane = interp2(alpha_two,S_over_C_2, Yp_beta_alpha_data(:,2:end),alpha_2,S_C_vane); % Yp beta=alpha for vane
end
if alpha_r_3 > 70
    Yp_beta_alpha_blade = interp1(40:0.01:70,interp2(alpha_two,S_over_C_2, Yp_beta_alpha_data(:,2:end),40:0.01:70,S_C_blade),alpha_r_3,'linear','extrap');
elseif alpha_r_3>80
    error=1;
    fprintf('!!!!!! alpha relative 3 is over 80 !!!!!!\n')
else
    Yp_beta_alpha_blade = interp2(alpha_two , S_over_C_2 , Yp_beta_alpha_data(:,2:end) , alpha_r_3 , S_C_blade); % Yp beta=alpha for blade
end
 
% Yp and Kp

Yp_AMDC_vane = Yp_beta_0_vane + abs((beta_1)/alpha_2) * ((beta_1)/alpha_2) * (Yp_beta_alpha_vane - Yp_beta_0_vane) * ((vane_max_thickness/vane_actual_chord)/0.2)^((beta_1)/alpha_2);
Yp_AMDC_blade = Yp_beta_0_blade + abs((beta_r_2)/alpha_r_3) * ((beta_r_2)/alpha_r_3) * (Yp_beta_alpha_blade - Yp_beta_0_blade) * ((blade_max_thickness/blade_actual_chord)/0.2)^((beta_r_2)/alpha_r_3);

Kp_vane = 0.914 * (2/3 * Yp_AMDC_vane * k_accel_vane + Ksh_vane);
Kp_blade = 0.914 * (2/3 * Yp_AMDC_blade * k_accel_blade + Ksh_blade);

% Secondary Losses
f_AR_vane = (1-0.25 * sqrt(2 - AR_v))/AR_v;
f_AR_blade = (1-0.25 * sqrt(2 - AR_b))/AR_b;

alpha_m_vane = atand(0.5 * (tand(alpha_1) - tand(alpha_2)));
alpha_m_blade = atand(0.5 * (tand(alpha_r_2) - tand(alpha_r_3)));

Cl_vane = S_C_vane * 2 * (tand(alpha_1)+tand(alpha_2)) * cosd(alpha_m_vane);
Cl_blade = S_C_blade * 2 * (tand(alpha_r_2)+tand(alpha_r_3)) * cosd(alpha_m_blade);

Ys_AMDC_vane = 0.0334*f_AR_vane * (cosd(alpha_2)/cosd(beta_1)) * (Cl_vane/S_C_vane)^2 * (cosd(alpha_2))^2/(cosd(alpha_m_vane))^3;
Ys_AMDC_blade = 0.0334*f_AR_blade * (cosd(alpha_r_3)/cosd(beta_r_2)) * (Cl_blade/S_C_blade)^2 * (cosd(alpha_r_3))^2/(cosd(alpha_m_blade))^3;

K3_vane = (1/(vane_height/vane_axial_chord))^2;
K3_blade = (1/(blade_height/blade_axial_chord))^2;

Ks_vane = 1-K3_vane * (1-k_accel_vane);
Ks_blade = 1-K3_blade * (1-k_accel_blade);

Ys_vane = 1.2 * Ys_AMDC_vane*Ks_vane; % Secondary flow coefficient
Ys_blade = 1.2 * Ys_AMDC_blade*Ks_blade;

% Trailing Edge Losses

throat_vane = vane_pitch*cosd(alpha_2); % Throat opening length
throat_blade = blade_pitch*cosd(alpha_r_3);

throat_vane_h = vane_pitch*cosd(alpha_h_2);
throat_blade_h = blade_pitch*cosd(alpha_r_h_2);

throat_vane_t = vane_pitch*cosd(alpha_t_2);
throat_blade_t = blade_pitch*cosd(alpha_r_t_2);

thick_open_vane = TE_v / throat_vane; % TE thickness to throat opening
thick_open_blade = TE_b / throat_blade;

tet_b_0_vane = interp1(delta_phi_data(:,1),delta_phi_data(:,2),thick_open_vane); % delta phi^2 for beta = 0
tet_b_a_vane = interp1(delta_phi_data(:,1),delta_phi_data(:,3),thick_open_vane);

tet_b_0_blade = interp1(delta_phi_data(:,1),delta_phi_data(:,2),thick_open_blade); % delta phi^2 for beta = alpha
tet_b_a_blade = interp1(delta_phi_data(:,1),delta_phi_data(:,3),thick_open_blade);

tet_vane = tet_b_0_vane + abs((beta_1)/alpha_2)*((beta_1)/alpha_2)*( tet_b_a_vane- tet_b_0_vane); % delta phi^2 net
tet_blade = tet_b_0_blade + abs((beta_r_2)/alpha_r_3)*((beta_r_2)/alpha_r_3)*( tet_b_a_blade - tet_b_0_blade);

KTE_vane = ((1 - 0.5 * (gamma-1) * M_2^2*((1/(1-tet_vane))-1))^(-gamma/(gamma-1)) - 1)/(1-(1+0.5 * (gamma-1)*M_2^2)^(-gamma/(gamma-1))); % trailing edge loss coefficient
KTE_blade = ((1 - 0.5 * (gamma-1) * M_r_3^2*((1/(1-tet_blade))-1))^(-gamma/(gamma-1)) - 1)/(1-(1+0.5 * (gamma-1)*M_r_3^2)^(-gamma/(gamma-1)));

% Reynolds Number Calculations

Re_vane = (rho_2*V_2*vane_actual_chord)/ interp1(kin_visc(:,1), kin_visc(:,2), T_2); % Reynolds Number
Re_blade = (rho_3*Vr_3*blade_actual_chord)/ interp1(kin_visc(:,1), kin_visc(:,2), T_3);

fre_vane = 0; % initialization
fre_blade = 0;

if Re_vane < 2e5
    fre_vane = (Re_vane/2e5)^-.4;
elseif Re_vane > 10e6
    fre_vane = (Re_vane/10e6)^-.2;
else
    fre_vane = 1;
end

if Re_blade < 2e5
    fre_blade = (Re_blade/2e5)^-.4;
elseif Re_blade > 10^6
    fre_blade = (Re_blade/10^6)^-.2;
else
    fre_blade = 1;
end

% Kt (overall loss before tip losses)

Kt_vane = fre_vane * Kp_vane + Ys_vane  + KTE_vane;
Kt_blade = fre_blade * Kp_blade + Ys_blade + KTE_blade;

% Tip clearance losses and total-to-total efficiency

zeta_vane = Kt_vane/(1+0.5*gamma*M_2^2);
zeta_blade = Kt_blade/(1+0.5*gamma*M_r_3^2);

eta_o = (1 + ((zeta_vane*V_2^2 + zeta_blade*Vr_3^2)/(2*Cp*(To_1-To_3))))^(-1);

delta_k = blade_tip_clearance*blade_height + 0.005*0.0254; 

de = 0.93 * (eta_o) * (r_t_3/r_m_3) * (delta_k/(blade_height*cosd(alpha_r_3)));
de_init = de;
eta_des = 1;
K_clr = 0;

% Iterate on the delta eta by adding K_clr and converging on delta eta
while abs(de - abs(eta_des - eta_o)) > 0.005
    K_clr = K_clr + 0.0001;
    Kt_blade = Kt_blade + 0.0001;
    zeta_blade = Kt_blade/(1+0.5*gamma*M_r_3^2);
    eta_des = (1 + ((zeta_vane*V_2^2 + zeta_blade*Vr_3^2)/(2*Cp*(To_1-To_3))))^(-1);
    de = 0.93 * (eta_des) * (r_t_3/r_m_3) * (delta_k/(blade_height*cosd(alpha_r_3)));
end

fprintf('Total-to-total efficiency of %4.3f \n',eta_des)
% fprintf('Desired efficiency is %4.3f \n', eta_i)

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%% OFF DESIGN PERFORMANCE %%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

U_od = 0.9 * U; % off design U speed
To_3_od = To_2 - .9*Ws/Cp;

% Fix Absolute Velocity Triangle

T_3_od = To_3_od / (1 + 0.5 * (gamma - 1) * M_3^2); % Temperature at Blade Exit [K]
V_3_od = M_3 * sqrt(gamma * Rg * T_3_od); % Blade Exit Inlet Velocity [m/s]
Va_3_od = sqrt(V_3_od^2 - Vu_3^2);
alpha_3_od = atand(Vu_3/Va_3_od);
rho_3_od = m_3/(A_3*Va_3_od);
P_3_od = rho_3_od*Rg*T_3_od;

% Station 2 Relative Velocity Triangle
Vru_2_od = Vu_2 - U_od; % Relative swirl velocity [m/s]
Vr_2_od = sqrt(Va_2^2 + Vru_2_od^2); % Relative velocity [m/s]
alpha_r_2_od = atand(Vru_2_od / Va_2); % Relative swirl angle [deg]
To_r_2_od = T_2 + Vr_2_od^2/(2*Cp);
Po_r_2_od = P_2 * (T_2 / To_r_2_od)^(-gamma / (gamma - 1)); % Relative Total Pressure[Pa]
M_r_2_od = Vr_2_od/sqrt(gamma*Rg*T_2); % Relative Mach

% Station 3 Relative Velocity Triangle
Vru_3_od = U_od + Vu_3; % Relative swirl velocity [m/s]
Vr_3_od = sqrt(Vru_3_od^2 + Va_3_od^2); % Relative Velocity [m/s]
alpha_r_3_od = atand(Vr_3_od/Va_3_od); % Relative swirl angle [deg]
To_r_3_od = T_3_od + Vr_3_od^2/(2*Cp);
Po_r_3_od = P_3_od * (T_3_od / To_r_3_od)^(-gamma / (gamma - 1)); % Relative Total Pressure[Pa]
M_r_3_od = Vr_3_od/sqrt(gamma*Rg*T_3); % Relative Mach

% Station 2 hub velocity triangle, for Ksh losses
Vru_h_2_od = Vu_h_2 -0.9*U_h; % Relative swirl velocity [m/s]
Vr_h_2_od = sqrt(Vru_h_2_od^2 + Va_h_2^2); % Relative velocity [m/s]
alpha_r_h_2_od = atand(Vru_h_2_od/Va_h_2); % Relative swirl angle [m/s]
Tor_h_2_od = T_h_2 + Vr_h_2_od^2/(2*Cp);
Mr_h_2_od = Vr_h_2_od/sqrt(gamma*Rg*T_h_2);

inc_2_od = alpha_r_2_od - beta_r_2;

%%% losses for off-design blade

% k_Accel (component of profile losses)

if M_r_3_od<=0.2
    K1_blade_od = 1;
else
    K1_blade_od = 1-1.25*(M_r_3_od-0.2);
end
K2_blade_od = (M_r_2_od/M_r_3_od)^2;

% Yp and Kp - Moustapha method
LE_dia = 0.00254; % Leading edge diameter [m] (assumed for now)

Chi_p = (LE_dia/blade_pitch)^-1.6 * (cosd(beta_r_2) / cosd(alpha_r_3))^-2 * (alpha_r_2_od - alpha_r_2);

if Chi_p > 0 && Chi_p < 800
    phi_p_od = 0.788e-5 * Chi_p + 0.56e-7 * Chi_p^2 + 0.4e-10 * Chi_p^3 + 2.054e-19 * Chi_p^6;
elseif Chi_p > -800 && Chi_p < 0
    phi_p_od = -5.1734e-6 * Chi_p + 7.6902e-9 * Chi_p^2;
else
    error = 1;
    fprintf('Chi_p out of bounds \n');
    phi_p_od = 9999;
end

phi_p0 = (1 + Kp_blade/(K1_blade + Kp_blade*K2_blade))^-1;

phi_p = phi_p0+phi_p_od;

Kp_blade_od = (K1_blade_od * (1-phi_p))/(phi_p-K2_blade_od*(1-phi_p));

% Secondary Losses - Moustapha method
Chi_s = ((alpha_r_2_od - beta_r_2)/(beta_r_2 + alpha_r_3)) * (cosd(beta_r_2) / cosd(alpha_r_3))^-1.5 * (LE_dia/blade_actual_chord)^-0.3;

if Chi_s > -0.4 && Chi_s < 0
    Ys_blade_od = Ys_blade * (exp(0.9*Chi_s) + 13*Chi_s^2 + 400*Chi_s^4);
elseif Chi_s >0 && Chi_s <0.3
    Ys_blade_od = Ys_blade * exp(0.9*Chi_s);
else
    error = 1;
    fprintf('Chi_s out of bounds \n');
    Ys_blade_od = 9999;
end

% Trailing Edge Losses

tet_blade_od = tet_b_0_blade + abs((beta_r_2)/alpha_r_3_od)*((beta_r_2)/alpha_r_3_od)*( tet_b_a_blade - tet_b_0_blade);

KTE_blade_od = ((1 - 0.5 * (gamma-1) * M_r_3_od^2*((1/(1-tet_blade_od))-1))^(-gamma/(gamma-1)) - 1)/(1-(1+0.5 * (gamma-1)*M_r_3_od^2)^(-gamma/(gamma-1)));

% Reynolds Number Calculations

Re_blade_od = (rho_3*Vr_3_od*blade_actual_chord)/ interp1(kin_visc(:,1), kin_visc(:,2), T_3_od);

fre_blade_od = 0;

if Re_blade_od < 2e5
    fre_blade_od = (Re_blade_od/2e5)^-.4;
elseif Re_blade_od > 10^6
    fre_blade_od = (Re_blade_od/10^6)^-.2;
else
    fre_blade_od = 1;
end

% Kt (overall loss before tip losses)

Kt_blade_od = fre_blade_od * Kp_blade_od + Ys_blade_od + KTE_blade_od;

% Tip clearance losses and total-to-total efficiency

zeta_blade_od = Kt_blade_od/(1+0.5*gamma*M_r_3_od^2);

eta_o_od = (1 + ((zeta_vane*V_2^2 + zeta_blade_od*Vr_3_od^2)/(2*Cp*(To_1-To_3_od))))^(-1);

de2 = 0.93 * (eta_o_od) * (r_t_3/r_m_3) * (delta_k/(blade_height*cosd(alpha_r_3_od)));

de_init_od = de2;
eta_des_od = 1;
K_clr_od = 0;

% Iterate on the delta eta by adding K_clr and converging on delta eta
while abs(de2 - abs(eta_des_od - eta_o_od)) > 0.005
    K_clr_od = K_clr_od + 0.0001;
    Kt_blade_od = Kt_blade_od + 0.0001;
    zeta_blade_od = Kt_blade_od/(1+0.5*gamma*M_r_3_od^2);
    eta_des_od = (1 + ((zeta_vane*V_2^2 + zeta_blade_od*Vr_3_od^2)/(2*Cp*(To_1-T_3_od))))^(-1);
    de2 = 0.93 * (eta_des_od) * (r_t_3/r_m_3) * (delta_k/(blade_height*cosd(alpha_r_3_od)));
end

fprintf('Total-to-total off-design efficiency of %4.3f \n',eta_des_od)
%                 if error_p == 0 && error == 0
%                     pause
%                 end
%                 end
%             end
%         end
%     end
% end
% 
% end