clc;clear all;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% Author: Noah Sadaka
% Date: February 17 2019
% Course: Gas Turbine Design

% Purpose: Create a model for the Turbine Design

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% Notes: Calculate the absolute triangles and determine U from there

% Stage Nomenclature
% Station 1 = Vane Inlet
% Station 2 = Vane Exit and Blade Inlet
% Station 3 = Blade Exit

% Cycle Analysis Variables

eta_i = 0.85; % desired efficiency

error = 0; % flag for errors

To_1 = 1147.98; % Turbine Inlet Temperature [K]
Po_1 = 315940.86; % Turbine Inlet Pressure [Pa]
m_1 = 4.842; % Turbine Inlet Massflow [kg/s]

To_2 = 1126.67; % Total temperature after bleed air addition [K]
%To_2 = To_1;

m_2 = 5.00; % Mass flow through blade [kg/s]

To_3 = 972.74; % Total temperature after station [K]
Po_3 = 149453.77; % Total pressure exiting blade, before ITD [Pa]
m_3 = m_2; % Mass flow exiting blade [kg/s]

W = 883830.57; % Stage work [W]

% Gas Constants and Coefficients
Cp = 1148; % Gas specific heat capacity [kg/JK]
gamma = 1.333; % Gas Gamma
Rg = 287; % Gas Constant

% RAYMOND: values to put in the input data csv table:
% alpha_3, M_3, R, AN2_max, the inc and dev, zweif blade and vane
% Thanks!

% Assumptions
alpha_3 = 0; % Blade exit swirl angle [deg] Range: -5 to 30
M_3 = 0.417; % Blade exit mach number. Range: 0.3-0.45
R = 0.54; % Reaction at the meanline.
max_U_h = 1100*0.3048; % Blade speed at meanline [m/s]
AN2 = 4.5E10; % AN2 [in2 rpm 0.5]
inc_1 = 0; % Vane Incidence [deg]
dev_2 = 0; % Vane Deviation [deg]
inc_2 = 0; % Blade Incidence [deg]vane_axial_chord
dev_3 = 0; % Blade Deviation [deg]
zweif_vane = 0.8;
zweif_blade = 0.9;
blade_tip_clearance = 0.009;

% Given Vane Variables
AR_v = 0.7; % Vane Aspect Ratio
TE_v = 0.045*0.0254; % Vane TE Thickness [m]

% Given Blade Variables
AR_b = 1.45; % Blade Aspect Ratio
TE_b = 0.025*0.0254; % Blade TE Thickness [m]

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
Vu_1 = sqrt(V_1^2 - Va_1^2); % Vane Inlet Tangential Velocity
A_1 = m_1/(rho_1 * Va_1); % Annulus Area at station 1 [m^2]

% Station 3 Absolute Velocity Triangle
T_3 = To_3 / (1 + 0.5 * (gamma - 1) * M_3^2); % Temperature at Blade Exit [K]
P_3 = Po_3 * (T_3 / To_3)^(gamma / (gamma - 1)); % Pressure at Blade Exit [Pa]
rho_3 = P_3 / (Rg * T_3); % Density at Blade Exit [kg/m3]
V_3 = M_3 * sqrt(gamma * Rg * T_3); % Blade Exit Inlet Velocity [m/s]
Va_3 = V_3 * cosd(alpha_3); % Axial velocity at station 3 [m/s]
A_3 = m_3 / (rho_3 * Va_3); % Annulus Area at station 3 [m^2]
Vu_3 =  sqrt(V_3^2 - Va_3^2); % Swirl velocity [m/s]
alpha_3 = atand(Vu_3/Va_3); % Swirl angle [deg]

% Station 2 Absolute Velocity Triangle
Ws = W/m_2; % Specific work [W/kg]
A_2 = A_3; % Assumption: rotor area is constant
T_2 = R * (T_1 - T_3) + T_3; % Temperature at station 2 [K]
V_2 = sqrt(2*Cp*(To_2 - T_2)); % Absolute Velocity [m/s]
M_2 = V_2/sqrt(T_2 * gamma * Rg); % Mach
Po_2 = Po_1 * (To_2/To_1)^(gamma/(gamma-1)); % Total pressure [Pa]
P_2 = Po_2 * (1+0.5 * (gamma-1)*M_2^2)^(-gamma/(gamma-1)); % Static pressure [Pa]
rho_2 = P_2/(Rg * T_2); % Density [kg/m3]
Va_2 = m_2 / (rho_2 * A_2); % Axial Velocity [m/s]
Vu_2 = sqrt(V_2^2 - Va_2^2); % Swirl Velocity [m/s]
alpha_2 = atand(Vu_2 / Va_2); % Swirl angle [deg]
U = Ws / (Vu_2 + Vu_3); % Blade Velocity [m/s]

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
alpha_r_3 = atand(Vr_3/Va_3); % Relative swirl angle [deg]
To_r_3 = T_3 + Vr_3^2/(2*Cp);
Po_r_3 = P_3 * (T_3 / To_r_3)^(-gamma / (gamma - 1)); % Relative Total Pressure[Pa]
M_r_3 = Vr_3/sqrt(gamma*Rg*T_3); % Relative Mach

% double check the velocity reaction
R_check = (Vr_3^2 - Vr_2^2)/(Vr_3^2 - Vr_2^2 + V_2^2 - V_3^2);
R_check2 = (Vr_3^2 - Vr_2^2)/(2*U*(Vu_2-Vu_3));

% if R_check ~= R
%     error = 1;
%     fprintf('Reactions do not match fuuuuuck\n')
% end

if Po_1 > Po_2 && Po_r_2 > Po_r_3
else
    error = 1;
    fprintf('Pressure physics is not respected Tabarnak\n')
end

% if T_3 < T_2 && T_2 > T_1 && T_r_2 > T_r_3
%     error = 1;
%     fprintf('Temperature physics is not respected fuck\n');
% end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%% RADII AND STRUCTURAL LIMITATIONS %%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

A_rpm = 1550*A_2;
N_rpm = sqrt(AN2/A_rpm); % Rotation speed [rpm]
N_rpm_max = sqrt(AN2/A_rpm);
%N_rpm = 15000;
N_rads = N_rpm * (1/60) * (2*pi); % Rotation speed [rad/s]
N_rads_max = N_rpm_max * (1/60) * (2*pi); % Rotation speed [rad/s]

r_h_min = max_U_h/N_rads;

r_m_2 = U/N_rads; % Mean radius at 2 [m]

r_h_2 = (-A_2/pi + 4*r_m_2^2)/(4 * r_m_2); % Hub radius at 2 [m]

% Iterate to find an r_h that works

interval = 0.0001:0.0001:0.5;
check_rh = 0;

if r_h_2 > r_h_min % if rhub is violating the max rim speed
    fprintf('Using RPM does not work. Testing with max hub speed\n');
    for i=1:length(interval)
        if abs(interval(i)-((max_U_h/U)*0.5*(interval(i)+ sqrt(A_2/pi + interval(i)^2)))) < 0.0001
            r_h_2 = interval(i);
            r_m_2 = r_h_2*U/max_U_h;
            check_rh = 1;
        end
    end
end

if check_rh == 1
    fprintf('Successful r_h calculated using max hub speed\n')
    fprintf('rotational speed is %6f\n',(U/r_m_2)*(30/pi))
end

if (U/r_m_2)*(30/pi) > N_rpm_max
    error = 1;
end

r_m_3 = r_m_2; % Mean radius at 3 [m]

r_h_3 = r_h_2; % Hub radius at 3 [m]
r_h_1 = r_h_2; % Hub radius at 1 [m] Assuming hub radius is cst

r_t_2 = 2*r_m_2 - r_h_2; % Tip radius at 2 [m]
r_t_3 = r_t_2; % Tip radius at 3 [m]

r_t_1 = sqrt(A_1/pi + r_h_1^2); % Tip radius at 1 [m]
r_m_1 = 0.5*(r_t_1 + r_h_1); % Mean radius at 1 [m]

if r_m_2 < r_h_2
    error = 1;
    fprintf('Mean radius is less than hub, fuck!')
else
    fprintf('Hub radius is %5.4f inches\n',r_h_2*39.37)
    fprintf('Vane inlet mean radius is %5.4f inches\n',r_m_1*39.37)
    fprintf('Vane inlet tip radius is %5.4f inches\n',r_t_1*39.37)
    fprintf('Mean blade radius is %5.4f inches\n',r_m_2*39.37)
    fprintf('blade tip radius is %5.4f inches\n',r_t_2*39.37)
end


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
    fprintf('Hub reaction is negative\n')
end

if U_h > max_U_h
    error = 1;
    fprintf('fuck velocities\n')
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

R_tip = (Vr_t_3^2 - Vr_t_2^2)/(Vr_t_3^2 - Vr_t_2^2 + V_t_2^2 - V_t_3^2); % Reaction at tip

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

stagger_data = csvread('stagger_angle.csv',1); % Stagger angle data
beta_1 = stagger_data(:,1); % beta before vane or blade
beta_2 = stagger_data(:,2:end); % beta after vane or blade
beta_2_vector = [80,75,70,65,60,55,50]; % set of curves that define the plot

stagger_vane = interp2(beta_2_vector,beta_1,beta_2,alpha_2+dev_2,alpha_1+inc_1); % Vane stagger angle [deg]
stagger_blade = interp2(beta_2_vector,beta_1,beta_2,alpha_r_3+dev_3,alpha_r_2+inc_2); % Blade stagger angle [deg]

%%%%%%%%%%%%%%%%%%%%%%%%
%%% NUMBER OF MORTYS %%%
%%%%%%%%%%%%%%%%%%%%%%%%

vane_axial_chord = vane_actual_chord * cosd(stagger_vane); % Axial chord of the vane [m]
blade_axial_chord = blade_actual_chord * cosd(stagger_blade); % Axial chord of the blad [m]

vane_pitch = (zweif_vane*vane_axial_chord*0.5)/((tand(alpha_1) + tand(alpha_2))* ((cosd(alpha_2))^2)); % Vane pitch [m]
vane_number = pi*(r_m_2+r_m_1)/vane_pitch; % Number of vanes

blade_pitch = (zweif_blade*blade_axial_chord*0.5)/((tand(alpha_r_2) + tand(alpha_r_3))* ((cosd(alpha_r_3))^2)); % Blade pitch [m]
blade_number = 2*pi*r_m_2/blade_pitch; % Number of blades

fprintf('%3.2f vanes \n',vane_number);
fprintf('%3.2f blades \n',blade_number);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%% THICKNESS OVER CHORD %%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%

tmaxc_data = csvread('tmaxc_v_betas.csv',1);

vane_max_thickness = vane_actual_chord*interp1(tmaxc_data(:,1),tmaxc_data(:,2),alpha_1 + inc_1 + alpha_2 + dev_2);
blade_max_thickness = blade_actual_chord*interp1(tmaxc_data(:,1),tmaxc_data(:,2),alpha_r_2 + inc_2 + alpha_r_3 + dev_3);

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

Yp_beta_0_data = csvread('Yp_beta1_0.csv',1);
S_over_C = Yp_beta_0_data(:,1);
alpha_two = [40,50,60,65,70,75,80];

Yp_beta_0_vane = interp2(alpha_two,S_over_C, Yp_beta_0_data(:,2:end),alpha_2,S_C_vane); % Yp beta=0 for vane
Yp_beta_0_blade = interp2(alpha_two , S_over_C , Yp_beta_0_data(:,2:end) , alpha_r_3 , S_C_blade); % Yp beta=0 for blade

clear alpha_two S_over_C Yp_beta_0_data

% Yp beta = alpha
Yp_beta_alpha_data = csvread('Yp_beta1_alfa1.csv',1);
S_over_C = Yp_beta_alpha_data(:,1);
alpha_two = [40 , 50 , 55 , 60 , 65 , 70];

Yp_beta_alpha_vane = interp2(alpha_two,S_over_C, Yp_beta_alpha_data(:,2:end),alpha_2,S_C_vane); % Yp beta=alpha for vane
Yp_beta_alpha_blade = interp2(alpha_two , S_over_C , Yp_beta_alpha_data(:,2:end) , alpha_r_3 , S_C_blade); % Yp beta=alpha for blade


% Yp and Kp

Yp_AMDC_vane = Yp_beta_0_vane + abs((alpha_1+inc_1)/alpha_2) * ((alpha_1+inc_1)/alpha_2) * (Yp_beta_alpha_vane - Yp_beta_0_vane) * ((vane_max_thickness/vane_actual_chord)/0.2)^((alpha_1+inc_1)/alpha_2);
Yp_AMDC_blade = Yp_beta_0_blade + abs((alpha_r_2+inc_2)/alpha_r_3) * ((alpha_r_2+inc_2)/alpha_r_3) * (Yp_beta_alpha_blade - Yp_beta_0_blade) * ((blade_max_thickness/blade_actual_chord)/0.2)^((alpha_r_2+inc_2)/alpha_r_3);

Kp_vane = 0.914 * (2/3 * Yp_AMDC_vane * k_accel_vane + Ksh_vane);
Kp_blade = 0.914 * (2/3 * Yp_AMDC_blade * k_accel_blade + Ksh_blade);


% Secondary Losses
f_AR_vane = (1-0.25 * sqrt(2 - AR_v))/AR_v;
f_AR_blade = (1-0.25 * sqrt(2 - AR_b))/AR_b;

alpha_m_vane = atand(0.5 * (tand(alpha_1) - tand(alpha_2)));
alpha_m_blade = atand(0.5 * (tand(alpha_r_2) - tand(alpha_r_3)));

Cl_vane = S_C_vane * 2 * (tand(alpha_1)+tand(alpha_2)) * cosd(alpha_m_vane);
Cl_blade = S_C_blade * 2 * (tand(alpha_r_2)+tand(alpha_r_3)) * cosd(alpha_m_blade);

Ys_AMDC_vane = 0.0334*f_AR_vane * (cosd(alpha_2)/cosd(alpha_1 + inc_1)) * (Cl_vane/S_C_vane)^2 * (cosd(alpha_2))^2/(cosd(alpha_m_vane))^2;
Ys_AMDC_blade = 0.0334*f_AR_blade * (cosd(alpha_r_3)/cosd(alpha_r_2 + inc_2)) * (Cl_blade/S_C_blade)^2 * (cosd(alpha_r_3))^2/(cosd(alpha_m_blade))^2;

K3_vane = (1/(vane_height/vane_axial_chord))^2;
K3_blade = (1/(blade_height/blade_axial_chord))^2;

Ks_vane = 1-K3_vane * (1-k_accel_vane);
Ks_blade = 1-K3_blade * (1-k_accel_blade);

Ys_vane = 1.2 * Ys_AMDC_vane*Ks_vane; % Secondary flow coefficient
Ys_blade = 1.2 * Ys_AMDC_blade*Ks_blade;

% Trailing Edge Losses

delta_phi_data = csvread('fig_14_enecoef.csv',1);

throat_vane = vane_pitch*cosd(alpha_2 + dev_2); % Throat opening length
throat_blade = blade_pitch*cosd(alpha_3+dev_3);

thick_open_vane = TE_v / throat_vane; % TE thickness to throat opening
thick_open_blade = TE_b / throat_blade;

tet_b_0_vane = interp1(delta_phi_data(:,1),delta_phi_data(:,2),thick_open_vane); % delta phi^2 for beta = 0
tet_b_a_vane = interp1(delta_phi_data(:,1),delta_phi_data(:,3),thick_open_vane);

tet_b_0_blade = interp1(delta_phi_data(:,1),delta_phi_data(:,2),thick_open_blade); % delta phi^2 for beta = alpha
tet_b_a_blade = interp1(delta_phi_data(:,1),delta_phi_data(:,3),thick_open_blade);

tet_vane = tet_b_0_vane + abs((alpha_1+inc_1)/alpha_2)*((alpha_1+inc_1)/alpha_2)*( tet_b_a_vane- tet_b_0_vane); % delta phi^2 net
tet_blade = tet_b_0_blade + abs((alpha_r_2+inc_2)/alpha_r_3)*((alpha_r_2+inc_2)/alpha_r_3)*( tet_b_a_blade - tet_b_0_blade);

KTE_vane = ((1 - 0.5 * (gamma-1) * M_2^2*((1/(1-tet_vane))-1))^(-gamma/(gamma-1)) - 1)/(1-(1+0.5 * (gamma-1)*M_2^2)^(-gamma/(gamma-1))); % trailing edge loss coefficient
KTE_blade = ((1 - 0.5 * (gamma-1) * M_r_3^2*((1/(1-tet_blade))-1))^(-gamma/(gamma-1)) - 1)/(1-(1+0.5 * (gamma-1)*M_r_3^2)^(-gamma/(gamma-1)));

% Reynolds Number Calculations

kin_visc = csvread('kinematic_viscosity.csv',1);

Re_vane = (V_2*vane_actual_chord)/ interp1(kin_visc(:,1), kin_visc(:,2), T_2); % Reynolds Number
Re_blade = (Vr_3*blade_actual_chord)/ interp1(kin_visc(:,1), kin_visc(:,2), T_3); 

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
elseif Re_blade > 10e6
    fre_blade = (Re_blade/10e6)^-.2;
else
    fre_blade = 1;
end


% Kt (overall loss before tip losses)

Kt_vane = fre_vane * Kp_vane + Ys_vane  + KTE_vane;
Kt_blade = fre_blade * Kp_blade + Ys_blade + KTE_blade;

% Tip clearance losses and total-to-total efficiency

eta_o_vane = 1-Kt_vane;
eta_o_blade = 1-Kt_blade;

delta_k = blade_tip_clearance*blade_height; 

eta_blade = eta_o_blade - 0.93 * eta_o_blade * (r_t_3/r_m_3) * (delta_k/blade_height*cosd(alpha_r_3));

efficiency = eta_blade*eta_o_vane;
fprintf('Total-to-total efficiency of %4.3f \n',efficiency)
fprintf('Desired efficiency is %4.3f \n', eta_i)

if efficiency < eta_i
    error = 1;
end





