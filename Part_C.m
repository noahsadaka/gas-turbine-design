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

% eta_i = 0.87; % desired efficiency

error = 0; % flag for errors

To_1 = 1147.98; % Turbine Inlet Temperature [K]
Po_1 = 315940.86; % Turbine Inlet Pressure [Pa]
m_1 = 4.842; % Turbine Inlet Massflow [kg/s]

To_2 = 1126.67; % Total temperature after bleed air addition [K]
m_2 = 5.00; % Mass flow through blade [kg/s]

To_3 = 972.74; % Total temperature after station [K]
Po_3 = 142188.33; % Total pressure exiting blade, before ITD [Pa]
m_3 = m_2; % Mass flow exiting blade [kg/s]

W = 883830.57; % Stage work [W]

% Gas Constants and Coefficients
Cp = 1148; % Gas specific heat capacity [kg/JK]
gamma = 1.333; % Gas Gamma
Rg = 287; % Gas Constant

% Assumptions 
alpha_3 = 5; % Blade exit swirl angle [deg] Range: -5 to 30
M_3 = 0.3; % Blade exit mach number. Range: 0.3-0.45
R = 0.4; % Reaction at the meanline. 
max_U_h = 1100*0.3048; % Blade speed at meanline [m/s]
inc_1 = 0; % Vane Incidence [deg]
dev_2 = 0; % Vane Deviation [deg]
inc_2 = 0; % Blade Incidence [deg]
dev_3 = 0; % Blade Deviation [deg]


% Given Vane Variables
AR_v = 0.7; % Vane Aspect Ratio
TE_v = 0.045; % Vane TE Thickness [in]

% Given Blade Variables
AR_b = 1.45; % Blade Aspect Ratio
TE_b = 0.025; % Blade TE Thickness [in]

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
A_2 = A_3; % Assumption: roton area is constant
T_2 = R * (T_1 - T_3) + T_3; % Temperature at station 2 [K]
P_2 = R * (P_1 - P_3) + P_3;
Po_2 = P_2*(To_2/T_2)^(gamma/(gamma-1));
rho_2 = P_2/(Rg * T_2); % Density [kg/m3]
V_2 = sqrt(2*Cp*(To_2 - T_2)); % Absolute Velocity [m/s]
M_2 = V_2/sqrt(T_2 * gamma * Rg); % Mach 
Va_2 = m_2 / (rho_2 * A_2); % Axial Velocity [m/s]
Vu_2 = sqrt(V_2^2 - Va_2^2); % Swirl Velocity [m/s]
alpha_2 = atand(Vu_2 / Va_2); % Swirl angle [deg]
U = Ws / (Vu_2 - Vu_3); % Blade Velocity [m/s]

% Station 2 Relative Velocity Triangle
Vru_2 = Vu_2 - U; % Relative swirl velocity [m/s]
Vr_2 = sqrt(Va_2^2 + Vru_2^2); % Relative velocity [m/s]
alpha_r_2 = atand(Vru_2 / Va_2); % Relative swirl angle [deg]
T_r_2 = To_2 - (Vr_2^2)/(2*Cp); % Relative temperature [K]
Po_r_2 = P_2 * (T_r_2 / To_2)^(gamma / (gamma - 1)); % Relative Total Pressure[Pa]
M_r_2 = Vr_2/sqrt(gamma*Rg*T_r_2); % Relative Mach 

% Station 3 Relative Velocity Triangle
Vru_3 = U + Vu_3; % Relative swirl velocity [m/s]
Vr_3 = sqrt(Vru_3^2 + Va_3^2); % Relative Velocity [m/s]
alpha_r_3 = atand(Vr_3/Va_3); % Relative swirl angle [deg]
T_r_3 = To_3 - (Vr_3^2)/(2*Cp); % Relative temperature [K]
Po_r_3 = P_3 * (T_r_3 / To_3)^(-gamma / (gamma - 1)); % Relative Total Pressure [Pa]
M_r_3 = Vr_3/sqrt(gamma*Rg*T_r_3); % Relative Mach 

% double check the velocity reaction
R_check = (Vr_3^2 - Vr_2^2)/(Vr_3^2 - Vr_2^2 + V_2^2 - V_3^2);


if R_check ~= R
    error = 1;
    fprintf('Reactions do not match fuuuuuck\n')
end

if Po_1 > Po_2 && Po_r_2 > Po_r_3
else
    error = 1;
    fprintf('Pressure Physics is not respected Tabarnak\n')
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%% RADII AND STRUCTURAL LIMITATIONS %%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

A_rpm = 1550*A_2; 
AN2_max = 4.5E10; % Max AN2 [in2 rpm 0.5]
N_rpm = sqrt(AN2_max/A_rpm); % Rotation speed [rpm]
%N_rpm = 30000;
N_rads = N_rpm * (1/60) * (2*pi); % Rotation speed [rad/s]

r_m_2 = U/N_rads; % Mean radius at 2 [m]
r_m_3 = r_m_2; % Mean radius at 3 [m]

r_h_2 = (-A_2/pi + 4*r_m_2^2)/(4 * r_m_2); % Hub radius at 2 [m]
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

U_h = U * (r_m_2/r_h_2);

% Station 1
Vu_h_1 = (r_m_1/r_h_1) * Vu_1; % hub swirl vel [m/s]
Va_h_1 = Va_1; % hub axial vel [m/s]
V_h_1 = sqrt(Va_h_1^2 + Vu_h_1^2); % hub velocity [m/s]
alpha_h_1 = atand(Vu_h_1/Va_h_1); % Swirl angle [deg]

% Station 2
Va_h_2 = Va_2;
Vu_h_2 = (r_m_2/r_h_2) * Vu_2; % hub swirl vel [m/s]
V_h_2 = sqrt(Va_h_2^2 + Vu_h_2^2); % hub velocity [m/s]
alpha_h_2 = atand(Vu_h_2/Va_h_2); % Swirl angle [deg]
Vru_h_2 = Vu_h_2 - U_h; % Relative swirl velocity [m/s]
Vr_h_2 = sqrt(Vru_h_2^2 + Va_h_2^2); % Relative velocity [m/s]
alpha_r_h_2 = atand(Vu_h_2/Va_h_2); % Relative swirl angle [m/s]


% Station 3
Va_h_3 = Va_3;
Vu_h_3 = (r_m_3/r_h_3) * Vu_3; % hub swirl vel [m/s]
V_h_3 = sqrt(Va_h_3^2 + Vu_h_3^2); % hub velocity [m/s]
alpha_h_3 = atand(Vu_h_3/Va_h_3); % Swirl angle [deg]
Vru_h_3 = U_h + Vu_h_3; % Relative swirl velocity [m/s]
Vr_h_3 = sqrt(Va_h_3^2 + Vru_h_3^2); % Relative velocity [m/s]
alpha_r_h_3 = atand(Vu_h_3/Va_h_3); % Relative swirl angle [m/s]

R_hub = (Vr_h_3^2 - Vr_h_2^2)/(Vr_h_3^2 - Vr_h_2^2 + V_h_2^2 - V_h_3^2);

if U_h > max_U_h
    fprintf('fuck velocities\n')
    
end
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%% TIP VELOCITY TRIANGLES %%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

U_t = U * (r_m_2/r_t_2);

% Station 1
Vu_t_1 = (r_m_1/r_t_1) * Vu_1; % tip swirl vel [m/s]
Va_t_1 = Va_1; % tip axial vel [m/s]
V_t_1 = sqrt(Va_t_1^2 + Vu_t_1^2); % tip velocity [m/s]
alpha_t_1 = atand(Vu_t_1/Va_t_1); % Swirl angle [deg]

% Station 2
Va_t_2 = Va_2;
Vu_t_2 = (r_m_2/r_t_2) * Vu_2; % tip swirl vel [m/s]
V_t_2 = sqrt(Va_t_2^2 + Vu_t_2^2); % tip velocity [m/s]
alpha_t_2 = atand(Vu_t_2/Va_t_2); % Swirl angle [deg]
Vru_t_2 = Vu_t_2 - U_t; % Relative swirl velocity [m/s]
Vr_t_2 = sqrt(Vru_t_2^2 + Va_t_2^2); % Relative velocity [m/s]
alpha_r_t_2 = atand(Vu_t_2/Va_t_2); % Relative swirl angle [m/s]


% Station 3
Va_t_3 = Va_3;
Vu_t_3 = (r_m_3/r_t_3) * Vu_3; % hub swirl vel [m/s]
V_t_3 = sqrt(Va_t_3^2 + Vu_t_3^2); % hub velocity [m/s]
alpha_t_3 = atand(Vu_t_3/Va_t_3); % Swirl angle [deg]
Vru_t_3 = U_t + Vu_t_3; % Relative swirl velocity [m/s]
Vr_t_3 = sqrt(Va_t_3^2 + Vru_t_3^2); % Relative velocity [m/s]
alpha_r_t_3 = atand(Vu_t_3/Va_t_3); % Relative swirl angle [m/s]

R_tip = (Vr_t_3^2 - Vr_t_2^2)/(Vr_t_3^2 - Vr_t_2^2 + V_t_2^2 - V_t_3^2);



% if error == 1
%     fprintf('Critical Error Somewhere!\n')
% end







