clear all; clc;

% initial conditions
theta0 = 0;
W      = 2.83176e6; % N
g      = 9.81; % m/s^2
m      = W/g; % Kg
Cw0    = 0.516;
u0     = 235.9; % m/s
rho    = 0.3045; % kg/m^3
S      = 511; % m^2
cbar   = 8.324; % m
Iy     = 0.449e8; % kg m^2

% table 6.1
Cxu        = -0.108;
Czu        = -0.106;
Cmu        = 0.1043;
Cxalpha    = 0.2193;
Czalpha    = -4.92;
Cmalpha    = -1.023;
Cxq        = 0;
Czq        = -5.921;
Cmq        = -23.92;
Cxalphadot = 0;
Czalphadot = 5.896;
Cmalphadot = -6.314;

% table 4.4
% Xu    = rho*u0* S * Cw0 * sind(theta0) + 0.5 * rho * u0 * S * Cxu;
% Xw    = 0.50 * rho * u0 * S * Cxalpha;
% Xq    = 0.25 * rho * u0 * cbar * S * Cxq;
% Xwdot = 0.25 * rho * cbar * S * Cxalphadot;

% Zu    = -rho * u0 * S * Cw0 * cosd(theta0) + 0.5 * rho * u0 * S * Czu;
% Zw    = 0.50 * rho * u0 * S * Czalpha;
% Zq    = 0.25 * rho * u0 * cbar * S * Czq;
% Zwdot = 0.25 * rho * cbar * S * Czalphadot;
%
% Mu    = 0.50 * rho * u0 * cbar * S * Cmu;
% Mw    = 0.50 * rho * u0 * cbar * S * Cmalpha;
% Mq    = 0.25 * rho * u0 * (cbar^2) * S * Cmq;
% Mwdot = 0.25 * rho * (cbar^2) * S * Cmalphadot;

% table 6.2
Xu    = -1.982e3;
Xw    =  4.025e3;
Xq    =  0;
Xwdot =  0;

Zu    = -2.595e4;
Zw    = -9.030e4;
Zq    = -4.524e5;
Zwdot =  1.909e3;

Mu    =  1.593e4;
Mw    = -1.563e5;
Mq    = -1.521e7;
Mwdot = -1.702e4;


A = [ Xu/m,                         Xw/m,                         0,                                    -g*cos(theta0);
      Zu/(m-Zwdot),                 Zw/(m-Zwdot),                (Zq + m*u0)/(m-Zwdot),                 -m*g*sin(theta0)/(m-Zwdot);
     (Mu + Mwdot*Zu/(m-Zwdot))/Iy, (Mw + Mwdot*Zw/(m-Zwdot))/Iy, (Mq + Mwdot*(Zq + m*u0)/(m-Zwdot))/Iy, -Mwdot*m*g*sin(theta0)/(Iy*(m-Zwdot));
      0                             0                             1                                      0]

[~, val] = eig(A);
% disp(val);
phugiod = val(4, 4)
wn0   = sqrt(real(phugiod)^2 + imag(phugiod)^2)
zeta0 = -real(phugiod)/wn0
tau0  = 1/(zeta0*wn0)

% page 229
Cxde = -3.818e-6;
Czde = -0.3648;
Cmde = -1.444;
Xde = Cxde * 0.5 * rho * u0^2 * S;
Zde = Czde * 0.5 * rho * u0^2 * S;
Mde = Cmde * 0.5 * rho * u0^2 * S * cbar;

B = [(Xde/m - Xw*Mde/(m*Mw));
     ((Mde*Zw - Mw*Zde)/(m*u0*Mw))]

Apwd = [Xu/m     -g;
       -Zu/(m*u0) 0]


[~, val] = eig(Apwd);
phugiod = val(2, 2)
wn   = sqrt(real(phugiod)^2 + imag(phugiod)^2)
zeta = -real(phugiod)/wn
tau  = 1/(zeta*wn)

fprintf('wn: %f %% change\n',   100*(wn-wn0)/wn0)
fprintf('zeta: %f %% change\n', 100*(zeta-zeta0)/zeta0)
fprintf('tau: %f %% change\n',  100*(tau-tau0)/tau0)

% Lanchester
wn_l   = sqrt(-g*(Zu - (Mu*Zw/Mw))/(m*u0))
zeta_l = -0.5*((g*(Mu*Zw/Mw - Zu)/(m*u0))^-0.5)*(g*Mu/(u0*Mw) + (Xu - Mu*Xw/Mw)/m)

fprintf('wn_l: %f %% change\n',   100*(wn_l-wn0)/wn0)
fprintf('zeta_l: %f %% change\n', 100*(zeta_l-zeta0)/zeta0)
