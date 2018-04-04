%clear all; clc;
function deltas = table_4_4(t, state, initial_conds)

    % initial conditions
    % theta0 = 0;
    W      = 2.83176e6; % N
    g      = 9.81; % m/s^2
    m      = W/g; % Kg
    Cw0    = 0.516;
    % u0     = 235.9; % m/s
    rho    = 0.3045; % kg/m^3
    S      = 511; % m^2
    cbar   = 8.324; % m
    Iy     = 0.449e8; % kg m^2

    u0     = initial_conds(1);
    theta0 = initial_conds(4);

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
    Xu    = rho*u0*S*Cw0*sind(theta0) + 0.5*rho*u0*S*Cxu;
    Xw    = 0.5*rho*u0*S*Cxalpha;
    Xq    = 0.25*rho*u0*cbar*S*Cxq;
    Xwdot = 0.25*rho*cbar*S*Cxalphadot;

    Zu    = -rho*u0*S*Cw0*cosd(theta0) + 0.5*rho*u0*S*Czu;
    Zw    = 0.5*rho*u0*S*Czalpha;
    Zq    = 0.25*rho*u0*cbar*S*Czq;
    Zwdot = 0.25*rho*cbar*S*Czalphadot;

    Mu    = 0.5*rho*u0*cbar*S*Cmu;
    Mw    = 0.5*rho*u0*cbar*S*Cmalpha
    Mq    = 0.25*rho*u0*(cbar^2)*S*Cmq;
    Mwdot = 0.25*rho*(cbar^2)*S*Cmalphadot;

    A = [ Xu/m,                         Xw/m,                         0,                                    -g*cos(theta0);
          Zu/(m-Zwdot),                 Zw/(m-Zwdot),                (Zq + m*u0)/(m-Zwdot),                 -m*g*sin(theta0)/(m-Zwdot);
         (Mu + Mwdot*Zu/(m-Zwdot))/Iy, (Mw + Mwdot*Zw/(m-Zwdot))/Iy, (Mq + Mwdot*(Zq + m*u0)/(m-Zwdot))/Iy, -Mwdot*m*g*sin(theta0)/(Iy*(m-Zwdot));
          0                             0                             1                                      0];

    deltas = A*state;
end
