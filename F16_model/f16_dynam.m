% [ x_dot, an, alat, qbar, amach, q, alpha ] = f16_dynam ( time, x, control, xcg )
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%                                                      %%
%%%      6-DOF full nonlinear model for F-16 aircraft    %%
%%%      [  rigid-body equations referenced to           %%
%%%           body-fixed axis coordinate system  ]       %%
%%%                                                      %%
%%%      Author : Ying Huo                               %%
%%%                                                      %%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%                                                                    %%
%%% ---- State Variables --------                                      %%
%%% x = [ vt    ( ft/sec )    - velocity                               %%  
%%%       alpha ( rad )       - angle of attack                        %%
%%%       beta  ( rad )       - sideslip angle                         %%
%%%       phi   ( rad )       - Euler angle                            %%
%%%       theta ( rad )       - Euler angle                            %%
%%%       psi   ( rad )       - Euler angle                            %%
%%%       P     ( rad/sec )   - roll rate                              %%
%%%       Q     ( rad/sec )   - pitch rate                             %%
%%%       R     ( rad/sec )   - yaw rate                               %%
%%%       integral of north speed    ( ft )  - north displacement      %%
%%%       integral of east speed     ( ft )  - east displacement       %%    
%%%       integral of vertical spped ( ft )  - altitude                %%
%%%       pow  ( percent, 0 <= pow <= 100 )   - power ];               %%
%%%                                                                    %%  
%%% ---- control Variables --------                                    %%   
%%% u = [ thtl ( 0 <= thtl <= 1.0 ) - throttle                         %%
%%%       el   ( deg )              - elevator                         %%
%%%       ail  ( deg )              - aileron                          %% 
%%%       rdr  ( deg )              - rudder ];                        %%
%%%                                                                    %%
%%% ---- parameters ---------                                          %% 
%%% xcg - center of gravity position as                                %%
%%%       fraction of mean aerodynamic chord                           %% 
%%%                                                                    %%
%%% ---- output Variables --------                                     %% 
%%% output = [ x_dot                - 1st order derivative of state x  %%
%%%            an    ( ft/sec^2 )   - normal acceleration              %% 
%%%            alat  ( ft/sec^2 )   - lateral acceleration in y-axis   %%
%%%            qbar  ( psf )        - dynamic pressure                 %%
%%%            amach                - mach number                      %% 
%%%            q     ( rad/sec )    - pitch rate                       %%
%%%            alpha ( rad )        - angle of attack ];               %%  
%%%                                                                    %%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%



function  [ x_dot, an, alat, qbar, amach, q, alpha ] = f16_dynam ( time, x, control, xcg )

%% ---- constant variable ----------------
s = 300;
b = 30;
cbar = 11.32;
rm = 1.57e-3; % 1 / mass
xcgr = 0.35;
he = 160.0;

%% ---- Inertia constants ----------------
c1 = -0.770;
c2 = 0.02755;
c3 = 1.055e-4;
c4 = 1.642e-6;
c5 = 0.9604;
c6 = 1.759e-2;
c7 = 1.792e-5;
c8 = -0.7336;
c9 = 1.587e-5;

rtod = 180 / pi; %% radians to degrees
g = 32.174;

%% ---- Assign state variables ---------------
vt = x(1);
alpha = x(2) * rtod; %% x(2) in radians, alpha in degrees.
beta = x(3) * rtod;  %% x(3) in radians, beta in degrees.
phi = x(4);
theta = x(5);
psi = x(6);
p = x(7);
q = x(8);
r = x(9);
alt = x(12);
pow = x(13); % power

%% ---- Assign state & control variables ---------------
thtl = control(1);
el   = control(2);
ail  = control(3);
rdr  = control(4);

%% ---- Air Data computer and engine model ------------
[tfac, t, rho, amach, qbar, ps ] = adc ( vt, alt ); 
cpow = tgear (thtl);
x_dot (13) = pdot ( pow, cpow );  %% x_dot(13) = power derivative
t = thrust ( pow, alt, amach );

%% ---- Look-up tables and component biuldup ------------
cxt = cx ( alpha, el );
cyt = cy ( beta, ail, rdr );
czt = cz ( alpha, beta, el );
dail = ail / 20.0;
drdr = rdr / 30.0;
dlda_value = dlda( alpha, beta );
dldr_value = dldr( alpha, beta );
clt = cl( alpha, beta ) + dlda_value * dail + dldr_value * drdr;
cmt = cm( alpha, el );
dnda_value = dnda( alpha, beta );
dndr_value = dndr( alpha, beta );
cnt = cn( alpha, beta ) + dnda_value * dail + dndr_value * drdr;

%% ---- Add damping derivatives ------------
tvt = 0.5 / vt;
b2v = b * tvt;
cq = cbar * q * tvt;
d = damping( alpha );
cxt = cxt + cq * d(1);
cyt = cyt + b2v * ( d(2) * r + d(3) * p );
czt = czt + cq * d(4);
clt = clt + b2v * ( d(5) * r + d(6) * p );
cmt = cmt + cq * d(7) + czt * ( xcgr - xcg );
cnt = cnt + b2v * ( d(8)* r + d(9) * p ) - cyt * ( xcgr - xcg ) * cbar / b;

%% ---- Get ready for state equations ------------
cbta = cos( x(3) );
u = vt * cos( x(2) ) * cbta;
v = vt * sin( x(3) );
w = vt * sin( x(2) ) * cbta;
sth = sin( theta );
cth = cos( theta );
sph = sin( phi );
cph = cos( phi );
spsi = sin( psi );
cpsi = cos( psi );
qs = qbar * s;
qsb = qs * b;
rmqs = rm *qs;
gcth = g * cth;
qsph = q * sph;
ay = rmqs * cyt;
az = rmqs * czt;

%% ---- Force equations ------------
udot = r * v - q * w - g * sth + rm * ( qs * cxt + t );
vdot = p * w - r * u + gcth * sph + ay;
wdot = q * u - p * v + gcth * cph + az;
dum = ( u * u + w * w );
x_dot(1) = ( u * udot + v * vdot + w * wdot ) / vt;
x_dot(2) = (u * wdot - w * udot ) / dum;
x_dot(3) = (vt * vdot - v * x_dot(1) ) * cbta / dum;

%% ---- Kinematics ------------
x_dot(4) = p + ( sth / cth ) * (qsph + r * cph );
x_dot(5) = q * cph - r * sph;
x_dot(6) = ( qsph + r * cph ) / cth;

%% ---- Moments ------------
x_dot(7) = ( c2 * p + c1 * r + c4 * he ) * q + qsb * ( c3 * clt + c4 * cnt );
x_dot(8) = ( c5 * p -c7 * he )* r + c6 * ( r^2 - p^2 ) + qs * cbar * c7 * cmt;
x_dot(9) = ( c8 * p - c2 * r + c9 * he ) * q + qsb * ( c4 * clt + c9 * cnt );

%% ---- Mavigation ------------
t1 = sph * cpsi;
t2 = cph * sth;
t3 = sph * spsi;
s1 = cth * cpsi;
s2 = cth * spsi;
s3 = t1 * sth - cph * spsi;
s4 = t3 * sth + cph * cpsi;
s5 = sph * cth;
s6 = t2 * cpsi + t3;
s7 = t2 * spsi - t1;
s8 = cph * cth;

x_dot(10) = u * s1 + v * s3 + w * s6;  %% north speed
x_dot(11) = u * s2 + v * s4 + w * s7;  %% east speed
x_dot(12) = u * sth - v * s5 - w * s8; %% vertical speed

x_dot = x_dot'; % transfer the vector to a column vector

%% ---- Outputs ------------
an = (-1) * az / g;
alat = ay / g;


