% [ A, B, C, D ] = JACOBIAN_F16( velocity, altitude, xcg ) 
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%                    Jacobian_F16.m                       %%
%%                                                         %%
%%  Author : Ying Huo                                      %%
%%                                                         %%
%%  This is to calculate the A, B, C, D Jacobian Matrices  %%
%%  for 6-DOF nonlinear state equations using numerical    %%
%%  linearization. The linearized model is derived         %%
%%  for straight & level flight at  sepecified velocity,   %%
%%  altitude and c.g. position with zero banking angle     %%
%%  where the longitudinal and the lateral subsystems      %%
%%  are decoupled.                                         %%
%%                                                         %%
%%   linearized model : x_dot = A * x + B * u              %%
%%                        y   = C * x + D * u              %%
%%                                                         %%                           
%% ---- State Variables --------                           %%
%%    x = [ vt    ( ft/sec )    - velocity                 %%
%%          h     ( ft )        - altitude                 %%
%%          alpha ( rad )       - angle of attack          %%
%%          theta ( rad )       - Euler angle              %%
%%          Q     ( rad/sec )   - pitch rate               %%
%%          pow                 - power                    %%
%%          beta  ( rad )       - sideslip angle           %%
%%          phi   ( rad )       - Euler angle              %%
%%          P     ( rad/sec )   - roll rate                %%
%%          R     ( rad/sec )   - yaw rate ];              %%
%%                                                         %%
%% ---- Control Variables --------                         %%
%%    u = [ thtl ( 0 ~ 1.0 )    - throttle setting         %%
%%          el   ( deg )        - elevon deflection        %%
%%          ail  ( deg )        - aileron deflection       %%
%%          rdr  ( deg )        - rudder deflection ];     %%
%%                                                         %%
%% ---- Output Variables --------                          %%
%%    y = [ an    ( ft/sec^2 )  - normal acceleration;     %%
%%          q     ( rad/sec )   - pitch rate;              %%
%%          alpha ( rad )       - angle of attack ]        %%
%%                                                         %%
%% ---- Function Inputs --------                           %%
%% velocity ( ft/sec )- true velocity                      %%
%% altitude ( ft )    - altitude                           %%
%% xcg                - center of gravity position as      %%
%%                      fraction of mean aerodynamic chord %% 
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function [ A, B, C, D ] = jacobian_f16( velocity, altitude, xcg ) 

%---- trim the flight ----
[ state_trim, control_trim ] = find_trim( velocity, altitude, xcg );
state_dot_trim = zeros( length( state_trim ), 1 );
state_dot_trim(10) = velocity;

disp('  ');
disp( [ ' trimmed angle of attack    ( rad )     = ', num2str( state_trim(2) ) ] );
disp( [ ' trimmed sideslip angle     ( rad )     = ', num2str( state_trim(3) ) ] );
disp('  ');
disp( [ ' trimmed throttle           ( 0-1 )     = ', num2str( control_trim(1) ) ] );
disp( [ ' trimmed elevator           ( deg )     = ', num2str( control_trim(2) ) ] );
disp( [ ' trimmed aileron            ( deg )     = ', num2str( control_trim(3) ) ] );
disp( [ ' trimmed rudder             ( deg )     = ', num2str( control_trim(4) ) ] );
disp('  ');

%-----------------------------------------
% ---- Compute Jacobian Matrix    ----
%-----------------------------------------
% state = [ Vt; h; alpha; theta; q; pow; beta; phi; p; r ];
% control = [ delta_T; delta_E; delta_A; delta_R ];
% output  = [ an - normal acceleration;
%              q - pitch rate;
%             alpha - angle of attack ]

%---------------------------------------------------------
%   The linerization algorithm chooses smaller and      %
%   smaller perturbations in the independent variable   %
%   and compares three successive approximations to     %
%   the particular partial derivative. If they agree    %
%   within a certain tolerance, then the size of the    %
%   perturbation is reduced to determine if an even     %
%   smaller tolerance can be satisfied. The algorithm   %
%   terminates successfully when a tolerance TOLMIN     %
%   is reached.                                         %
%---------------------------------------------------------

disp('----------------------------------------------------');
disp('Please wait while computing jacobian matrices ......');
disp('----------------------------------------------------');

A = jacobian_A( state_trim, state_dot_trim, control_trim, xcg );
B = jacobian_B( state_trim, state_dot_trim, control_trim, xcg );
C = jacobian_C( state_trim, state_dot_trim, control_trim, xcg );
D = jacobian_D( state_trim, state_dot_trim, control_trim, xcg );