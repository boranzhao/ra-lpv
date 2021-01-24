% [ state_trim, control_trim ] = find_trim( velocity, altitude, xcg )
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%                  find_trim.m                           %%
%%                                                        %%
%%  Author : Ying Huo                                     %%
%%                                                        %%
%%  This is the function to find the trim condition for   %%
%%  straight&level flight at given velocity and altitude. %%
%%  The trim data are derived by optimizing the cost      %%
%%  function cost_f16.m using simplex search.             %%
%%                                                        %%
%% ---- Input Variables --------                          %%
%% velocity (ft/sec)- true velocity                       %%
%% altitude (ft)    - altitude                            %%
%% xcg              - center of gravity position as       %%
%%                    fraction of mean aerodynamic chord  %% 
%%                                                        %%
%% ---- Output Variables --------                         %%
%% state_trim   - trim data for the state vector where    %%
%%                                                        %%
%%                state vector =                          %%
%%                                                        %%
%%                [ velocity             ( ft/sec )       %%  
%%                  angle of attack      ( rad )          %%
%%                  sideslip angle       ( rad )          %%
%%                  phi   - Euler angle  ( rad )          %%
%%%                 theta - Euler angle  ( rad )          %%
%%%                 psi   - Euler angle  ( rad )          %%
%%                  roll rate            ( rad/sec )      %%
%%                  pitch rate           ( rad/sec )      %%
%%                  yaw rate             ( rad/sec )      %%
%%                  north displacement   ( ft )           %%
%%                  east displacement    ( ft )           %%    
%%                  altitude             ( ft )           %%
%%                  power  ( percent, 0 <= pow <= 100 ) ] %%
%%                                                        %%
%% control_trim - trim data for the control vector where  %%
%%                                                        %%
%%                control vector =                        %%
%%                                                        %%
%%                [ throttle setting ( 0 - 1 )            %%
%%                  elevon           ( deg )              %%
%%                  aileron          ( deg )              %% 
%%                  rudder           ( deg ) ];           %%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function [ state_trim, control_trim ] = find_trim( velocity, altitude, xcg )

%---- flight constraints --------------------------
global climb_angle; % rate-of-climbing constraint
% climb_angle = input( 'Input the climb angel(deg)' );
climb_angle = 0.0;
global coordi_turn; % coordinate turn constraint
coordi_turn = 0;
global stab; % stability-axis roll constraint  
stab = 0;
global skid_turn; % skidding turn constraint
skid_turn = 0;
global rad_gamma; % flight path angle gamma in radian
rad_gamma = 0;
global phi_r; % reference phi
phi_r = 0;
global roll_rate; % reference roll rate
roll_rate = 0;
global pitch_rate; % reference pitch rate
pitch_rate = 0;

    
%---- data ---------------------------
rtod = 57.29577951; % radian to degree

no_step = 5000; % no. of iteration steps for trimming
disp('  ');
read_no = input('Please input the # of trim iterations ( default = 5000 ):  ');
if read_no == []
    no_step = read_no;
end
disp('  ');
disp('----------------------------------------------------');
disp('Please wait while computing the trim data ......');
disp('----------------------------------------------------');
epsilon = -1.0;

%---- initial condition ----
x0 = [ velocity
       0.0
       0.0
       0.0
       0.0
       0.0
       0.0
       0.0
       0.0 
       0.0 
       0.0 
       altitude
       90 ];
    
 u0 = [  0.73
        -1.0
         0.0  
         0.0 ]; 

%---- define initial simplex --------------------
s = [ u0(1); u0(2); x0(2); u0(3); u0(4); x0(3) ];
ds = [ 0.2; 1.0; 0.02; 1.0; 1.0; 0.02 ];

%---- simplex algorithm -----------------
init_cost =  cost_f16( x0, u0, s, xcg );
[ s_trim, f_final ] = simplex( s, ds, x0, u0, no_step, epsilon, xcg );

%---- output the trim result ------------
control_trim = u0;
control_trim(1) = s_trim(1);
control_trim(2) = s_trim(2);
control_trim(3) = s_trim(4);
control_trim(4) = s_trim(5);
state = x0;
state(2) = s_trim(3);
state(3) = s_trim(6);
final_cost = cost_f16( state, control_trim, s_trim, xcg );
state(13) = tgear( control_trim(1) );
state_trim = constraint( state );
control_trim;

