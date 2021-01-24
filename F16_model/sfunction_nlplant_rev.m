
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%                sfunction_nlplant.m  --  s function            %%
%%                                                               %%
%%  Author : Ying Huo                                            %%
%%                                                               %%
%%  6-DOF full nonlinear dynamics of F-16 aircraft used for      %%
%%  simulink simulation of full nonlinear plant                  %%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function [ sys, x0, str, ts ] = sfunction_nlplant_rev( t, x , u, flag, state_trim, xcg)

switch flag
case 0
   [ sys, x0, str, ts ] = mdlInitializeSizes( state_trim ); % Initialization
   
case 1
   sys = mdlDerivatives( t, x, u, xcg ); % Calculation of derivatives     

case 3
   sys = mdlOutputs(x); % Calculate outputs     
  
case {2,4,9}
   sys = []; % Unused flags
   
otherwise
   error(['Unhandled flag = ', num2str(flag)]); % error handling
end

% end of function sfunct_nlplant   

%=============================================================================
% mdlInitializeSizes
% Return the sizes, initial conditions, and sample times for the S-function.
%=============================================================================
%
function [ sys, x0, str, ts ] = mdlInitializeSizes( state_trim )

sizes = simsizes;

sizes.NumContStates  = 13;
sizes.NumDiscStates  = 0;
sizes.NumOutputs     = 15;  % qbar, amach  
sizes.NumInputs      = 4;   % 4 controls
sizes.DirFeedthrough = 0;   % cganged by Pan from 1 to 0 on Dec 26
sizes.NumSampleTimes = 1;   % at least one sample time is needed

sys = simsizes(sizes);

% initialize the initial conditions
x0  = state_trim;

% str is always an empty matrix for s-function as m-files
str = [];

% initialize the array of sample times to be continuous sample time
ts  = [0 0]; 

% end mdlInitializeSizes

%=============================================================================
% mdlDerivatives
% Return the derivatives for the continuous states.
%=============================================================================
%
function sys = mdlDerivatives( t, x, u, xcg )

% to calculate the state derivatives
[ x_dot, ~, ~, ~,~, ~,~] = f16_dynam ( t, x, u, xcg );
sys = x_dot;

% add the uncertainties to angle of attack and pitch rate
% sys([2 7])= sys([2 7]) +u(5:6); 


% end mdlDerivatives

%=============================================================================
% mdlOutputs
% Return the block outputs.
%=============================================================================
%
function sys = mdlOutputs(x)
% if t==0
%     u = control_trim;
% end
% to calculate the state derivatives
[qbar,amach] = f16_dynam_only_ouput(x);
sys = [  x; qbar; amach ];

% end mdlOutputs
