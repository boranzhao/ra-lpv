%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%              constraint.m                   %%
%%                                             %%
%%  Author : Ying Huo                          %%
%%                                             %%
%%  This is to define & compute the constraint %%
%%  states for different flight conditions.    %%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function  state_constr = constraint( state )

global climb_angle;
global coordi_turn;   % coordinate turn constraint
global stab;          % stability-axis roll constraint  
global skid_turn;     % skidding turn constraint
global rad_gamma;     % flight path angle gamma in radian
global phi_r;
global roll_rate;     % reference roll rate
global pitch_rate;    % reference pitch rate

c_alpha = cos( state(2) );
s_alpha = sin( state(2) );
c_beta  = cos( state(3) );
s_beta  = sin( state(3) );

if coordi_turn
     % coordinated turn logic here
elseif skid_turn ~= 0.0
     % skidding turn logic here
else % non-turning flight here
    state(4) = phi_r;
    d = state( 2 );
    if phi_r ~= 0.0
        d = (-1) * state( 2 ); % inverted
    end
    s_gamma = sin( rad_gamma );
    if s_gamma ~= 0.0 % climbing
        sgocb = sin( rad_gamma ) / c_beta; 
        state(5) = d + atan( sgocb / sqrt( 1.0 - sgocb * sgocb ) ); % rate-of-climbing constraint
    else
        state(5) = d; % level flight
    end
    state(7) = roll_rate;
    state(8) = pitch_rate;
    if stab % stability-axis roll
        state = roll_rate * s_alpha / c_alpha;
    else    % body-axis roll
        state(9) = 0.0;
    end
end
state_constr = state;