%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%               jacobian_B.m                            %%
%%                                                       %%
%%  Author : Ying Huo                                    %%
%%                                                       %%
%% This is to compute the Jacobian matrix B of the       %%
%% nonlinear state equations of 6-DOF aircraft model.    %%
%%   linear model : x_dot = A * x + B * control          %%
%%                    y   = C * x + D * control          %%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function B_reduce = jacobian_B( state_eq, state_dot_eq, control_eq, xcg )

%-----------------------------------------
% ---- Compute Jacobian Matrix    ----
%-----------------------------------------

%---- size of the Jacobian matrix ----
row = 1 : length(state_eq);      % array corresponding to rows of the Jacobian matrix
column = 1 : length(control_eq); % array corresponding to columns of the Jacobian matrix
no_row = length( row );          % No. of rows of the Jacobian matrix
no_col = length( column );       % No. of columns of the Jacobian matrix

%---- parameters ----
del = 0.01;
dmin = 0.5;
tolmin = 3.3e-5;
oktol = 8.1e-4;
min_dv = 0.0001;

ij= 1;
answer = [];
count = 0;
for j = 1 : no_col
    dv0 = max( abs(del * state_eq(column(j)) ), dmin );
    dv1 = 1.1 * dv0;
    dv2 = 1.2 * dv0;
    abc = [];
    for i = 1 : no_row
        tol = 0.1;
        a0 = fdu( state_eq, control_eq, xcg, row(i), column(j), dv0 );
        a1 = fdu( state_eq, control_eq, xcg, row(i), column(j), dv1 );
        a2 = fdu( state_eq, control_eq, xcg, row(i), column(j), dv2 );
        b0 = min( abs(a0), abs(a1) ); % b0 = min( |a0|, |a1| )
        b1 = min( abs(a1), abs(a2) );
        d0 = abs( a0 - a1 );          % d0 = | a1 - a0 |
        d1 = abs( a2 - a1 );
        k = 0;
        while k <= 100 
            new_dv = 0.6 * dv0;
            if new_dv <= min_dv * state_eq( column(j) )
                break;
            else
                a2 = a1;
                a1 = a0;
                a0 = fdu( state_eq, control_eq, xcg, row(i), column(j), new_dv );
                b1 = b0;
                b0 = min( abs(a0), abs(a1) );
                d1 = d0;
                d0 = abs( a0 - a1 );
                if d0 <= tol * b0 & d1 <= tol * b1
                    tol = tol * 0.8;
                    if tol <= tolmin
                        tol = tolmin;
                        break;
                    end
                end
                k = k + 1;                    
            end
        end
        
        ans = a1;
        count = count + 1;
        abc = [ abc, ans ];
    end
    answer = [ answer;  abc];
end
ij = ij + 1;
answer = answer';
B_org = answer;

%-----------------------------------------
%----      reorganize the matrix      ----
%-----------------------------------------

%-- delete uninterested states --
%% rows %%
answer(6, :)=[];     % euler angle psi
answer(10-1, :)=[];  % north displacement
answer(11-2, :)=[];  % east displacement

% -----------------------------------------
%   organize the matrix to decouple the
%    longitudinal and lateral motions 
% -----------------------------------------
% state = [ Vt; h; alpha; theta; q; pow; beta; phi; p; r ];

%% rows %%
[ m, n ] = size( answer );
answer_B = zeros( m, n );
answer_B( 1, : ) = answer( 1, : );
answer_B( 2, : ) = answer( 9, : );
answer_B( 3, : ) = answer( 2, : );
answer_B( 4, : ) = answer( 5, : );
answer_B( 5, : ) = answer( 7, : );
answer_B( 6, : ) = answer( 10, : );
answer_B( 7, : ) = answer( 3, : );
answer_B( 8, : ) = answer( 4, : );
answer_B( 9, : ) = answer( 6, : );
answer_B( 10, : ) = answer( 8, : );

%% ---- output of the jacobian matrix ----
B_reduce = answer_B;
% count;

