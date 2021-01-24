%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%                   ydu.m                          %%
%%                                                  %%
%%  Author : Ying Huo                               %%
%%                                                  %%
%% This is to approximate the partial derivative    %%
%% of output y w.r.t. control by computing the      %%
%% difference of the output y(state, control)       %%
%% w.r.t. control at (control+du) and (control-du). %% 
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function  dydu = ydu( x, control, xcg, i, j, du )

% y = [ an; q; alpha ];

time = 0.0;
t = control(j);
control(j) = t - du;  
[ x_dot, an, alat, qbar, amach, q, alpha ] = f16_dynam ( time, x, control, xcg );
y = [ an; q; alpha ];
y1 = y(i);            
control(j) = t + du;  
[ x_dot, an, alat, qbar, amach, q, alpha ] = f16_dynam ( time, x, control, xcg );
y = [ an; q; alpha ];
y2 = y(i);            
dydu = ( y2 - y1 ) / ( du + du ); 
control(j) = t;  % return to original control(j)

