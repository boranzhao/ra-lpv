% [ an, qbar, amach ] = f16_dynam ( time, x, control, xcg )
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%                                                      %%
%%%      6-DOF full nonlinear model for F-16 aircraft    %%
%%%      [  rigid-body equations referenced to           %%
%%%           body-fixed axis coordinate system  ]       %%
%%%                                                      %%
%%%      Author : Ying Huo (edited by Pan Zhao)          %%
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
%%% output = [ an    ( ft/sec^2 )   - normal acceleration              %% 
%%%            qbar  ( psf )        - dynamic pressure                 %%
%%%            amach                - mach number                      %% 
%%%                                                                    %%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%



function  [qbar, amach] = f16_dynam_only_ouput(x)


%% ---- Assign state variables ---------------
vt = x(1);
alt = x(12);


%% ---- Air Data computer and engine model ------------
[tfac, t, rho, amach, qbar, ps ] = adc ( vt, alt ); 
% cpow = tgear (thtl);
% t = thrust ( pow, alt, amach );
