clc;
disp('  ');
disp(' You have chosen to derive the linear model of F-16 aircraft. ');
disp('  ');
disp(' linearized model : x_dot = A * x + B * u' );
disp('                      y   = C * x + D * u');
disp('  ');
disp('  ----State Variables --------                          ');
disp('     x = [ vt    ( ft/sec )    - velocity               ');
disp('           h     ( ft )        - altitude               ');
disp('           alpha ( rad )       - angle of attack        ');  
disp('           theta ( rad )       - Euler angle            ');  
disp('           Q     ( rad/sec )   - pitch rate             ');  
disp('           pow                 - power                  ');  
disp('           beta  ( rad )       - sideslip angle         ');  
disp('           phi   ( rad )       - Euler angle            ');  
disp('           P     ( rad/sec )   - roll rate              ');  
disp('           R     ( rad/sec )   - yaw rate ];            ');  
disp('                                                        ');  
disp('  ---- Control Variables --------                       ');  
disp('     u = [ thtl ( 0 ~ 1.0 )    - throttle setting       ');  
disp('           el   ( deg )        - elevon deflection      ');  
disp('           ail  ( deg )        - aileron deflection     ');  
disp('           rdr  ( deg )        - rudder deflection ];   ');  
disp('                                                        ');  
disp('  ---- Output Variables --------                        ');  
disp('     y = [ an    ( ft/sec^2 )  - normal acceleration;   ');
disp('           q     ( rad/sec )   - pitch rate;            ');
disp('           alpha ( rad )       - angle of attack ]      ');
disp('  ');
disp('---------------------------------------------------------------');
disp('  ');
disp(' Flight parameters to get the trim point:');

altitudes = 5000:5000:40000;
velocities = 350:50:900;  % the model is not accurate when velocity is equal to 350 and altitudes >=35000
[altitudes,velocities] = meshgrid(altitudes,velocities);
[m,n] = size(altitudes);
As = zeros(2,2,m,n);
Bs = zeros(2,1,m,n);
Cs = zeros(1,2,m,n);
Ds = zeros(m,n);
state_trims = zeros(13,m,n);
control_trims = zeros(4,m,n);
pbars = zeros(m,n);

for i = 1:m
    for j = 1:n
        alti = altitudes(i,j);
        veloc = velocities(i,j);
        
%         % for test       
%         alti = 35000;
%         veloc = 350;

        xcg = 0.3;
        [ A, B, C, D,state_trim,control_trim] = jacobian_f16_rev( veloc, alti, xcg );
        fprintf(1,"alti: %5.0f ft, veloc: %3.0f ft/s., trimmed alpha: %5.4f rad, trimmed delta_e: %6.4f deg  \n", alti,veloc,state_trim(2),control_trim(2));
        % extract the short-period (sp) dynamics 
        Asp = A([3,5],[3,5]);
        Bsp = B([3,5],2); % only elevator deflection
        Csp = [1 0];       % angle of attack in radian
        Dsp = 0; 
        state_trim_sp = state_trim([2,8],1);
        control_trim_sp = control_trim(2,1);
        pbar = 0.5*rho_fcn(alti)*veloc^2;
        
        As(:,:,i,j) = Asp;
        Bs(:,:,i,j) = Bsp;
        Cs(:,:,i,j) = Csp;
        Ds(i,j) = Dsp;
        state_trims(:,i,j) = state_trim;
        control_trims(:,i,j)= control_trim;
        pbars(i,j)= pbar;
    end
end 
save('lin_results_tmp.mat')


