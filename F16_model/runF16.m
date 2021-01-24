%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%                           runF16.m                            %%
%%  Author : Ying Huo                                            %%
%%
%%  This is the main entry of the software package. Users can    %%
%%  select the task he/she like and input appropriate flight     %%
%%  parameters. The results will be shown on the screen promtly. %%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

close all;
clear control_trim state_trim;
% 
exit_code = 1;
while exit_code == 1
    clc;
    disp('        _____         ______                              ');
    disp('        | : \         |    \                              ');
    disp('        | :  `\______|______\_______                      '); 
    disp('         \______              \_____\_____                ');
    disp('          \____/-)_,---------,_____________>              ');
    disp('                     \       /                            ');
    disp('                      |     /                             ');
    disp('                      |____/                              ');
    disp('  ');
    disp('  ');
    disp( ' Welcome to the software package of 6-Degree-of-Freedom F-16 fighter aircraft model !' );
    disp( '  ' );
    disp('  ');
    disp( ' Please make a selection :' );
    disp( ' 1. To get the trim flight data for different flights' );
    disp( ' 2. To derive the linearized model of the aircraft ' );
    disp( ' 3. To simulate the nonlinear response from full nonlinear model' );
    disp(' 4. Quit');
    disp( '  ' );
    index = input( ' which task you want to accomplish? ' );
    disp( '  ' );

    switch index
        
        case 1
            clc;
            disp('  ');
            disp(' You have chosen to compute the trim flight data for flights.');
            disp('  ');
            disp(' Please input the following flight parameters:');
            disp('  ');
            velocity = input(' true velocity ( 350 ~ 1000 ft/sec ) : ');
            disp('  ');
            altitude = input(' altitude ( 0 ~ 40000 ft ) : ');
            disp('  ');
            xcg = input(' center of gravity position ( 0.2 ~ 0.5, reference xcg = 0.35 )  ');
            disp('  ');
            [ state_trim, control_trim ] = find_trim( velocity, altitude, xcg );
            disp('  ');
            disp(' state_trim  = ');
            disp('  ');
            disp( [ '    velocity           ( ft/sec )  = ', num2str( state_trim(1) ) ] );
            disp( [ '    angle of attack    ( rad )     = ', num2str( state_trim(2) ) ] );
            disp( [ '    sideslip angle     ( rad )     = ', num2str( state_trim(3) ) ] );
            disp( [ '    Euler angle phi    ( rad )     = ', num2str( state_trim(4) ) ] );
            disp( [ '    Euler angle theta  ( rad )     = ', num2str( state_trim(5) ) ] );
            disp( [ '    Euler angle psi    ( rad )     = ', num2str( state_trim(6) ) ] );
            disp( [ '    roll rate P        ( rad/sec ) = ', num2str( state_trim(7) ) ] );
            disp( [ '    pitch rate Q       ( rad/sec ) = ', num2str( state_trim(8) ) ] );
            disp( [ '    yaw rate R         ( rad/sec ) = ', num2str( state_trim(9) ) ] );
            disp( [ '    north displacement ( ft )      = ', num2str( state_trim(10) ) ] );
            disp( [ '    east displacement  ( ft )      = ', num2str( state_trim(11) ) ] );
            disp( [ '    altitude           ( ft )      = ', num2str( state_trim(12) ) ] );
            disp( [ '    power       ( percent, 0-100 ) = ', num2str( state_trim(13) ) ] );
            disp('  ');
            disp(' control_trim  = ');
            disp( [ '    throttle           ( 0-1 )     = ', num2str( control_trim(1) ) ] );
            disp( [ '    elevator           ( deg )     = ', num2str( control_trim(2) ) ] );
            disp( [ '    aileron            ( deg )     = ', num2str( control_trim(3) ) ] );
            disp( [ '    rudder             ( deg )     = ', num2str( control_trim(4) ) ] );
            disp('  ');
            exit_code = input('Do you want to continue ( 1 = yes ; 0 = no ) ?  ');
            while exit_code ~= 1 & exit_code ~= 0  
                disp('please input "1" or "0" ' );
                exit_code = input('Do you want to continue ( 1 = yes ; 0 = no ) ?  ');
            end
            if exit_code == 0
                break;
            end;
        
            
        case 2
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
            disp(' You need to input the flight parameter first to get the trim point.');
            disp(' Please input the following flight parameters:');
            disp('  ');
            velocity = input(' true velocity ( 350 ~ 900 ft/sec ) : ');
            disp('  ');
            altitude = input(' altitude ( 0 ~ 40000 ft ) : ');
            disp('  ');
            xcg = input(' center of gravity position ( 0.2 ~ 0.5, reference xcg = 0.35 )  ');
            disp('  ');
            [ A, B, C, D ] = jacobian_f16( velocity, altitude, xcg )
            disp('  ');            
            exit_code = input('Do you want to continue ( 1 = yes ; 0 = no ) ?  ');
            while exit_code ~= 1 & exit_code ~= 0  
                disp('please input "1" or "0" ' );
                exit_code = input('Do you want to continue ( 1 = yes ; 0 = no ) ?  ');
            end
            if exit_code == 0
                break;
            end;

            
        case 3
            clc;
            disp('  ');
            disp(' You have chosen to simulate the nonlinear response from full nonlinear model. ');
            disp('  ');
            disp(' To simulate the F-16 aircraft, you need to choose the flight conditions and get');
            disp(' the trim flight data. The control variables are left as inputs to the block of ');
            disp(' F-16 aircraft model. The software package can do the simulation via simulink model');
            disp(' simF16.mdl and the trajectories of the dynamic responses can be observed in multiple');
            disp(' figure windows.');
            disp('  ');
            disp(' As an example, you will see the response without external control input signals,');
            disp(' i.e., the responses of state and output variables with trimmed control signals.');
            disp(' You need to input the flight parameter first to get the trim point.');
            disp(' Please input the following flight parameters:');
            disp('  ');
            velocity = input(' true velocity ( 350 ~ 1000 ft/sec ) : ');
            disp('  ');
            altitude = input(' altitude ( 0 ~ 40000 ft ) : ');
            disp('  ');
            xcg = input(' center of gravity position ( 0.2 ~ 0.5, reference xcg = 0.35 )  ');
            disp('  ');            
            [ state_trim, control_trim ] = find_trim( velocity, altitude, xcg );
            disp('  ');
            disp( [ ' trimmed angle of attack    ( rad )     = ', num2str( state_trim(2) ) ] );
            disp( [ ' trimmed sideslip angle     ( rad )     = ', num2str( state_trim(3) ) ] );
            disp('  ');
            disp( [ ' trimmed throttle           ( 0-1 )     = ', num2str( control_trim(1) ) ] );
            disp( [ ' trimmed elevator           ( deg )     = ', num2str( control_trim(2) ) ] );
            disp( [ ' trimmed aileron            ( deg )     = ', num2str( control_trim(3) ) ] );
            disp( [ ' trimmed rudder             ( deg )     = ', num2str( control_trim(4) ) ] );
            disp('  ');
            disp('  ');            
            disp('------------------------------------------------');
            disp(' Please wait while the simulation is running ...');
            disp('------------------------------------------------');
            disp(' Note : Because some state variables are almost constants during the trimmed flight,');
            disp(' Matlab will render the axes limits with minimum range allowed by machine precision.');
            disp(' Some warning messages may appear. It does not affect the accuracy of the result. ')
            disp('  ');
            sim( 'simF16_u0', [ 0, 20 ]);
            %% plot the figures
            figure(1);
            subplot(2,2,1); plot(Time, simout(:,1)); title('velocity  vs.  time'); grid;
            ylabel(' velocity (ft/sec)'); xlabel('time (sec)');
            subplot(2,2,2); plot(Time, simout(:,2)); title('angle of attack  vs.  time'); grid;
            ylabel(' angle of attack (rad)'); xlabel('time (sec)');
            subplot(2,2,3); plot(Time, simout(:,3)); title('sideslip angle  vs.  time'); grid;
            ylabel(' sideslip angle (rad)'); xlabel('time (sec)');
            subplot(2,2,4); plot(Time, simout(:,16)); title('mach number  vs.  time'); grid;
            ylabel(' mach number '); xlabel('time (sec)');
            figure(2);   
            subplot(3,1,1); plot(Time, simout(:,4)); title('euler angle phi  vs.  time'); grid;
            ylabel(' phi (rad)'); xlabel('time (sec)');
            subplot(3,1,2); plot(Time, simout(:,5)); title('euler angle theta  vs.  time'); grid;
            ylabel(' theta (rad)'); xlabel('time (sec)');
            subplot(3,1,3); plot(Time, simout(:,6)); title('euler angle psi  vs.  time'); grid;
            ylabel(' psi (rad)'); xlabel('time (sec)');
            figure(3);   
            subplot(3,1,1); plot(Time, simout(:,7)); title('roll rate  vs.  time'); grid;
            ylabel(' roll rate (rad/sec)'); xlabel('time (sec)');            
            subplot(3,1,2); plot(Time, simout(:,8)); title('pitch rate  vs.  time'); grid;
            ylabel(' pitch rate (rad/sec)'); xlabel('time (sec)');
            subplot(3,1,3); plot(Time, simout(:,9)); title('yaw rate  vs.  time'); grid;
            ylabel(' yaw rate (rad/sec)'); xlabel('time (sec)');
            figure(4);  
            subplot(2,2,1); plot(Time, simout(:,10)); title('north distance  vs.  time'); grid;
            ylabel(' north distance (ft)'); xlabel('time (sec)');
            subplot(2,2,2); plot(Time, simout(:,11)); title('east distance  vs.  time'); grid;
            ylabel(' east distance (ft)'); xlabel('time (sec)');
            subplot(2,2,3); plot(Time, simout(:,12)); title('altitude  vs.  time'); grid;
            ylabel(' altitude (ft)'); xlabel('time (sec)');
            subplot(2,2,4); plot(Time, simout(:,13)); title('normal acceleration  vs.  time'); grid;
            ylabel(' normal accelaration (ft/sec^2)'); xlabel('time (sec)');            
%             disp('  ');
%             disp(' To help users get familiar with the how to simulate the aircraft dynamics, the');
%             disp(' software package provides an function to help you implement one-time-step computation.');
%             disp('  ');
%             disp(' For more detailed information, please refer to the manuscript file model.pdf ');
%             disp(' or type "help f16_dynam" in Matlab command window.');
%             disp('  ');
            
            exit_code = input('Do you want to continue ( 1 = yes ; 0 = no ) ?  ');
            while exit_code ~= 1 & exit_code ~= 0  
                disp('please input "1" or "0" ' );
                exit_code = input('Do you want to continue ( 1 = yes ; 0 = no ) ?  ');
            end
            if exit_code == 0
                break;
            end;

        case 4
            exit_code = 0;
    end
end

disp( '  ' );
disp( ' For more information please refer to "readme" and "model.pdf" in the package.' );
disp( '  ' );
disp( ' Thank you. Good bye. ' );
disp( '  ' );


