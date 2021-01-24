close all;
clear;
%% simulation settings
NumParas = 2;
PltType = 2; % Type of plant used for simulation: 1 for LPV, 2 for full nonlinear, 3 for nonlinear (comparing no compensation and partial compensation)
unknown_input_gain = .7;   % note that under unknown input gain, f and C(s)\hat{sigma}(s) will not be close as the later also need to compensate for the unknown input gain. 
enable_L1 = 1;              % whether to enable L1
enable_um_dist_comp = 1;    % whether to enable unmatched disturbance compensation
print_file = 0;             % whether to print the file

% load baseline LPV controller 
load parameters_twoparas_ss_large_omega_range.mat 

% load the settings for unmatched uncertainty conpensation
load unmatched_comp_filter.mat; 
%% sim parameters
sim_t = 6; % Sim time
Ts = 1e-3; % Sim step size passed to Simulink

%% Fit an LPV model for the linearized models
% fit_lpv_model;

%% design a baseline state-feedback LPV controller 
% design_lpv_controller;

%% simulation condition
% veloc: row, alti: colum
% one way to get the trim condition and linearized matrices. 

% idx_row = find(velocities(:,1)==veloc_0,1);
% idx_col = find(altitudes(1,:)==alti_0,1);
% 
% state_trim = state_trims(:,idx_row,idx_col);
% control_trim = control_trims(:,idx_row,idx_col);
% 
% Alti = As(:,:,idx_row,idx_col); Blti = Bs(:,:,idx_row,idx_col);
% control_trims_sorted = reshape(control_trims,4, m*n);
% control_trims_sorted = control_trims_sorted(:,Index);
% state_trims_sorted = reshape(state_trims,13,m*n);
% state_trims_sorted = state_trims_sorted(:,Index);

%% configuration for L1
% fitler 
s = tf('s');
bandwidth = 20;
kDs = bandwidth/s; % which will give a filter kDs/(1+kDs)
[C_num, C_den] = tfdata(kDs,'v');
filter = tf(bandwidth, [1 bandwidth]);
kDs_bar = 200/s;  % 
% prediction error dynamics
Atilde = -10*eye(2); 

% Adaptive Law
sample_time = Ts;

Mat_expm = expm(Atilde*sample_time); % e^(Am*Ts)
Phi = Atilde\(Mat_expm - eye(2));
adapt_gain_no_Bm = -inv(Phi)*Mat_expm;

%% simulate
r2d = 180/pi; % gain for converting radian to degree. 

if PltType == 1
    ref_steps = [1 2 5];
    simdata = cell(1,length(ref_steps));
%         open('simF16_lpvplant.slx');
%         open('simF16_lpvplant_unmatched.slx');
    open('simF16_lpvplant_unmatched.slx');
    duration = 2;
    for  i = 1:length(ref_steps)
        ref_step = ref_steps(i);    
%             sim('simF16_lpvplant.slx');
        sim('simF16_lpvplant_unmatched.slx');
        simdata{i} = simout;        
    end
else
    sim_t = 6;
    ref_step = 2;
    duration = 3.01;
%         
    if PltType == 2
        sim_scenarios = {[40000 400], [20000 600], [5000, 900]}; % different scenarios defined by the altitude and velocity 
        enable_L1_opts =1;
    else
        sim_scenarios = {[20000 600]};            
        enable_um_dist_comp  = 0;    
        enable_L1_opts = [0,1];
    end
    simdata = cell(1,length(sim_scenarios));
    open("simF16_nlinplant_unmatched.slx");
    simu_index = 1;
    for i = 1:length(sim_scenarios)
        alti_0 = sim_scenarios{i}(1);
        veloc_0 = sim_scenarios{i}(2);

        % get trim condition, and linearized matrices (which are only
        % needed for LTI simulation)
        [ A, B, C, D,state_trim,control_trim] = jacobian_f16_rev(veloc_0, alti_0, xcg );
        fprintf(1,"alti: %5.0f ft, veloc: %3.0f ft/s., trimmed alpha: %5.4f rad, trimmed delta_e: %6.4f deg  \n", alti_0,veloc_0,state_trim(2),control_trim(2));
        % % extract the short-period (sp) dynamics 
        Alti = A([3,5],[3,5]);
        Blti = B([3,5],2); % only elevator deflection
        Csp = [1 0];       % angle of attack in radian
        Dsp = 0; 
        for enable_L1 = enable_L1_opts
            sim("simF16_nlinplant_unmatched.slx");
            simdata{simu_index} = simout_nlin;
            simu_index = simu_index+1;
        end
        if PltType == 3
            % add the case for complete compensation
            enable_L1 = 1; enable_um_dist_comp  = 1;
            sim("simF16_nlinplant_unmatched.slx");
            simdata{simu_index} = simout_nlin;
        end

    end   
end

%% plot the results for the simulation
%  set(0,'DefaultAxesFontName','CMU Serif Roman')
% set(0,'defaulttextinterpreter','latex')
linestyle = {'-','-.','--'};
linecolor = ['k', 'r', 'b'];
linewd0 = 1.5;
close all;
if PltType == 1
    for hide= 1
    Time = simdata{1}.pbar_s.Time;
    figure(1)
    % scheduling parameter
    plot(Time,simdata{1}.pbar_s.Data,'k','Linewidth',linewd0); hold on;
    plot(Time,simdata{1}.veloc_s.Data,'-.r','Linewidth',linewd0);
    
    for i=1:length(simdata)  
        % x(1) & r(t) & x_des
        figure(2)
        hold on;
        h_ref = plot(Time,simdata{i}.ref_rad.Data*r2d,':k','Linewidth',1);
        h_x1 = plot(Time,simdata{i}.delta_xsp.Data(:,1)*r2d,'r','Linewidth',linewd0);
        h_x1d = plot(Time,simdata{i}.delta_xsp_des.Data(:,1)*r2d,'--k','Linewidth',linewd0);

        % x(2)
        figure(3)
        hold on;
        h_x2 = plot(Time,simdata{i}.delta_xsp.Data(:,2),'r','Linewidth',linewd0);
        h_x2d = plot(Time,simdata{i}.delta_xsp_des.Data(:,2),'--k','Linewidth',linewd0);

        % ul1(t)
        figure(4)
        hold on;
        plot(Time,simdata{i}.uL1.Data,[linestyle{i},'k'],'Linewidth',linewd0);
        
        figure(6)
        hold on;
        plot(Time,simdata{i}.u_bl.Data,[linestyle{i},'k'],'Linewidth',linewd0);
    end

     % sigma & C sigma
    figure(5)

%     plot(Time,simdata{2}.f.Data,'k-',Time,simdata{2}.uL1_sigma.Data,'b-.','Linewidth',linewd0);
    plot(Time,simdata{2}.sigma.Data(:,1),'k-',Time,simdata{2}.sigmahat.Data(:,1),'r--','Linewidth',linewd0);
    hold on;
    plot(Time,simdata{2}.sigma.Data(:,2),'b-',Time,simdata{2}.sigmahat.Data(:,2),'m--','Linewidth',linewd0);

    % common name
    %  common_str_print= sprintf('lpv_.pdf',veloc_0,alti_0);
    common_str_print = 'lpv_';
    if enable_L1 ==0
        common_str_print = 'lpv_no_comp_';
    elseif enable_um_dist_comp == 0
        common_str_print = 'lpv_matched_comp_';
    end
    

    % add the caption
    figure(1);
    legend({'$\theta_1(t)$', '$\theta_2(t)$'},'Interpreter','latex')
    xlabel('Time (s)')
    grid on;
    goodplot;
    
%     if print_file        
%         print_file_name = [common_str_print 'theta'];
%         savefig(print_file_name);
%         print(print_file_name,'-painters','-dpdf', '-r150')
%     end

    figure(2);
    legend([h_ref, h_x1,h_x1d], {'$r_i(t)$','$x_1(t)$','$x_{1,\mathrm{id}}(t)$'},'Interpreter','latex')
    xlabel('Time (s)')
    ylabel('$\tilde{\alpha}$ (deg)','Interpreter','latex')    
    grid on;
    goodplot;
    if print_file
        print_file_name = [common_str_print 'x1'];
        savefig(print_file_name);
        print(print_file_name,'-painters','-dpdf', '-r150')
    end

    figure(3);
    legend([h_x2,h_x2d], {'$x_2(t)$','$x_{2,\mathrm{id}}(t)$'},'Interpreter','latex')
    xlabel('Time (s)')
    ylabel('$\tilde{q}$ (rad/s)','Interpreter','latex')    
    grid on;
    goodplot;
    if print_file
        print_file_name = [common_str_print 'x2'];
        savefig(print_file_name);
        print(print_file_name,'-painters','-dpdf', '-r150')
    end


    figure(4)
    legend({'$r_1(t)$','$r_2(t)$','$r_3(t)$'},'Interpreter','latex','Location','best')
    xlabel('Time (s)')
    ylabel('$\tilde{\delta}_e: \ u_{ad}$ (deg)','Interpreter','latex')
    grid on;
%     ylim([-9 4]);
    grid on;
    goodplot;
    if print_file        
        print_file_name = [common_str_print 'uL1'];
        savefig(print_file_name);
        print(print_file_name,'-painters','-dpdf', '-r150')
    end
    
    
    figure(6)
    legend({'$r_1(t)$','$r_2(t)$','$r_3(t)$'},'Interpreter','latex','Location','best')
    xlabel('Time (s)')
    ylabel('$\tilde{\delta}_e:\ u_{bl}$ (deg)','Interpreter','latex')
%     ylim([-9 4]);
    grid on;
    goodplot;
    if print_file
        print_file_name = [common_str_print 'ubl'];
        savefig(print_file_name);
        print(print_file_name,'-painters','-dpdf', '-r150')
    end

    figure(5)
    legend({'$\sigma_1(t,x,u)$','$\hat{\sigma}_1(t)$','$\sigma_2(t,x,u)$','$\hat{\sigma}_2(t)$'},'Orientation','Horizontal','Interpreter','latex')
    xlabel('Time (s)')
    ylabel('Uncertainty');
%     ylim([-4 4]);
    grid on;
    goodplot;
%     if print_file
%         print_file_name = [common_str_print 'sigmahat'];
%         savefig(print_file_name);
%         print(print_file_name,'-painters','-dpdf', '-r150')
%     end    
    end
elseif PltType == 2    
    Time = simdata{1}.pbar_s.Time;
    figure(1)
    for i=1:length(simdata)  
        % x(1) & r(t) & x_des
        figure(1)
       
%         h_ref = plot(Time,(simdata{i}.alpha_r.Data-simdata{i}.xsp_trim.Data(:,1))*r2d,':k','Linewidth',0.8);
        h_x1 = plot(Time,simdata{i}.delta_xsp.Data(:,1)*r2d, ['-' linecolor(i)] ,'Linewidth',linewd0*0.8);
        hold on;
        h_x1d = plot(Time,simdata{i}.delta_xsp_des.Data(:,1)*r2d,['--' linecolor(i)],'Linewidth',2);
        
        
         % x(2)
        figure(5)
%         h_ref = plot(Time,(simdata{i}.alpha_r.Data-simdata{i}.xsp_trim.Data(:,1))*r2d,':k','Linewidth',0.8);
        h_x1 = plot(Time,simdata{i}.delta_xsp.Data(:,2), ['-' linecolor(i)] ,'Linewidth',linewd0);
        hold on;
        h_x1d = plot(Time,simdata{i}.delta_xsp_des.Data(:,2),['--' linecolor(i)],'Linewidth',2);
        

        % uL1(t)
        figure(2)
        hold on;
        plot(Time,simdata{i}.uL1_sigma.Data+simdata{i}.uL1_ref.Data,[linestyle{i},linecolor(i)],'Linewidth',linewd0);
        
        figure(3)
        % uBL(t)
        hold on;
        plot(Time,simdata{i}.uBL.Data, [linestyle{i},linecolor(i)],'Linewidth',linewd0);
    end

     % f(t) & C sigma
    figure(4)
    hold on;
%     plot(Time,simdata{1}.sigma.Data(:,1),[linecolor(1) '-.'],Time,-simdata{1}.uL1_sigma.Data,linecolor(i),'Linewidth',linewd0);
%     plot(Time,simdata{1}.sigma.Data(:,1),linecolor(1),Time,simdata{1}.sigmahat.Data(:,1),[linecolor(1) '-.'],'Linewidth',linewd0);
%     plot(Time,simdata{1}.sigma.Data(:,2),linecolor(2),Time,simdata{1}.sigmahat.Data(:,2),[linecolor(2) '-.'],'Linewidth',linewd0);
%     
    % common name
    %  common_str_print= sprintf('lpv_.pdf',veloc_0,alti_0);

    common_str_print = 'nonlinear_';

    figure(1);
    ylim(ylim+[0.2 0]);
    legend({'low $\bar{q}$: $x_1$','low $\bar{q}$: $x_{1,\mathrm{id}}$','med. $\bar{q}$: $x_1$','med. $\bar{q}$: $x_{1,\mathrm{id}}$', 'high $\bar{q}$: $x_1$','high $\bar{q}$: $x_{1,\mathrm{id}}$'},'Interpreter','latex','Location','north','NumColumns',3)
    xlabel('Time (s)')
    ylabel('$\tilde{\alpha}$ (deg)','Interpreter','latex')
    grid on;
    goodplot;
    if print_file
        print_file_name = [common_str_print 'x1'];
        savefig(print_file_name);
        print(print_file_name,'-painters','-dpdf', '-r150')
    end
    
    
    figure(5);
    legend({'low $\bar{q}$: $x_2$','low $\bar{q}$: $x_{2,\mathrm{id}}$','med. $\bar{q}$: $x_2$','med. $\bar{q}$: $x_{2,\mathrm{id}}$', 'high $\bar{q}$: $x_2$','high $\bar{q}$: $x_{2,\mathrm{id}}$'},'Interpreter','latex','Location','best','NumColumns',3)
    xlabel('Time (s)')
    ylabel('$\tilde{q}$ (rad/s)','Interpreter','latex')
    ylim(ylim+[0 0.02]);
    grid on;
    goodplot;
    if print_file
        print_file_name = [common_str_print 'x2'];
        savefig(print_file_name);
        print(print_file_name,'-painters','-dpdf', '-r150')
    end
     
    figure(2);
    legend({'low $\bar{q}$','med. $\bar{q}$', 'high $\bar{q}$'},'Interpreter','latex','Location','northeast')
    xlabel('Time (s)')
    ylabel('$\tilde{\delta}_e: u_{ad}$ (deg)','Interpreter','latex')
    ylim([-12 4]);
    grid on;
    goodplot;
    if print_file
        print_file_name = [common_str_print 'uL1'];
        savefig(print_file_name);
        print(print_file_name,'-painters','-dpdf', '-r150')
    end

       
    figure(3)
    legend({'low $\bar{q}$','med. $\bar{q}$', 'high $\bar{q}$'},'Interpreter','latex','Location','best')
    xlabel('Time (s)')
    ylabel('$\tilde{\delta}_e: u_{bl}$ (deg)','Interpreter','latex')
%     ylim([-9 4]);
    grid on;
    goodplot;
    if print_file
        print_file_name = [common_str_print 'uBL'];
        savefig(print_file_name);
        print(print_file_name,'-painters','-dpdf', '-r150')
    end
    
    figure(4)
%     legend({'$f_1(x,t)$','$\hat{\sigma}_1$','$f_2(x,t)$','$\hat{\sigma}_2$'},'Interpreter','latex')
%     xlabel('Time (s)')
% %     ylim([-4 4]);
%     grid on;
%     goodplot;   
%     if print_file
%         print_file_name = [common_str_print 'sigmahat'];
%         savefig(print_file_name);
%         print(print_file_name,'-painters','-dpdf', '-r150');    
%     end
elseif PltType == 3
    Time = simdata{1}.pbar_s.Time;
    figure(1)
    linecolor = ['k','b','r'];
    
    h_x1d = plot(Time,simdata{1}.delta_xsp_des.Data(:,1)*r2d,[':' 'k'],'Linewidth',2);
    for i=1:length(simdata)  
        % x(1) & r(t) & x_des
        figure(1)
        hold on;
        h_x1 = plot(Time,simdata{i}.delta_xsp.Data(:,1)*r2d, [linestyle{i}, linecolor(i)] ,'Linewidth',linewd0);
        

        % uL1(t)
        figure(2)
        hold on;
        plot(Time,simdata{i}.uL1_sigma.Data+simdata{i}.uL1_ref.Data,[linestyle{i},linecolor(i)],'Linewidth',linewd0);
        
        figure(3)
        % uBL(t)
        hold on;
        plot(Time,simdata{i}.uBL.Data, [linestyle{i},linecolor(i)],'Linewidth',linewd0);
    end

    
    % common name
    %  common_str_print= sprintf('lpv_.pdf',veloc_0,alti_0);

    common_str_print = 'nonlinear_no_comp_';

    figure(1);
    legend({'ideal response', 'no adaptive comp.','comp. for $\hat{\sigma}_m$ only', 'comp. for both $\hat{\sigma}_m$ \& $\hat{\sigma}_{um}$'},'Interpreter','latex','Location','best')
    xlabel('Time (s)')
    ylabel('$\tilde{\alpha}$ (deg)','Interpreter','latex')
    grid on;
    goodplot;
    if print_file
        print_file_name = [common_str_print 'x1'];
        savefig(print_file_name);
        print(print_file_name,'-painters','-dpdf', '-r150')
    end
     
    figure(2);
    legend({'no adaptive comp.','comp. for $\hat{\sigma}_m$ only', 'comp. for both $\hat{\sigma}_m$ \& $\hat{\sigma}_{um}$'},'Interpreter','latex','Location','northeast')
    xlabel('Time (s)')
    ylabel('$\tilde{\delta}_e: u_{ad}$ (deg)','Interpreter','latex')
    ylim(ylim+[0 1]);
    grid on;
    goodplot;
    if print_file
        print_file_name = [common_str_print 'uL1'];
        savefig(print_file_name);
        print(print_file_name,'-painters','-dpdf', '-r150')
    end
       
    figure(3)
    legend({ 'no adaptive comp.','comp. for $\hat{\sigma}_m$ only', 'comp. for both $\hat{\sigma}_m$ \& $\hat{\sigma}_{um}$'},'Interpreter','latex','Location','best')
    xlabel('Time (s)')
    ylabel('$\tilde{\delta}_e: u_{bl}$ (deg)','Interpreter','latex')
%     ylim([-9 4]);
    grid on;
    goodplot;
    if print_file
        print_file_name = [common_str_print 'uBL'];
        savefig(print_file_name);
        print(print_file_name,'-painters','-dpdf', '-r150')
    end   
end


