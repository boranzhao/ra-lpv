%% State-feedback LPV controller design with pole placement constraints

% version history,
% 01/02/2020: edited for two scheduling parameters
% 12/23/2019: added the pole placement approach

% note: the natural frequency of the nominal closed-loop system is
% determined by scaled dynamic pressue.

% weighting functions
Wp = 10;
Cp = [1 0];

%% 
load lpv_model.mat
zeta = 0.7;  % damping ratio
omega_lower_bnd = 2;  
omega_upper_bnd = 4;
omega = omega_lower_bnd+ (pbar_s + 1)/2*(omega_upper_bnd-omega_lower_bnd);% natural frequency 


Gplt.A = blkdiag(Ap);
Gplt.B1 = [0;0];
Gplt.B2 = Bp;
Gplt.C1 = -Wp*Cp;
Gplt.D1 = Wp;
Gplt.D2 = 0;

Theta_grid = -1:0.2:1;   % range of scaled dynamic pressure)
Theta2_grid = -1:0.4:1;  % range of scaled aircraft speed. 
d_Theta = [-0.1 0.1];
d_Theta2 = [-0.1 0.1];
F_theta = @(x) [1 x(1) x(2)]; %function for PD matrices

d_F_theta = @(x) [0 1 1];
FthetaNum = [1 1 1];

% % for parameter-dependent X
F_thetaX = F_theta; %function for PD matrices
d_F_thetaX = d_F_theta; %function for PD matrices

% for constant X
% F_thetaX = @(x) [1 0]; 
% d_F_thetaX = @(x) [0 0];

% for avoiding numerical problem
vareps1 = 0; % used for BRL condition and positivity of Lyapunov function

%% Generalized plant parameter  
n =  size(Gplt.A,1);
nw = size(Gplt.B1,2);
nu = size(Gplt.B2,2);
nz = size(Gplt.C1,1);
I = eye(n); 

for nomeaning = 1
lminum = 0;
%% Defining opt variables
setlmis([]);
Gam = lmivar(1,[1 1]);
for Id_Ftheta=1:sum(FthetaNum)
   X(Id_Ftheta)=lmivar(1,[n 1]);
   L(Id_Ftheta)=lmivar(2,[nu n]);
%    Kff(Id_Ftheta)=lmivar(2,[nu nw]); % feedforward term
end 

%% Positivity of Lyapunov function    
for theta1  = Theta_grid % linspace(min(thetaT),max(thetaT),5) % Gridding may not be necessary because use of affine Laypunov function
    for theta2 = Theta2_grid  
        theta=[theta1;theta2];
        FthetaX = F_thetaX(theta);                    
        lminum = lminum +1;
        for Id_Ftheta=1:sum(FthetaNum)
            lmiterm([lminum 1 1 X(Id_Ftheta)],FthetaX(Id_Ftheta),-1);
        end
    %                                  X = X+X_theta(:,:,Id_Ftheta)*Ftheta(Id_Ftheta);                    
        lmiterm([lminum 1 1 0],vareps1); % improving numerical property for simulation

    end 
end

for id_theta1 = 1:length(Theta_grid)
    theta1 = Theta_grid(id_theta1);
    for id_theta2 = 1:length(Theta2_grid)        
        theta = [theta1; Theta2_grid(id_theta2)];

        Ga = AugPltEvalSF(Gplt,theta);
        A = Ga.A;
        B1 = Ga.B1; B2 = Ga.B2;
        C1 = Ga.C1; 
        D1 = Ga.D1; D2 = Ga.D1;                    

        Ftheta = F_theta(theta); 
        FthetaX = F_thetaX(theta);   

        d_Ftheta = d_F_theta(theta);
        d_FthetaX = d_F_thetaX(theta);
        for id_d_theta1 = 1:length(d_Theta)
            d_theta1 = d_Theta(id_d_theta1);  
            for id_d_theta2 = 1:length(d_Theta2)
                d_theta2 = d_Theta2(id_d_theta2);                  
                lminum = lminum+1;  
                for id_Ftheta=1:sum(FthetaNum)
                    lmiterm([lminum 1 1 X(id_Ftheta)],FthetaX(id_Ftheta)*A,1,'s');
                    lmiterm([lminum 1 1 L(id_Ftheta)],Ftheta(id_Ftheta)*B2,1,'s');
                    lmiterm([lminum 3 1 X(id_Ftheta)],FthetaX(id_Ftheta)*C1,1);  
                    lmiterm([lminum 3 1 L(id_Ftheta)],Ftheta(id_Ftheta)*D2,1);

        %             lmiterm([lminum 2 1 -Kff(id_Ftheta)],1,Ftheta(id_Ftheta)*B2');
        %             lmiterm([lminum 3 2 Kff(id_Ftheta)],1,Ftheta(id_Ftheta)*D2);
                end
                for id_Ftheta = 2:1+FthetaNum(2)
                    lmiterm([lminum 1 1 X(id_Ftheta)],-d_FthetaX(id_Ftheta)*d_theta1,1);  % -d_X = -(d_X+X_theta(:,:,id_Ftheta)*d_Ftheta(id_Ftheta)*d_theta1); 
                end        
                for id_Ftheta = 1+FthetaNum(2)+1:sum(FthetaNum)
                    lmiterm([lminum 1 1 X(id_Ftheta)],-d_FthetaX(id_Ftheta)*d_theta2,1);                                       
                end    
                
                lmiterm([lminum 2 1 0],B1');            


                lmiterm([lminum 2 2 Gam],-1,1);
                lmiterm([lminum 3 2 0],D1);
                lmiterm([lminum 3 3 Gam],-1,1);          
            end % id_d_theta2
        end %id_d_theta1


        % for pole placement, the poles are constrained by the following
        % constraints

        % Constraints to make sure zeta*omega_n (real part of poles) is in 
        % a interval by using 
        % <AX+B2*L>+ 2*alpha*X <0;   (equation (20) in [Chilali 1996])
       
        % desired poles only depends on the dynamic pressure.
        % offset_low and scale are tuning parameters to generate the
        % desired closed-loop poles (confirmed after ploting the poles of
        % the resulting closed-loop system) controller design).
        offset_low = 4.7;
        scale = 2.5;
        alpha= zeta*(scale*theta1+offset_low ); 
        lminum = lminum + 1;
        for id_Ftheta=1:sum(FthetaNum)
             lmiterm([lminum 1 1 X(id_Ftheta)],FthetaX(id_Ftheta)*A,1,'s');
             lmiterm([lminum 1 1 L(id_Ftheta)],Ftheta(id_Ftheta)*B2,1,'s');
             lmiterm([lminum 1 1 X(id_Ftheta)],2*alpha*FthetaX(id_Ftheta),1);
        end
        % <AX+B2*L>+ 2*alpha*X >0
        offset_upp = offset_low+1;
        alpha = zeta*(scale*theta1+offset_upp);
        lminum = lminum + 1;
        for id_Ftheta=1:sum(FthetaNum)
             lmiterm([-lminum 1 1 X(id_Ftheta)],FthetaX(id_Ftheta)*A,1,'s');
             lmiterm([-lminum 1 1 L(id_Ftheta)],Ftheta(id_Ftheta)*B2,1,'s');
             lmiterm([-lminum 1 1 X(id_Ftheta)],2*alpha*FthetaX(id_Ftheta),1);
        end


        % constraints to make sure omega_n in a circle
        % not a LMI region; thus cannot be enforced
    %     (-rX (A+B2*K)X; X(A+B2*K)^T -rX) < 0 (equation (21) in [Chilali 1996])
    %     r = theta+3+10;
    %     lminum = lminum + 1;
    %     for id_Ftheta=1:sum(FthetaNum)
    %          lmiterm([lminum 1 1 X(id_Ftheta)],-r*FthetaX(id_Ftheta),1);
    %          lmiterm([lminum 2 1 -L(id_Ftheta)],Ftheta(id_Ftheta),B2');
    %          lmiterm([lminum 2 1 X(id_Ftheta)],FthetaX(id_Ftheta),A');
    %          lmiterm([lminum 2 2 X(id_Ftheta)],-r*FthetaX(id_Ftheta),1);
    %     end

       

        % constraints to enforce damping ratio
        % [sin(theta)(AclX+XAcl') cos(theta)(AclX-XAcl');
        %  cos(theta)(XAcl'-AclX) sin(theta)(AclX+XAcl')]<0 (equation (22)
        %  in [Chilali 1996])
        beta = atan(sqrt(1/zeta^2-1)); %(the angle about real axis corresponding to the desired damping ratio
        lminum = lminum + 1;
        for id_Ftheta=1:sum(FthetaNum)
             lmiterm([lminum 1 1 X(id_Ftheta)],A,FthetaX(id_Ftheta)*sin(beta),'s');
             lmiterm([lminum 1 1 L(id_Ftheta)],B2,Ftheta(id_Ftheta)*sin(beta),'s');

             lmiterm([lminum 2 1 X(id_Ftheta)],FthetaX(id_Ftheta)*cos(beta),A');
             lmiterm([lminum 2 1 -L(id_Ftheta)],Ftheta(id_Ftheta)*cos(beta),B2');
             lmiterm([lminum 2 1 X(id_Ftheta)],A, -FthetaX(id_Ftheta)*cos(beta));
             lmiterm([lminum 2 1 L(id_Ftheta)],B2,-Ftheta(id_Ftheta)*cos(beta));

             lmiterm([lminum 2 2 X(id_Ftheta)],A,FthetaX(id_Ftheta)*sin(beta),'s');
             lmiterm([lminum 2 2 L(id_Ftheta)],B2,Ftheta(id_Ftheta)*sin(beta),'s');
        end
%         
    end        
end % id_theta1            

lmisys = getlmis;
% formulating the cost function
nvar = decnbr(lmisys); %number of decision variable.
c = zeros(nvar,1);
c(1)=1; 
% LMI Lab solver setting
options(1)= 1e-4; % relative accuary on the optimal value
options(2)= 1000; %Number of Iteration
options(4) =  10; % J, the code terminates when the objective has not decreased by more than the desired relative accuracy during the last J iterations
options(5)= 0; % 1 for not showing the process

% solve the optimization problem
[Gam,xopt] = mincx(lmisys,c,options);

% check whether all constraints are satisfied
% lmifail = 0;
%         evals = evallmi(lmisys,xopt);        
%         for i = 1: lminum
%             [lhs,rhs] = showlmi(evals,i);
%             if max(real(eig(lhs-rhs))) > 5e-7 
%                 eig_max = max(real(eig(lhs-rhs)))
%                 disp ('Not all LMI constranits are satisfied')
%                 lmifail = 1;
% %                 break;
%             end
%         end 
%         lmifail
if ~isempty(Gam)        

    for id_Ftheta = 1:sum(FthetaNum)
        X1(:,:,id_Ftheta) = dec2mat(lmisys,xopt,X(id_Ftheta));
        L1(:,:,id_Ftheta) = dec2mat(lmisys,xopt,L(id_Ftheta));
%         Kff1(:,id_Ftheta) = dec2mat(lmisys,xopt,Kff(id_Ftheta));
        
    end
    
    X = X1;
    L = L1;
%     Kff = Kff1;
else
    X = [];
    gam = Inf;
    L = [];
    lmifail = 1;
end
end % nomeaning, just for hiding the code


%% check the pole location of the closed-loop system
% system 
close all
figure;
theta1_grid = -1:0.2:1;
theta2= 1;
poles_all_cl = zeros(2,length(theta1_grid));
poles_all_ol = zeros(2,length(theta1_grid));

for ind = 1:length(theta1_grid)
     theta1 = theta1_grid(ind);
     theta = [theta1;theta2];
     Ftheta = F_theta(theta);   
     
    
     X_eval = zeros(2);
     L_eval = zeros(1,2);
     for id_Ftheta=1:sum(FthetaNum)
        X_eval = X_eval + Ftheta(id_Ftheta)*X(:,:,id_Ftheta);
%         Kff_eval = Kff_eval + Ftheta(id_Ftheta)*Kff(:,:,id_Ftheta);
        L_eval = L_eval + Ftheta(id_Ftheta)*L(:,:,id_Ftheta);
     end
    
     % Get the close-loop A matrix      
     Ap_eval = double(subs(Ap,{pbar_s,veloc_s},{theta(1),theta(2)}));
     Bp_eval = double(subs(Bp,{pbar_s,veloc_s},{theta(1),theta(2)}));
     Acl = Ap_eval + Bp_eval*L_eval/X_eval;   
     
     subplot(2,1,1);
     TF_op = ss(Ap_eval,Bp_eval,Cp,0);
     poles = pole(TF_op);
     poles_all_ol(:,ind) = poles;
     scatter(real(poles),imag(poles)); hold on;
     xlabel('Real Axis');
     ylabel('Imaginary Axis');
     title('Open-loop system poles vs $\bar{q}_s$','interp','latex');
     
     
     subplot(2,1,2)
     TF_cl = ss(Acl,Bp_eval,Cp,0);   

     poles = pole(TF_cl);
     poles_all_cl(:,ind) = poles;
     fprintf(1,'ps_bar: %.3f,omega: %.3f pole bounds: [%.3f  %.3f],actual poles: %.3f, %.3f\n',theta1, abs(poles(1)), -0.7*(scale*theta1+offset_upp), -0.7*(scale*theta1+offset_low),poles); 
     scatter(real(poles),imag(poles)); hold on;
     xlabel('Real Axis');
     ylabel('Imaginary Axis');
     title('Closed-loop system poles vs $\bar{q}_s$','interp','latex');
end

% plot the lines for the desired damping ratio
beta = atan(sqrt(1/zeta^2-1));

x = -4:0.1:0;
y = x*tan(beta);

subplot(2,1,1)
plot(x,y,'k--',x,-y,'k--');
goodplot([7 5]);
legend off
subplot(2,1,2)
plot(x,y,'k--',x,-y,'k--');
goodplot([7 5]);
legend off

% print('OL_CL_poles','-painters','-dpdf', '-r150');  
%% A different way to plot;
close all;
colormap(jet);
c_tmp = jet;
indices = ceil(linspace(1,256,length(theta1_grid)));
c =  c_tmp(indices,:);
% c = jet(9);
figure; 
subplot(2,1,1)
hold on;
for ind = 1:9
% plot(real(poles_all_ol(:,ind)),imag(poles_all_ol(:,ind)),'o','LineWidth',2,'Color',c(ind,:))
scatter(real(poles_all_ol(:,ind)),imag(poles_all_ol(:,ind)),'o','MarkerEdgeColor',c(ind,:))
end
xlabel('Real Axis');
ylabel('Imaginary Axis');
% title('Open-loop system poles vs $\bar{q}_s$','interp','latex');
% sgrid;
hcb = colorbar('Ticks',[0 0.5 1],'TickLabels',{'-1','0','1'});%,'location','northoutside'
title(hcb,'$\bar{q}_s$','interp','latex');

subplot(2,1,2)
hold on;
for ind = 1:9
% scatter(real(poles_all_cl(:,ind)),imag(poles_all_cl(:,ind)),'o','LineWidth',2,'Color',c(ind,:))
scatter(real(poles_all_cl(:,ind)),imag(poles_all_cl(:,ind)),'o','MarkerEdgeColor',c(ind,:))
end
xlabel('Real Axis');
ylabel('Imaginary Axis');
% title('Open-loop system poles vs $\bar{q}_s$','interp','latex');
% sgrid;
hcb = colorbar('Ticks',[0 0.5 1],'TickLabels',{'-1','0','1'}); %,'Direction','Reverse'
colormap(jet);
title(hcb,'$\bar{q}_s$','interp','latex');

beta = atan(sqrt(1/zeta^2-1));

x = -5:0.1:0;
y = x*tan(beta);

subplot(2,1,1)
plot(x,y,'k--',x,-y,'k--');
xlim([-3 0]);
goodplot([7 5]);
legend off
subplot(2,1,2)
plot(x,y,'k--',x,-y,'k--');
xlim([-5 0]);
goodplot([7 5]);
legend off
colormap(jet);
%%
% print('OL_CL_poles','-painters','-dpdf', '-r150');  


