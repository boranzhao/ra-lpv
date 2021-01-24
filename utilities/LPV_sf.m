% State-feedback LPV controller design. 
% reference system matrices 
Ar = [0 1; -omega^2  -2*zeta*omega];
Br = [0; omega^2];
Cr = [1 0];


% weighting functions
Wp = 10;

Gplt.A = blkdiag(Ap,Ar);
Gplt.B1 = [0;0;Br];
Gplt.B2 = [Bp;0;0];
Gplt.C1 = [-Wp*Cp Wp*Cr];
Gplt.D1 = 0;
Gplt.D2 = 0;

Theta1_grid = -1:0.1:1;
Theta2_grid = 0;
d_Theta = {[-0.01 0.01],0};
F_theta = @(x) [1 x(1)]; %function for PD matrices
d_F_theta = @(x) [0 1];
FthetaNum = [1 1];

% for avoiding numerical problem
vareps = 1e-5^2; % used for expressing equality constraint
vareps1 = 0; % used for BRL condition and positivity of Lyapunov function

 
if isempty(Theta2_grid) || norm(Theta2_grid) == 0 %% only one GS para theta1
   GSParaNum = 1; 
   d_Theta{2} = 0;
else
    GSParaNum = 2;
end
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
   Ktilde2(Id_Ftheta)=lmivar(2,[nu n]);
   Kff(Id_Ftheta)=lmivar(2,[nu nw]); % feedforward term
end 
    

%% Positivity of Lyapunov function    
for theta1  = Theta1_grid % linspace(min(thetaT),max(thetaT),5) % Gridding may not be necessary because use of affine Laypunov function
    for theta2 = Theta2_grid
        theta = [theta1;theta2];
        Ftheta = F_theta(theta);                    
        lminum = lminum +1;
        for Id_Ftheta=1:sum(FthetaNum)
            lmiterm([lminum 1 1 X(Id_Ftheta)],Ftheta(Id_Ftheta),-1);
        end
%                                  X = X+X_theta(:,:,Id_Ftheta)*Ftheta(Id_Ftheta);                    
        lmiterm([lminum 1 1 0],vareps1); % improving numerical property for simulation
    end
end 

for Id_theta1 = 1:length(Theta1_grid)
    theta1 = Theta1_grid(Id_theta1);
    for Id_theta2 = 1:length(Theta2_grid)
        theta2 = Theta2_grid(Id_theta2);    
        if GSParaNum == 1
            theta = theta1;
        elseif GSParaNum == 2
            theta = [theta1;theta2]; 
        end

        Ga = AugPltEvalSF(Gplt, theta);
        A = Ga.A;
        B1 = Ga.B1; B2 = Ga.B2;
        C1 = Ga.C1; 
        D1 = Ga.D1; D2 = Ga.D1;                    

        Ftheta = F_theta(theta);   
        d_Ftheta = d_F_theta(theta);             
        for Id_d_theta1 = 1:length(d_Theta{1})
            d_theta1 = d_Theta{1}(Id_d_theta1);
            for Id_d_theta2 = 1:length(d_Theta{2}) % d_theta{2} = 0, for one GS para case
                d_theta2 = d_Theta{2}(Id_d_theta2);              
                lminum = lminum+1;  
                for Id_Ftheta=1:sum(FthetaNum)
                    lmiterm([lminum 1 1 X(Id_Ftheta)],Ftheta(Id_Ftheta)*A,1,'s');
                    lmiterm([lminum 1 1 Ktilde2(Id_Ftheta)],Ftheta(Id_Ftheta)*B2,1,'s');
                    lmiterm([lminum 3 1 X(Id_Ftheta)],Ftheta(Id_Ftheta)*C1,1);  
                    lmiterm([lminum 3 1 Ktilde2(Id_Ftheta)],Ftheta(Id_Ftheta)*D2,1);
                    
                    lmiterm([lminum 2 1 -Kff(Id_Ftheta)],1,Ftheta(Id_Ftheta)*B2');
                    lmiterm([lminum 3 2 Kff(Id_Ftheta)],1,Ftheta(Id_Ftheta)*D2);
                end
                for Id_Ftheta = 2:1+FthetaNum(2)
                    lmiterm([lminum 1 1 X(Id_Ftheta)],-d_Ftheta(Id_Ftheta)*d_theta1,1);
%                                         -d_X = -(d_X+X_theta(:,:,Id_Ftheta)*d_Ftheta(Id_Ftheta)*d_theta1); 
                end                                 
                for Id_Ftheta = 2+FthetaNum(2):sum(FthetaNum)%
                    lmiterm([lminum 1 1 X(Id_Ftheta)],-d_Ftheta(Id_Ftheta)*d_theta2,1);
                end
                lmiterm([lminum 2 1 0],B1');            
                
                
                lmiterm([lminum 2 2 Gam],-1,1);
                lmiterm([lminum 3 2 0],D1);
                lmiterm([lminum 3 3 Gam],-1,1);                                
            end %Id_d_theta2
        end %Id_d_theta2
    end % Id_theta2
end % Id_theta1            

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

    for Id_Ftheta = 1:sum(FthetaNum)
        X1(:,:,Id_Ftheta) = dec2mat(lmisys,xopt,X(Id_Ftheta));
        Ktilde1(:,:,Id_Ftheta) = dec2mat(lmisys,xopt,Ktilde2(Id_Ftheta));
        Kff1(:,Id_Ftheta) = dec2mat(lmisys,xopt,Kff(Id_Ftheta));
        
    end
    
    X = X1;
    Ktilde2 = Ktilde1;
    Kff = Kff1;
else
    X = [];
    gam = Inf;
    Ktilde2 = [];
    lmifail = 1;
end
end % nomeaning, just for hiding the code



%% compare the the differences between the reference system and the acutal
% system 
figure;
theta_grid = -1:0.5:1;
for theta = theta_grid
     Ftheta = F_theta(theta);   
     d_Ftheta = d_F_theta(theta); 
     
     Kff_eval = 0;
     X_eval = zeros(2);
     Ktilde2_eval = zeros(1,2);
     for Id_Ftheta=1:sum(FthetaNum)
        X_eval = X_eval + Ftheta(Id_Ftheta)*X(:,:,Id_Ftheta);
        Kff_eval = Kff_eval + Ftheta(Id_Ftheta)*Kff(:,:,Id_Ftheta);
        Ktilde2_eval = Ktilde2_eval + Ftheta(Id_Ftheta)*Ktilde2(:,:,Id_Ftheta);
     end
    
     % Get the close-loop A matrix 
     
     Ap_eval = subs(Ap, theta);
     Bp_eval = subs(Bp,theta);
     Acl = Ap_eval + Bp_eval*Ktilde2_eval/X_eval;
     
     s = tf('s');
     
     TF_actual = Cp/(s*eye(2)-Acl)*Bp_eval*Kff_eval;
     
     
     % evalue the reference system 
     Ar_eval = subs(Ar,theta);
     Br_eval = subs(Br,theta);
     
     
     TF_ref = Cr/(s*eye(2)-Ar_eval)*Br_eval;
     
     subplot(length(theta_grid),1,1);
     bode(TF_actual,TF_ref);  

end


