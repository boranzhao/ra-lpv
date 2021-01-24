
function [Opt,stop] = inverse_lpv(Gasym,Ga_ctmat,theta_info, design_para, opt_yalmip)
%Inverse an LPV system based on the following paper.
% Sato, Masayuki. "Gain-scheduled inverse system and filtering system without derivatives of scheduling parameters." In Proc. ACC, vol. 5, pp. 4173-4178. 2003.

% Yalmip is used to formulate the optimization problem, which is solved by sedumi 
% The default GS parameter number is 2, but can also be used for 1 GS parameter, in this case, the following settings should be imposed:
%   Theta2T = {0}, Delta2T = {0}, d_Thetah ={XX,0}, and Fcn_theta is only function of theta1, for instance Fcn_theta = @(x) [1 x(1)]  
% Modified on 10/07/2015 to use same constant X or Y, so that X(i) = X(j) need not be imposed. In this case, under average dwell time switching,
% constant terms are also imposed to be equal at switching surfaces. The contents for ASWCond_Ceq = 0 does not work now. 
% Modifed on 10/04/2015 to includes codes for using LMI Lab
% created on 10/02/2015

% Gasym, generalized plant, created using Matlab symbols
% theta_info a struct includes the following element:
%       ThetaT: cell array, gridded points of theta1 for all regions, ThetaT{1}  for subregion 1, Theta{2} for subregion 2,...
%       Theta2T: cell arry,gridded points of theta2 for all regions
%       Theta1: the variation range of theta1
%       Theta2: the variation range of theta2
%       d_Theta: cell array, bounds for derivatives of theta, d_Theta{1} for theta1, d_Theta{2} for theta2
% design_para a struct includes the following element:
%       XY_PD: 0 for constant X and Y, 1 for PDX and constant Y, 2 for constant X
%           and PD Y, 3 for PD X and Y
%       FthetaNum: a vector, each element show the number of constant matrices and matrices as a function of GS parameters theta1, theta2, ... 
%            for instance, [1 1 1] means one constant matrix, one matrix as a function of theta1, one matrix as a function of theta2,
%            while [1 2] means one constant matrix, two matrices as a function of theta1, no theta2.
%       Fcn_theta: a function handle, Ftheta(theta) will give all the scalar functions for
%                   the matrix variables, for instance in X = X0+
%                   f1(theta)X1+f2(theta)X2, Fcn_theta(theta) = [1 f1(theta) f2(theta)]
%       d_Fcn_theta: a function handle for the derivative of the scalar functions
% opt_yalmip: settings for yalmip
%% Paras 
Theta1 = theta_info.Theta1;
Theta2 = theta_info.Theta2;
d_Theta = theta_info.d_Theta;
gs_num = theta_info.gs_num;

FthetaNum = design_para.FthetaNum;
F_theta= design_para.F_theta;
d_F_theta= design_para.d_F_theta;
LMI_relax = design_para.LMI_relax(2); % use to relax LMIs to guarantee strict feasiblity and avoid numerical problems. 
lminum = 0;
% eps_ri = 1e-5; % to get a solution that is in the interior of the solution set 
% epsi = 1e-; % for imposing equality between two matrices using LMIs 

if norm(Theta2) == 0 || gs_num == 1 %% only one GS para theta1
   gs_num = 1; % Gain scheduling parameter number
   d_Theta{2} = 0; 
else
    gs_num = 2;
end

if XY_PD == 0
    d_Theta = {0,0};
    disp('Derivative of Thetah is set to 0 due to use of constant Lyapunov function');
end
%% Generalized plant parameter
n =  size(Gasym.A,1);
nw = size(Gasym.B1,2);
nu = size(Gasym.B2,2);
nz = size(Gasym.C1,1);
ny = size(Gasym.C2,1); 
nwz = n+nw+nz;
In = eye(n); 
%% subregion parameter 
regnum1 = max(REGID(:,2));
regnum2 = max(REGID(:,3));
regnum = max(REGID(:,1)); % Total num. of subsets. S shape from the bottom to order them, an example:
% 4 5 6
% 1 2 3
% just for positivity of Lyapunov function 
if gs_num == 1
    theta2gridnum = 1;
else
    theta2gridnum = 10;
end    

%% Solve an optimization problem with fixed SSs using YALMIP
% common optimization variables
Pf_theta = sdpvar(n,n,sum(FthetaNum),regnum);
Y_theta = sdpvar(n,n,sum(FthetaNum),regnum);
sig_tmp = sdpvar(1,regnum); %note that using the same sig for different subsets is enough; but using different sigs may improve numerical stability

Gam = sdpvar(1);
gam_tmp = sdpvar(1);
M = sdpvar(gs_num*(n+nw+nz),(gs_num+1)*(n+nw+nz),regnum);
N = sdpvar(gs_num*(n+nz+nw),(gs_num+1)*(n+nz+nw),regnum);   


    Consts = [];
switch XY_PD
    case 0
        X_coeff = [1 zeros(1,sum(FthetaNum)-1)];
        Y_coeff = [1 zeros(1,sum(FthetaNum)-1)];
    case 1
        X_coeff = ones(1,sum(FthetaNum));
        Y_coeff = [1 zeros(1,sum(FthetaNum)-1)];
    case 2
        X_coeff = [1 zeros(1,sum(FthetaNum)-1)];
        Y_coeff = ones(1,sum(FthetaNum));
    case 3
        X_coeff = ones(1,sum(FthetaNum));
        Y_coeff = ones(1,sum(FthetaNum));
end    
    for regid = 1:regnum  
        % note that to impose equal constant Xi or Yi, only the constant
        % from region 1 is utilized
        gam = gam_tmp(regid);  sig = sig_tmp(regid);       
        for i = 1:gs_num
            X_theta(:,:,i+1,regid) = X_theta(:,:,i+1,regid)*X_coeff(i+1);
            Y_theta(:,:,i+1,regid) = Y_theta(:,:,i+1,regid)*Y_coeff(i+1);            
        end
        if XY_PD == 1
            Y_theta(:,:,1,regid) = Y_theta(:,:,1,1)*Y_coeff(1); % use the same constant matrix for all subsets
        elseif XY_PD == 2
            X_theta(:,:,1,regid) = X_theta(:,:,1,1)*X_coeff(1);
        end        
        X0 = X_theta(:,:,1,regid); X1 = X_theta(:,:,2,regid); 
        Y0 = Y_theta(:,:,1,regid); Y1 = Y_theta(:,:,2,regid);        
        if gs_num > 1
            X2 = X_theta(:,:,3,regid);     
            Y2 = Y_theta(:,:,3,regid);
        end
%         if XY_PD == 1
%             Y0 = Y_theta(:,:,1,1)*Y_coeff(1); % use the same constant matrix for all subsets
%             X0 = X_theta(:,:,1,regid)*X_coeff(1);
%         elseif XY_PD == 2
%             X0 = X_theta(:,:,1,1)*X_coeff(1);
%             Y0 = Y_theta(:,:,1,regid)*Y_coeff(1);
%         end
        
        %% Get the grid points for admissible region   
        regid1 = REGID(regid,2); regid2 = REGID(regid,3);
        thetaT = Theta1(regid1,:); theta2T = Theta2(regid2,:); %check the vertices are enough due to parameterically affine LMIs
                
        % Positivity of const Lyapunov function. 
            if XY_PD == 0
                lminum = lminum+1;
                Consts = [Consts [X0 In; In Y0] >= LMI_relax(1)*eye(2*n)];  
            end 
        %% Main LMIs
        for theta1 = thetaT
            for theta2 = theta2T   
%                 theta = [theta1;theta2];               
%                 Ftheta = F_theta(theta);   
%                 d_Ftheta = d_F_theta(theta);                
%                 X= In*0; Y = X;
%                 for Id_Ftheta = 1:length(Ftheta)
%                     X = X+X_theta(:,:,Id_Ftheta,regid)*(Ftheta(Id_Ftheta)*X_coeff(Id_Ftheta)); % parentheses are used to redcue computational burden by implenting scalar multiplication first 
%                     Y = Y+Y_theta(:,:,Id_Ftheta,regid)*(Ftheta(Id_Ftheta)*Y_coeff(Id_Ftheta));
%                 end   
                
                X = X0 + theta1*X1;
                Y = Y0 + theta1*Y1;
                if gs_num > 1
                    X = X+ theta2*X2;
                    Y = Y+ theta2*Y2;
                end                
                % Positivity of parameter-dependent Lyapunov function. 
                if XY_PD ~= 0
                    lminum = lminum +1; 
                    Consts = [Consts [X In; In Y] >= LMI_relax(1)*eye(2*n)];  
                end
                for d_theta1 = d_Theta{1} % d_Thetah = {0,0} is imposed when  XZ_PD ==0, so no problem for using for loop
                    for  d_theta2 = d_Theta{2}% d_Thetah{2} = 0, for one GS para case                               
%                         d_X = In*0; 
%                         d_Y = In*0;                                 
%                         for Id_Ftheta = 2:1+FthetaNum(2)
%                             d_X = d_X+X_theta(:,:,Id_Ftheta,regid)*(d_Ftheta(Id_Ftheta)*d_theta1*X_coeff(Id_Ftheta)); % parentheses are used to redcue computational burden by implenting scalar multiplication first 
%                             d_Y = d_Y+Y_theta(:,:,Id_Ftheta,regid)*(d_Ftheta(Id_Ftheta)*d_theta1*X_coeff(Id_Ftheta)); 
%                         end                                 
%                         for Id_Ftheta = 2+FthetaNum(2):sum(FthetaNum)% if 
%                             d_X = d_X+X_theta(:,:,Id_Ftheta,regid)*(d_Ftheta(Id_Ftheta)*d_theta2*Y_coeff(Id_Ftheta)); % parentheses are used to redcue computational burden by implenting scalar multiplication first 
%                             d_Y = d_Y+Y_theta(:,:,Id_Ftheta,regid)*(d_Ftheta(Id_Ftheta)*d_theta2*Y_coeff(Id_Ftheta)); 
%                         end   
                        %% LMI1 
                        d_X = d_theta1*X1;
                        d_Y = d_theta1*Y1;
                        if gs_num == 2
                            d_X = d_X + d_theta2*X2;
                            d_Y = d_Y + d_theta2*Y2;
                        end
                        clear LMI1;
%                         LMI1 = zeros((gs_num+1)*(n+nw+nz)); % This
%                         line cannot be added otherwise LMI1 will become a
%                         double matrix, to which a sdpvar cannot be
%                         assigned to (otherwise, the matrix inequality
%                         will be taken as a element-wise inequality. 
%                         LMI1(1:nwz,1:nwz) = ...
                        M11=[d_X/2+X0*A_0-sig/2*(C2_0'*C2_0)  zeros(n,nw+nz);
                            B1_0'*X0-sig*D21_0'*C2_0       -gam/2*eye(nw)-sig/2*(D21_0'*D21_0)  zeros(nw,nz);
                            C1_0                              D11_0                           -gam/2*eye(nz)];                        
%                         LMI1(nwz+1:2*nwz,1:nwz) = ... 
                        M21=[X0*A_1+X1*A_0-sig*C2_0'*C2_1                        zeros(n,nw+nz);
                             B1_0'*X1+B1_1'*X0-sig*D21_0'*C2_1-sig*D21_1'*C2_0   -sig*D21_0'*D21_1   zeros(nw,nz);
                             C1_1                                                D11_1               zeros(nz,nz)];                        
%                         LMI1(nwz+1:2*nwz, nwz+1:2*nwz) = ...
                        M22=[X1*A_1-sig/2*(C2_1'*C2_1)        zeros(n,nw+nz);
                             B1_1'*X1-sig*D21_1'*C2_1      -sig/2*(D21_1'*D21_1)      zeros(nw,nz);
                             zeros(nz,nwz)];
                        if gs_num>1
%                             LMI1(2*nwz+1:3*nwz,1:nwz) = ... 
                            M31=[X0*A_2+X2*A_0-sig*C2_0'*C2_2                        zeros(n,nw+nz);
                                 B1_0'*X2+B1_2'*X0-sig*D21_0'*C2_2-sig*D21_2'*C2_0   -sig*D21_0'*D21_2   zeros(nw,nz);
                                 C1_2                                                D11_2               zeros(nz,nz)];
%                             LMI1(2*nwz+1:3*nwz, nwz+1:2*nwz) = ...
                            M32=[X1*A_2+X2*A_1-sig*C2_1'*C2_2        zeros(n,nw+nz);
                                 B1_1'*X2+B1_2'*X1-sig*D21_1'*C2_2-sig*D21_2'*C2_1   -sig*D21_1'*D21_2      zeros(nw,nz);
                                 zeros(nz,nwz)];
%                             LMI1(2*nwz+1:3*nwz, 2*nwz+1:3*nwz) = ...
                            M33=[X2*A_2-sig/2*(C2_2'*C2_2)        zeros(n,nw+nz);
                             B1_2'*X2-sig*D21_2'*C2_2      -sig/2*(D21_2'*D21_2)     zeros(nw,nz);
                             zeros(nz,nwz)];
                            LMI1 = [M11 zeros(nwz,2*nwz);
                                    M21 M22 zeros(nwz,nwz);
                                    M31 M32 M33]; 
                            theta = [theta1 theta2]';
                        else
                            LMI1 = [M11 zeros(nwz,nwz);
                                    M21 M22];
                            theta = theta1;
                                    
                        end
                        LMI1 = LMI1 + kron([theta';-eye(gs_num)],eye(nwz))*M(:,:,regid);                                             
                        LMI1 = LMI1'+LMI1;                         
                        %% LMI2 
%                         LMI2 = zeros((gs_num+1)*nwz);
%                         LMI2(1:nwz,1:nwz) = ...
                        M11=[-d_Y/2+Y0*A_0'-sig/2*(B2_0*B2_0')  zeros(n,nw+nz);
                            C1_0*Y0-sig*D12_0*B2_0'          -gam/2*eye(nz)-sig/2*(D12_0*D12_0')  zeros(nz,nw);
                            B1_0'                              D11_0'                           -gam/2*eye(nw)];                        
%                         LMI2(nwz+1:2*nwz,1:nwz) = ... 
                        M21=[Y0*A_1'+Y1*A_0'-sig*B2_0*B2_1'                        zeros(n,nw+nz);
                             C1_0*Y1+C1_1*Y0-sig*D12_0*B2_1'-sig*D12_1*B2_0'   -sig*D12_0*D12_1'   zeros(nz,nw);
                             B1_1'                                                D11_1'              zeros(nw)];                        
%                         LMI2(nwz+1:2*nwz, nwz+1:2*nwz) = ...
                        M22=[Y1*A_1'-sig/2*(B2_1*B2_1')        zeros(n,nw+nz);
                             C1_1*Y1-sig*D12_1*B2_1'      -sig/2*(D12_1*D12_1')      zeros(nz,nw);
                             zeros(nw,nwz)];
                        if gs_num>1
%                             LMI2(2*nwz+1:3*nwz,1:nwz) = ... 
                            M31=[Y0*A_2'+Y2*A_0'-sig*B2_0*B2_2'                       zeros(n,nw+nz);
                                 C1_0*Y2+C1_2*Y0-sig*D12_0*B2_2'-sig*D12_2*B2_0'   -sig*D12_0*D12_2'   zeros(nz,nw);
                                 B1_2'                                                D11_2'               zeros(nw)];
%                             LMI2(2*nwz+1:3*nwz, nwz+1:2*nwz) = ...
                            M32=[Y1*A_2'+Y2*A_1'-sig*B2_1*B2_2'        zeros(n,nw+nz);
                                 C1_1*Y2+C1_2*Y1-sig*D12_1*B2_2'-sig*D12_2*B2_1'   -sig*D12_1*D12_2'      zeros(nz,nw);
                                 zeros(nw,nwz)];
%                             LMI2(2*nwz+1:3*nwz, 2*nwz+1:3*nwz) = ...
                            M33=[Y2*A_2'-sig/2*(B2_2*B2_2')        zeros(n,nw+nz);
                             C1_2*Y2-sig*D12_2*B2_2'      -sig/2*(D12_2*D12_2')    zeros(nz,nw);
                             zeros(nw,nwz)]; 
                            LMI2 = [M11 zeros(nwz,2*nwz);
                                    M21 M22 zeros(nwz,nwz);
                                    M31 M32 M33]; 
                        else
                            LMI2 = [M11 zeros(nwz,nwz);
                                    M21 M22];
                        end
                        LMI2 = LMI2 + kron([theta';-eye(gs_num)],eye(nwz))*N(:,:,regid);                                             
                        LMI2 = LMI2'+LMI2;                         
                        lminum = lminum + 2;
                        Consts = [Consts LMI1<=-LMI_relax(1)*eye(size(LMI2)) LMI2<=-LMI_relax(1)*eye(size(LMI2))];      
                    end % d_thetah2
                end % d_thetah1       
            end %theta2
        end %theta1  
        lminum = lminum + 1;
        if gam_min_method == 1
            Consts = [Consts gam<= Gam];
        elseif gam_min_method == 2 % not necessary
            if ~isempty(gam_upper)
                Consts = [Consts gam<= gam_upper*1.1];
            end
        end
    end% regid
   
   
    obj1=0; obj2 = 0;
    for j = 1: regnum       
        obj1 = obj1 + (norm(X_theta(:,:,1,j),inf)+norm(Y_theta(:,:,1,j),inf));
        obj2 = obj2 + (norm(M(:,:,j),inf) + norm(N(:,:,j),inf)); %% to avoid ill-condtioned solution (i.e. too large M & N)
    end
    obj1 = obj1/regnum; obj2 = obj2/regnum; obj3 = sum(sig_tmp)/regnum;
    obj = Gam+num_cost_gain*[obj1 obj2 obj3]';
    

    sol = optimize(Consts,obj,opt_yalmip); %+0.1*(obj1+obj2+1e-6*sig) %
    
    if sol.problem ~= 0
        sol.info  
    end 
    if sol.problem ~=0 && sol.problem~=4
        stop = 1;
    end
    %% clear all the intermediate variables
    clear d_X d_Y 
    Opt.X = value(X_theta);
    Opt.Y = value(Y_theta);
    if ss_cond_lc == 1  && SS_DSGN == 1
        Opt.Y_k = Opt.Y;
    end
%     for i=2:regnum
%         if XY_PD == 1
%             Opt.Y(:,:,1,i) = Opt.Y(:,:,1,1); % use same constant matrices for all subsets
%         elseif XY_PD == 2
%             Opt.X(:,:,1,i) = Opt.X(:,:,1,1);
%         end            
%     end    
%     Opt.X(isnan(Opt.X)) = 0;
%     Opt.Y(isnan(Opt.Y)) = 0;     
    Opt.gam = value(gam_tmp);
    Opt.Gam = value(Gam);
    Opt.obj = value(obj);
    Opt.sig = value(sig_tmp);
    Opt.M = value(M);
    Opt.N = value(N);
    Opt.sol = sol;
    
    