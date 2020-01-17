%% =========================================================== 
%   MTOP-aIVG                                                         
%   (Multiresolution Topology OPtimization adaptive   
%   Isosurface Variable Grouping)
% 
%                                                        by Jaeeun Yoo
%
%                                                 ver. 10.2s ('20.1.17) 
%  
% --------------------------------------------------------------------
% Compliance minimization with volume fraction constraint
% Increasing penalization, Adaptive threshold projection filter
%
% Input parameter :
% MTOP_aIVG_R10_2s(XL,YL,ZL,elem_size,FE_order,rmin_m,volfrac,data_N,DL_ratio,filter_id,penal_min,elem_type,ex_type,solver_type,logon_id)
% XL : physical length (X), YL : Y length, ZL : Z length (meter)
% elem_size : physical FE mesh size (rectangular grid, meter)
% FE_order : shape function order (FE analysis basis function, 1:1st, 2:2nd)
% rmin_m : radius of regularization filter (meter)
% volfrac : volume fraction (0~1)
% data_N : IVG reducing variable number (0:No IVG)
%          0 = no grouping
%          1 = auto grouping
%          N = fixed grouping number (target number is fixed, effective
%          group number is changing )
% DL_ratio : MTOP ratio , ex) when MTOP ratio of 3, number of density element
%                                  is 3 times of FE elements
% filter_id : 1=Sensitivity filter, 2=Density filter, 
%             3=Adaptive heaviside projection filter (ATP)
% penal_min : minimal penalizae value for continuation scheme (max = 3.0)
% elem_type :  'QUAD4'  : 4 node 2D elem, static
%              'QUAD8'  : 8 node 2D elem, 2nd order, static
%              'QUAD4_HT': 4 node 2D elem, heat transfer
%              'BRICK8' : 8 node 3D brick elem, static
%              'BRICK27': 27 node 3D brick elem, 2nd order, static
% ex_type : numerical example name
% solver_type : FEA solver type, 'Direct' or 'PCG'
% logon_id : 0 = No display
%            1 = Plot all figures
%            2 = No plot, but save some figures
%            log on : 1, log off : 0, iteration history & iteration plot, 2 : short log / no plot
% =======================================================================

function MTOP_aIVG_R10_2s(XL,YL,ZL,elem_size,FE_order,rmin_m,volfrac,data_N,DL_ratio,filter_id,penal_min,elem_type,ex_type,solver_type,logon_id)
close all
tic


%% SW parameters
density_th=0.2; % display density threshold

%% Optimization parameters

% SIMP parameters
penal = penal_min;
penal_max=3.0;

group_variable_coef = 1;

penal_th=0.0001; % penalty threshold = default 0.05, update : 0.005 ('19.06.30)
penal_inc=1.1; % penalty value increment

data_N_high_lim = 1/2; % maximum group number, limit, ratio
data_N_low_lim_val = 100; % minimum group number, limit, value

rmin=(rmin_m/elem_size*DL_ratio);
volfrac_ini = volfrac; % initial volume fraction

% Optimizer parameters
tolx=0.01; % convergence criteria
kkt_tol=1e-4; % Convergence criteria for MMA, KKTnorm
tol_obj=0.0001; % convergenc criteria for object function
maxloop=2000;

exitflag=0; % stop condition, exitflag : 0=iteration limit, 1=design variable change tol. 2=objective function decrease done

%% Filter & MMA parameters

Mnd_c=0.0001; % gray level indicator variation criteria for projection filter parameter beta, % unit, '19.03.12
beta = 0.0; % projection filter parameter, initial, '19.02.17, 0.01
beta_max=5.0; % projection filter parameter max (beta),  '19.02.17
beta_inc=1.5; % projection filter parameter increment (beta = beta0 + beta_inc^(k-1)), '19.02.19
eta=0.5; % threshold projection filter parameter (eta)

%% Variable initialize

addpath([pwd,'\\FEA_func']); % FEA functions directory
addpath([pwd,'\\MTOP_func']); % MTOP functions directory
addpath([pwd,'\\IVG_func']); % IVG functions directory
addpath([pwd,'\\Optimizer']); % Optimizer functions directory
addpath([pwd,'\\Filter']); % Filter functions directory
addpath([pwd,'\\Plot_func']); % Plot functions directory

iter_ini=0;
t_total_ini=0;

%% FEA parameters

load([ex_type,'.mat']); % load FEA parameters
int_scale=FE_order+1; % integration point number per one density element, number of gauss point

nelx=round(XL/elem_size);
nely=round(YL/elem_size);
nelz=0;
nelx_d=round(nelx*DL_ratio);
nely_d=round(nely*DL_ratio);
nelz_d=0;

x(1:nely_d,1:nelx_d) = volfrac_ini;
nele=nelx_d*nely_d;

debug_data=zeros(maxloop,10); % debugging data for optimizer
time_ansys=zeros(12,1);


%% Design / non-design domain
[x,x_d,x_nd_s,x_nd_v,new_volfrac_ini]=Define_design_domain_r2(nelx_d,nely_d,nelz_d,x,volfrac_ini,elem_type,ex_type);

if data_N==0 % data_N==0 -> No isosurface variables grouping
    data_N2 = size(x_d,1);
elseif data_N==1 % data_N==1 -> Auto group variable
    data_N2 = round(size(x_d,1) * data_N_high_lim);
else 
    data_N2 = data_N; % fixed group variable (effective variable num fixed)
end



%% SW INITIALIZE

if filter_id==1
    filter_id_str='Sensitivity filter';
elseif filter_id==2
    filter_id_str='Density filter';
elseif filter_id==3
    filter_id_str='Heaviside projection filter(tanh)';
elseif filter_id==4
    filter_id_str='Heaviside projection filter(exp)';
else
    filter_id_str='No filter';
end


switch ex_type
    case 'HT_ex1'        
        DOF_num = 1; %'19.03.31, 2D heat transfer, DOFs per node
    case 'HT_ex2'
        DOF_num = 1; %'19.03.31, 2D heat transfer, DOFs per node
    otherwise
        DOF_num = 2; %'19.03.31, 2D static, DOFs per node
       
end

if FE_order == 1
    u_ini = zeros(DOF_num*(nelx+1)*(nely+1),1);
elseif FE_order == 2
    u_ini = zeros(DOF_num*(2*nelx+1)*(2*nely+1),1);
end
disp('== Parameters (2D topology MTOP) ==================')
disp(['FEM model : ',ex_type])
disp(['X length : ',num2str(XL),' m'])
disp(['Y length : ',num2str(YL),' m'])
disp(['FE element size : ',num2str(elem_size),' m'])
disp(['FE element order : ',num2str(FE_order),' (linear/quadratic) (',elem_type,')'])
disp(['Number of FE element (x) : ',num2str(nelx)])
disp(['Number of FE element (y) : ',num2str(nely)])
disp(['Number of Top element (x) : ',num2str(nelx_d)])
disp(['Number of Top element (y) : ',num2str(nely_d)])
disp(['Design Variables : ',num2str(size(x_d,1)),' -> ',num2str(data_N2)]);
disp(['Filter type : ',filter_id_str]);
disp(['Filter size : ',num2str(rmin),' element'])
disp(['Volume fraction : ',num2str(volfrac)])
disp(['Force : ',num2str(Fnew),' N / theta1 = ',num2str(theta1),' deg / theta2 = ',num2str(theta2),' deg'])
disp(['Stop criteria (max element density change) : ',num2str(tolx)])
disp(['Stop loop limit (max iteration) : ',num2str(maxloop)])
disp(['Initial penalize value : ',num2str(penal_min)])
disp('====================================================')



%% Input output variables

new_dir=['2D_',elem_type,'_Nel_',num2str(nelx),'_',num2str(nely),'_Sele_',num2str(elem_size),'_r_',num2str(rmin_m),'_Nd_',num2str(data_N),'_DLratio_',num2str(DL_ratio),'_Vf_',num2str(volfrac),'_P_',num2str(penal_min),'_to_',num2str(penal_max),'_p_inc_',num2str(penal_inc),'_F',num2str(filter_id),'_',solver_type];
list_dir = dir;
folder_exist_id = 0;
for j=1:size(list_dir,1)
    if strcmp(list_dir(j).name,'result') && (list_dir(j).isdir == 1)
        folder_exist_id = 1;
    end
end
if folder_exist_id == 1
else
    mkdir('result');
end
mkdir([pwd,'\result'],new_dir);
save_path=[[pwd,'\result\'],new_dir,'\'];

inp_var = struct('XL',XL,...
    'YL',YL,'elem_size',elem_size,'FE_order',FE_order,'rmin_m',rmin_m,...
    'data_N',data_N,'theta1',theta1,'theta2',theta2,'penal_min',penal_min,...
    'DL_ratio',DL_ratio,'int_scale',int_scale,...
    'volfrac',volfrac,...
    'E_min',E_min,'E',E,'filter_id',filter_id,...
    'tolx',tolx,'maxloop',maxloop,'Fnew',Fnew,'elem_type',elem_type);



%% PREPARE FILTER

if filter_id==1 || filter_id==2 || filter_id==3 || filter_id==4
    
    inp_filter_prep = struct('nele',nele,'rmin',rmin,...
        'nelx_d',nelx_d,'nely_d',nely_d);
    [H,Hs]=filter_prep_r2(1,inp_filter_prep,x_d);
    
end


%% MMA optimizer initialize

[m,n,xmin,xmax,xold1,xold2,low,upp,a0,a,c_MMA,d] = MMA_param_initialize(data_N2,new_volfrac_ini,1);
xold1_xd = new_volfrac_ini*ones(nele,1);
xold2_xd = new_volfrac_ini*ones(nele,1);
xold1_xd(x_nd_s) = 1;
xold2_xd(x_nd_s) = 1;
xold1_xd(x_nd_v) = 0;
xold2_xd(x_nd_v) = 0;

xval = new_volfrac_ini*ones(data_N2,1);

%% I matrix calculation
disp('I_mat calculation...')

% calculate I_mat 
I_mat = I_mat_func_r2(elem_size,DL_ratio,int_scale,nu,elem_type);

% filter initialize value
x_1d=x(:);
if filter_id == 1 || filter_id == 2
    x_phys_d = x_1d(x_d);
    xTilde = x_1d(x_d);
elseif filter_id == 3
    xTilde = x_1d(x_d);    
    if beta == 0
        x_phys_d = xTilde;
    else
        x_phys_d = (tanh(beta*eta)+tanh(beta*(xTilde-eta)))/(tanh(beta*eta)+tanh(beta*(1-eta))); % threshold projection filter
    end
    
elseif filter_id == 4
    xTilde = x_1d(x_d);
    x_phys_d = 1-exp(-beta*xTilde)+xTilde*exp(-beta); % heaviside step function    
else
    x_phys_d = x_1d(x_d);
    xTilde = x_1d(x_d);
end
    

disp('Filter prepared, I_mat prepared')

Mnd2_old=100; % gray level indicator initial 
Mnd2=100; % gray level indicator initial

loop=iter_ini;
change = 1;

beta_inc_cnt = 0;
beta_inc_cnt_old = 0;

penal_inc_cnt = 0;
penal_inc_cnt_old = 0;

if filter_id == 1 || filter_id == 2
    beta = beta_max;
end

beta_cnt = 0; % beta incement count
% loopbeta = 0;

change_obj = 1;
c_old = 1;

kktnorm = 1;
data_N_tg = data_N;

penal_cnt = 0;
dc_f_normalize_max_old = 0; % dc_f_normalize_max, previous iter. value 

loop_MMA = 0;

while ((kktnorm > kkt_tol && change_obj > tol_obj) || beta_inc_cnt_old > 0 ||...
        penal_inc_cnt_old > 0 || penal < penal_max || beta < beta_max) && loop < maxloop % '19.03.10, penal max criteria
        
  loop = loop + 1;
  loop_MMA = loop_MMA + 1;
 
  loop_s_time=toc;
 
  
  %% Increasing penalize method
  %penal_inc=1.1;
  penal_cnt=penal_cnt+1;
  if (change_obj < penal_th || penal_cnt>=50) && penal < penal_max
      penal_cnt=1;
      loop_MMA = 1;
      
      penal=penal*penal_inc;
      if penal>penal_max
          penal=penal_max;
      end
      disp(['Penalize value is changed : ',num2str(penal)])
            
      penal_inc_cnt = penal_inc_cnt + 1;
      
  end
  
  
  
  %% variable discretization, create sensitivity map
  % gray level value
  Mnd2=sum(4*x_phys_d(:).*(1-x_phys_d(:)))/length(x_phys_d(:))*100; % When x equals 0 or 100, Mnd goes to 0 value
  
  if filter_id==3 || filter_id==4
      
      if loop>=2 && penal>=penal_max    % modified, addition of penal condition, '18.07.19         
          if (abs(Mnd2_old-Mnd2)/abs(Mnd2_old) < Mnd_c) && beta<beta_max
            
              if beta == 0
                  beta = 1.0;
              else
                  beta = beta*beta_inc;
              end

              beta_cnt = beta_cnt + 1;

              if beta>beta_max
                  beta=beta_max;
              end
              disp(['beta parameter is changed to ',num2str(beta)])
              
              beta_inc_cnt = beta_inc_cnt+1;
          end
      end
      
  end
    
  Mnd2_old=Mnd2;
  
  %% data_N2 changes
  if loop > 1 && data_N > 0 
      
      if data_N == 1

          xmma_min = min(xmma_temp);
          xmma_max = max(xmma_temp);
          xmma_band = xmma_max-xmma_min;
          
          data_N_tg = group_variable_coef*(xmma_band*0.80)/tolx;  % '19.08.21
      end      
      
      if group_num_estimation < (data_N_tg/2) && data_N2 < data_N_high_lim*size(x_d,1)
          data_N_old = data_N2;
          data_N2 = 2 * data_N2;
          
      elseif group_num_estimation > (2*data_N_tg) && data_N2 > (data_N_low_lim_val)
          data_N_old = data_N2;
          data_N2 = fix(data_N2/2);
      else
          data_N_old = data_N2;
      end
      
      if data_N2 > data_N_old || data_N2 < data_N_old
          [m,n,xmin,xmax,xold1,xold2,low,upp,a0,a,c_MMA,d,xval] = MMA_param_update(data_N_old,data_N2,xold1,xold2,low,upp,xval,1);
          disp(['eff. data_N = ',num2str(size(eff_x_d,1))])
          disp(['data_N changed : ',num2str(data_N_old),' -> ',num2str(data_N2)])          
          
      end
  end
  
  
  %% Generate Global Stiffness (GK) matrix
  
  x_phys_FE_1d=x(:);
  x_phys_FE_1d(x_d)=x_phys_d;
  x_phys_FE_1d(x_nd_s)=1;
  x_phys_FE_1d(x_nd_v)=0;  
  
  x_phys_FE=reshape(x_phys_FE_1d,nely_d,nelx_d);
  nelz = 0;
  
  [KG]=Gen_GK_func_r3(nelx,nely,nelz,x_phys_FE,penal,E,E_min,DL_ratio,I_mat,elem_type);
    
  %disp('KG prepared')
  time_ansys(1,1)=time_ansys(1,1)+toc-loop_s_time;
  
  %% Solve FEA
  loop_s_time=toc;
    
  if FE_order==1
      switch ex_type          
          case 'HMBB'
              [U]=FE5_2D_lin_HMBB(nelx,nely,KG,Fnew,theta1,u_ini,solver_type); % '18.05.07, Half MBB ex         
              
      end  
  end
  
    %% u_ini estimation (1st order taylor approximation)
  if strcmp(solver_type,'PCG')
      u_ini = U;
  else
      u_ini = U;      
  end

  %disp('FE(U) prepared')
  time_ansys(2,1)=time_ansys(2,1)+toc-loop_s_time;
  
  loop_s_time=toc;
  
  
  %% calculate compliance (FE analysis)
  
  
  if FE_order==1
      [c,dc]=cal_compliance_2D_lin_r4(x_phys_FE,nelx,nely,penal,U,E,E_min,DL_ratio,I_mat); % '19.06.23, 2D linear
      
  elseif FE_order==2
      [c,dc]=cal_compliance_2D_quad_r2(x_phys_FE,nelx,nely,penal,U,E,E_min,DL_ratio,I_mat); % '19.02.15, 2D quad
  end
  
  
  if loop == 1 % initial compliance for MMA normalization / New penal iteration
      c_ini = abs(c);
      disp(['compliance initialized (C initial =  ',num2str(c_ini),')'])
  end
  
  time_ansys(3,1)=time_ansys(3,1)+toc-loop_s_time;
  
  loop_s_time=toc;
  
  % sensitivity filter
  dc_1d=dc(:);
  x_1d=x(:);
  %x_1d=x_phys_d(:);
  if filter_id==1
      dc_f = H*(x_1d(x_d).*dc_1d(x_d))./Hs./max(1e-3,x_1d(x_d));
  elseif filter_id==2
      dc_f = H*(dc_1d(x_d)./Hs);
  elseif filter_id==3
      %dx = beta*exp(-beta*xTilde)+exp(-beta); % heaviside step function
      if beta == 0
          dc_f = H*(dc_1d(x_d)./Hs);
%           dc_f_xTilde = dc_f;
      else          
          dx = beta*(1-tanh(beta*(xTilde-eta)).^2)/(tanh(beta*(1-eta))+tanh(eta*beta)); % threshold projection filter
          dc_f = H*(dc_1d(x_d).*dx(:)./Hs);
%           dc_f_xTilde = H*(dc_1d(x_d)./Hs);
      end
      
  elseif filter_id==4
      dx = beta*exp(-beta*xTilde)+exp(-beta); % heaviside step function
      %dx = beta*(1-tanh(beta*(xTilde-eta)).^2)/(tanh(beta*(1-eta))+tanh(eta*beta)); % threshold projection filter
      dc_f = H*(dc_1d(x_d).*dx(:)./Hs);
      
  elseif filter_id==5
      dc_f = dc_1d(x_d);
      
  end
  
  
    
  time_ansys(4,1)=time_ansys(4,1)+toc-loop_s_time;
  
  loop_s_time=toc;
 
  
%% variable group criteria

if data_N==0
        
else  
    
    if filter_id==1  % sensitivity filter
        dc_f_norm_temp = dc_f;        
      
    else 
        dc_f_norm_temp = H*(dc_1d(x_d)./Hs);

    end
    
    dc_f_norm_idx_minus = find(dc_f_norm_temp<0);
    dc_f_norm_idx_plus = find(dc_f_norm_temp>=0);
        
    
    dc_f_normalize_max = max(dc_f_normalize_max_old,max(abs(dc_f_norm_temp))); 
    dc_f_normalize_min = 0; 
    
    dc_f_normalize = dc_f_norm_temp;
    dc_f_normalize(dc_f_norm_idx_plus) = (dc_f_norm_temp(dc_f_norm_idx_plus)-...
        dc_f_normalize_min)/abs(dc_f_normalize_max-dc_f_normalize_min);

    dc_f_normalize(dc_f_norm_idx_minus) = (dc_f_norm_temp(dc_f_norm_idx_minus)+...
        dc_f_normalize_min)/abs(dc_f_normalize_max-dc_f_normalize_min);

    x_criteria = x_1d(x_d);           
    
    x_max = max(abs(x_criteria));
    x_min = 0;

    x_norm = (x_criteria-x_min)/(x_max-x_min); % minimum density value for grouping criteria, '19.02.17
   
  
    dc_f_normalize2 = dc_f_normalize;
    dc_f_normalize2(dc_f_norm_idx_plus)=dc_f_normalize(dc_f_norm_idx_plus).^(1/penal); % normalized dc_f map, '19.12.14
    dc_f_normalize2(dc_f_norm_idx_minus)=-1*(abs(dc_f_normalize(dc_f_norm_idx_minus)).^(1/penal)); % normalized dc_f map, '19.12.14
    
    mod_dc_f=dc_f_normalize2.*x_norm; % normalized dc_f map, '19.12.14
            
end

%%

time_ansys(5,1)=time_ansys(5,1)+toc-loop_s_time;

loop_s_time=toc;


if data_N==0
    
else
    sens_map_mat=sub_sens_map_r4_5(mod_dc_f(:),data_N2);
end


%%
time_ansys(6,1)=time_ansys(6,1)+toc-loop_s_time;

loop_s_time=toc;
if data_N==0
    df0dx_temp=dc_f(:);
else
    df0dx_temp=sub_sens_reduce_r1(dc_f(:),sens_map_mat,data_N2);    
end


df0dx_w=1/c_ini; % dc normalization weight, initial compliance;
df0dx=df0dx_temp*df0dx_w; % dc normalize
clear df0dx_temp;

f0val = c*df0dx_w; % f0val nomalize


time_ansys(7,1)=time_ansys(7,1)+toc-loop_s_time;

%% histogram plot (IVG)
if logon_id==1 && data_N>0
      
    x_hist=1:1:data_N2;
    figure(21);hist(sens_map_mat,x_hist);grid on;title(['Histogram for the number of elements in the group (',num2str(loop),' iter.)']);xlabel('group variable number');ylabel('number of design variable');

    f = getframe(gcf);
    imwrite(f.cdata,[save_path,'histgram_',num2str(loop),'.jpg']);
       
end
%%

loop_s_time=toc;
if data_N==0
    xval=x(x_d);
    eff_x_d = find(xval>=0);    
else
    
    xval_temp=sub_x_reduce_r2(x(x_d),sens_map_mat,data_N2);
    xold1_temp=sub_x_reduce_r2(xold1_xd(x_d),sens_map_mat,data_N2);
    xold2_temp=sub_x_reduce_r2(xold2_xd(x_d),sens_map_mat,data_N2);
    
    eff_x_d = find(xval_temp);

    xval(eff_x_d) = xval_temp(eff_x_d);
    xold1(eff_x_d) = xold1_temp(eff_x_d);
    xold2(eff_x_d) = xold2_temp(eff_x_d);    
    
end

   

time_ansys(8,1)=time_ansys(8,1)+toc-loop_s_time;

loop_s_time=toc;


[con,coneq,gcon,gconeq,volfrac_eff]=cal_vol_const7(volfrac,nelx_d,nely_d,x_phys_FE);


time_ansys(9,1)=time_ansys(9,1)+toc-loop_s_time;

loop_s_time=toc;


gcon_1d=gcon(:);

if filter_id==1
    gcon_f=gcon_1d(x_d);
    %     gcon_f = H*(x_1d(x_d).*gcon_1d(x_d))./Hs./max(1e-3,x_1d(x_d));
elseif filter_id==2
    gcon_f = H*(gcon_1d(x_d)./Hs);
elseif filter_id==3
    if beta == 0
        gcon_f = H*(gcon_1d(x_d)./Hs);
    else
        gcon_f = H*(gcon_1d(x_d).*dx(:)./Hs);
    end
elseif filter_id==4
    gcon_f = H*(gcon_1d(x_d).*dx(:)./Hs);
elseif filter_id==5
    gcon_f = gcon_1d(x_d);
end



if data_N==0
    dfdx_temp=gcon_f;
else
    dfdx_temp=sub_sens_reduce_r1(gcon_f,sens_map_mat,data_N2);
    
end

dfdx_w=1/nele; % normalization, initial volume
dfdx=dfdx_temp'*dfdx_w; % gcon normalize
clear dfdx_temp;
clear gcon;
clear gcon_f;

fval  = con*dfdx_w; % con normalize


time_ansys(10,1)=time_ansys(10,1)+toc-loop_s_time;

loop_s_time=toc;

% Penalization increment, MMA low/upp initialize
penal_inc_cnt_old = penal_inc_cnt;
if penal_inc_cnt > 0
            
    if penal_inc_cnt >= 2
        penal_inc_cnt = 0;
    else
        penal_inc_cnt = penal_inc_cnt+1;
    end
   
end

% Projection filter beta increment, MMA low/upp initialize
beta_inc_cnt_old = beta_inc_cnt;
if (filter_id==3 || filter_id==4) && beta_inc_cnt > 0 
    
        
    if beta_inc_cnt >= 2
        beta_inc_cnt = 0;
    else
        beta_inc_cnt = beta_inc_cnt+1;
    end
   
end


xmma = xval;

[xmma_temp,ymma,zmma,lam,xsi,eta_mma,mu,zet,s,low_temp,upp_temp] = mmasub(m, size(eff_x_d,1), loop_MMA,...
    xval(eff_x_d), xmin(eff_x_d), xmax(eff_x_d), xold1(eff_x_d), xold2(eff_x_d),...
    f0val,df0dx(eff_x_d),fval,dfdx(:,eff_x_d),low(eff_x_d),upp(eff_x_d),a0,a,c_MMA,d);


if loop_MMA >=1
    [residu,kktnorm,residumax] = ...
        kktcheck(m,size(eff_x_d,1),xmma_temp,ymma,zmma,lam,xsi,eta_mma,mu,zet,s, ...
        xmin(eff_x_d), xmax(eff_x_d),df0dx(eff_x_d),fval,dfdx(:,eff_x_d),a0,a,c_MMA,d);
    
end


xmma(eff_x_d) = xmma_temp;
low(eff_x_d) = low_temp;
upp(eff_x_d) = upp_temp;

% plot xval, xold1, xold2, xmma, for debug
if logon_id==1 && data_N>0
   
    %x_hist=1:1:data_N2;
    figure(41);
    plot(xmma(eff_x_d),'r');hold on;
    plot(xval(eff_x_d),'b');
    plot(xold1(eff_x_d),'g');
    plot(xold2(eff_x_d),'k');hold off;
    
    grid on;title(['Group variable value (',num2str(loop),' iter.)']);xlabel('level(eff.)');ylabel('value');
    legend('xmma','xval','xold1','xold2')
    f = getframe(gcf);
    imwrite(f.cdata,[save_path,'xmma_xval_xold_value_',num2str(loop),'.jpg']);
    
end


if data_N > 0 
    if data_N == 1
        
        xmma_max = max(xmma_temp);
        xmma_min = min(xmma_temp);
        xmma_band = xmma_max - xmma_min;
        
        
        group_num_estimation_band = intersect(find(xmma_temp<(xmma_min+xmma_band*0.90)),...
            find(xmma_temp>(xmma_min+xmma_band*0.10)));
        
        group_num_estimation = size(group_num_estimation_band,1);

      disp(['Proper group num. = ',num2str(group_num_estimation),...          
          '/ data_N_tg (',num2str(data_N_tg),') / Eff.group_num (',num2str(size(eff_x_d,1)),')'])
      
    else
        group_num_estimation = size(eff_x_d,1);
        
    end
        
end


%%

time_ansys(11,1)=time_ansys(11,1)+toc-loop_s_time;



loop_s_time=toc;

if data_N==0
    xnew_xd = xmma;
    xold2 = xold1;
    xold1 = xval;
else
    xnew_xd=sub_x_recover_r1(x(x_d),sens_map_mat,xmma);
    xold2_xd(x_d)=sub_x_recover_r1(x(x_d),sens_map_mat,xold1);
    xold1_xd(x_d)=sub_x_recover_r1(x(x_d),sens_map_mat,xval);    
      
end


time_ansys(12,1)=time_ansys(12,1)+toc-loop_s_time;

%% filter

if filter_id==1
    x_phys_temp = xnew_xd;
    
elseif filter_id==2 % density filter
    x_phys_temp = (H*xnew_xd)./Hs; % density filter
    
elseif filter_id==3  % density projection filter
    
    xTilde = (H*xnew_xd)./Hs; % density filter with heaviside step function
    if beta == 0
        x_phys_temp = xTilde;        
    else
        x_phys_temp = (tanh(beta*eta)+tanh(beta*(xTilde-eta)))/(tanh(beta*eta)+tanh(beta*(1-eta))); % threshold projection filter    
    end
    
elseif filter_id==4  % density projection filter
    xTilde = (H*xnew_xd)./Hs; % density filter with heaviside step function
    x_phys_temp = 1-exp(-beta*xTilde)+xTilde*exp(-beta); % density filter with heaviside step function
        
elseif filter_id==5      
    x_phys_temp = xnew_xd;
    
end

xnew_temp=x(:);
xnew_temp(x_d)=xnew_xd;
xnew_temp(x_nd_s)=1;
xnew_temp(x_nd_v)=0;

xnew     = reshape(xnew_temp,nely_d,nelx_d);
x_phys_d = x_phys_temp;

change = max(abs(xval(eff_x_d)-xmma(eff_x_d)));

x = xnew;

change_obj = abs(c - c_old)/abs(c_old);
c_old = c;

x_phys=x(:);
x_phys(x_d)=x_phys_d;
x_phys(x_nd_s)=1;
x_phys(x_nd_v)=0;

x_phys = reshape(x_phys,nely_d,nelx_d);
debug_data(loop,4)=sum(x_phys(:))/(nelx_d*nely_d);

% debugging data
debug_data(loop,1)=loop; % Iteration number
debug_data(loop,2)=c; % Compliance
debug_data(loop,3)=fval(1); % objective functon value
debug_data(loop,5)=penal; % penalization factor
debug_data(loop,6)=Mnd2; % gray level indicator
debug_data(loop,7)=beta; % projection filter beta
debug_data(loop,8)=change; % max design variable change value
debug_data(loop,9)=size(eff_x_d,1); % Effective design variable number for MMA solver
debug_data(loop,10)=kktnorm; % KKT norm for MMA solver


%% plot result

% PRINT RESULTS
  
disp([' It.: ' sprintf('%4i',loop) ' Obj.: ' sprintf('%10.4f',c) ...
    ' Vol.: ' sprintf('%6.3f',debug_data(loop,4)) ...
    ' ch.: ' sprintf('%4.3f',change ) ...
    ' Gray_indicator : ' sprintf('%6.3f',Mnd2) ...
    %'  max(U) @ load pt :' sprintf('%6.3f',u1*1000) ' mm'])
    ])

if logon_id==1
    
    save([save_path,'density_result_2d_',num2str(loop),'.mat'],'x_phys','nelx_d','nely_d','XL','YL','loop','Mnd2','beta','xold1','xold2','xTilde','KG','U','dc','dc_f','x','xval','xmma'); % save iteration result
    
elseif logon_id==2
    save([save_path,'density_result_2d_',num2str(loop),'.mat'],'x_phys','nelx_d','nely_d','XL','YL','loop','Mnd2','beta','U','xmma'); % save iteration result
    
end

% PLOT DENSITIES  
loop_3dplot_time=0;
if logon_id==1
    loop_s_time=toc;    
    
    plot_density(x_phys,nelx_d,nely_d,XL,YL,loop,Mnd2)
    
    f = getframe(gcf);
    imwrite(f.cdata,[save_path,'density_dist_result_2d_',num2str(loop),'.jpg']);
    
    loop_3dplot_time=toc-loop_s_time;
end


end

if kktnorm <= kkt_tol % '19.03.10, KKTnorm tolerance
    exitflag=0;
    
elseif loop >= maxloop 
    exitflag=1;
    
elseif change_obj <= tol_obj
    exitflag=2;
end


if logon_id==0 || logon_id==2
       
    plot_density(x_phys,nelx_d,nely_d,XL,YL,loop,Mnd2)
    
    f = getframe(gcf);
    imwrite(f.cdata,[save_path,'density_dist_result_2d_',num2str(loop),'.jpg']);
    save([save_path,'density_result_2d_',num2str(loop),'.mat'],'x_phys','nelx_d','nely_d','XL','YL','loop','Mnd2','beta','xold1','xold2','xTilde','KG','U','xmma','xval'); % save iteration result
    
end


toc
t_total=t_total_ini+toc-loop_3dplot_time*loop; % total time (without 3D plot time)

time_ansys=time_ansys/loop;

debug_data = debug_data(1:loop,:);
result = struct('It',       loop, ...
                   'Obj',  c, ...
                   'Vol',  debug_data(loop,4), ...
                   'X_change',  change, ...
                   't_total',t_total, ...
                   'avg_eff_x_size', mean(debug_data(:,9)), ...
                   'exitflag',exitflag)
              
               

% plot debug data

debug_data_N=max(debug_data(:,1));

figure(2)
plot(debug_data(1:debug_data_N,1),debug_data(1:debug_data_N,2));
title('Objective function')
ylabel('value')
xlabel('Iteration')
grid on

f = getframe(gcf);
imwrite(f.cdata,[save_path,'Obj_function_plot.jpg']);



figure(3)
plot(debug_data(1:debug_data_N,1),debug_data(1:debug_data_N,3));
title('Constraint function')
ylabel('value')
xlabel('Iteration')
grid on

f = getframe(gcf);
imwrite(f.cdata,[save_path,'Const_function_plot.jpg']);


figure(4)
plot(debug_data(1:debug_data_N,1),debug_data(1:debug_data_N,4));
title('Volume fraction')
ylabel('value')
xlabel('Iteration')
grid on

f = getframe(gcf);
imwrite(f.cdata,[save_path,'Volume_frac_plot.jpg']);


if data_N > 0 && (logon_id == 0 || logon_id == 2)
    figure(5);
    plot(xmma(eff_x_d),'r');hold on;
    plot(xval(eff_x_d),'b');
    plot(xold1(eff_x_d),'g');
    plot(xold2(eff_x_d),'k');hold off;
    
    grid on;title(['Group variable value (',num2str(loop),' iter.)']);xlabel('level(eff.)');ylabel('value');
    legend('xmma','xval','xold1','xold2')
    f = getframe(gcf);
    imwrite(f.cdata,[save_path,'xmma_xval_xold_value_',num2str(loop),'.jpg']);
end

% time analysis
disp('=== Elapsed time analysis ==========================')
disp(['#1 GK calc :',num2str(time_ansys(1,1)),' sec'])
disp(['#2 FE calc :',num2str(time_ansys(2,1)),' sec'])
disp(['#3 Compliance :',num2str(time_ansys(3,1)),' sec'])
disp(['#4 Filter :',num2str(time_ansys(4,1)),' sec'])
disp(['#5 Variable reduce criteria :',num2str(time_ansys(5,1)),' sec'])
disp(['#6 Variable reduce map :',num2str(time_ansys(6,1)),' sec'])
disp(['#7 Variable reduce sensitivity :',num2str(time_ansys(7,1)),' sec'])
disp(['#8 Variable reduce x :',num2str(time_ansys(8,1)),' sec'])
disp(['#9 Constraint function calc.:',num2str(time_ansys(9,1)),' sec'])
disp(['#10 Constraint derivative reduced calc:',num2str(time_ansys(10,1)),' sec'])
disp(['#11 MMA done :',num2str(time_ansys(11,1)),' sec'])
disp(['#12 Variable recover :',num2str(time_ansys(12,1)),' sec'])

time_ansys_reduce=time_ansys(5,1)+time_ansys(6,1)+time_ansys(7,1)+time_ansys(8,1)+time_ansys(10,1)+time_ansys(12,1);
disp(['---------------------------------------------'])
disp(['Variable reduce related calculations :',num2str(time_ansys_reduce),' sec'])


time_result = struct('time_ansys',time_ansys','time_ansys_reduce',time_ansys_reduce);

save([save_path,'debug_data.mat'],'debug_data','inp_var','result','time_result');


figure(6)
bar(time_ansys);
title('Elapsed time (during iteration)')
xlabel('module')
ylabel('time (sec)')
grid on

f = getframe(gcf);
imwrite(f.cdata,[save_path,'Time_ansys_plot.jpg']);


