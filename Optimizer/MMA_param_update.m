function [m,n,xmin,xmax,xold1_new,xold2_new,low_new,upp_new,a0,a,c_MMA,d,xval_new] = MMA_param_update(data_N_old,data_N,xold1,xold2,low,upp,xval,m)

% INITIALIZE MMA OPTIMIZER
% m     = 1;                % The number of general constraints.
n     = data_N;             % The number of design variables x_j.
xmin  = zeros(data_N,1);       % Column vector with the lower bounds for the variables x_j.
% xmin  = 10^(-3)*ones(data_N,1);       % Column vector with the lower bounds for the variables x_j.
xmax  = ones(data_N,1);        % Column vector with the upper bounds for the variables x_j.

%xold1 = volfrac*ones(data_N3,1);             % xval, one iteration ago (provided that iter>1).
%xold2 = volfrac*ones(data_N3,1);             % xval, two iterations ago (provided that iter>2).

%low   = zeros(data_N3,1);        % Column vector with the lower asymptotes from the previous iteration (provided that iter>1).
%upp   = ones(data_N3,1);        % Column vector with the upper asymptotes from the previous iteration (provided that iter>1).
a0    = 1;                % The constants a_0 in the term a_0*z.
a     = zeros(m,1);       % Column vector with the constants a_i in the terms a_i*z.
c_MMA = 1.0*10^5*ones(m,1);   % Column vector with the constants c_i in
% c_MMA = 1.0*10^4*ones(m,1);   % Column vector with the constants c_i in, % update '19.11.09
d     = zeros(m,1);       % Column vector with the constants d_i in the terms 0.5*d_i*(y_i)^2.


%x_new = zeros(data_N,1);
xold1_new = zeros(data_N,1);
xold2_new = zeros(data_N,1);
low_new = zeros(data_N,1);
upp_new = zeros(data_N,1);
xval_new = zeros(data_N,1);

if data_N_old < data_N
    for j=1:data_N_old
        %x_new(2*j-1,1) = x(j,1);
        %x_new(max(data_N,2*j),1) = x(j,1);
        xval_new(2*j-1,1) = xval(j,1);
        xval_new(min(data_N,2*j),1) = xval(j,1);
        xold1_new(2*j-1,1) = xold1(j,1);
        xold1_new(min(data_N,2*j),1) = xold1(j,1);
        xold2_new(2*j-1,1) = xold2(j,1);
        xold2_new(min(data_N,2*j),1) = xold2(j,1);
        low_new(2*j-1,1) = low(j,1);
        low_new(min(data_N,2*j),1) = low(j,1);
        upp_new(2*j-1,1) = upp(j,1);
        upp_new(min(data_N,2*j),1) = upp(j,1);
    end
else
    for j=1:data_N
        %x_new(j,1) = 0.5*(x(2*j-1,1)+x(max(data_N,2*j),1));
        xval_new(j,1) = 0.5*(xval(2*j-1,1)+xval(max(data_N,2*j),1));
        xold1_new(j,1) = 0.5*(xold1(2*j-1,1)+xold1(max(data_N,2*j),1));
        xold2_new(j,1) = 0.5*(xold2(2*j-1,1)+xold2(max(data_N,2*j),1));
        low_new(j,1) = 0.5*(low(2*j-1,1)+low(max(data_N,2*j),1));
        upp_new(j,1) = 0.5*(upp(2*j-1,1)+upp(max(data_N,2*j),1));
        
    end
end

