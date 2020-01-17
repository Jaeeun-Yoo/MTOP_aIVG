function [m,n,xmin,xmax,xold1,xold2,low,upp,a0,a,c_MMA,d] = MMA_param_initialize(data_N3,volfrac,m)

% INITIALIZE MMA OPTIMIZER
% m     = 1;                % The number of general constraints.
n     = data_N3;             % The number of design variables x_j.
xmin  = zeros(data_N3,1);       % Column vector with the lower bounds for the variables x_j.
% xmin  = 10^(-3)*ones(data_N3,1);       % Column vector with the lower bounds for the variables x_j.
xmax  = ones(data_N3,1);        % Column vector with the upper bounds for the variables x_j.

xold1 = volfrac*ones(data_N3,1);             % xval, one iteration ago (provided that iter>1).
xold2 = volfrac*ones(data_N3,1);             % xval, two iterations ago (provided that iter>2).

low   = zeros(data_N3,1);        % Column vector with the lower asymptotes from the previous iteration (provided that iter>1).
upp   = ones(data_N3,1);        % Column vector with the upper asymptotes from the previous iteration (provided that iter>1).
a0    = 1;                % The constants a_0 in the term a_0*z.
a     = zeros(m,1);       % Column vector with the constants a_i in the terms a_i*z.
c_MMA = 1.0*10^5*ones(m,1);   % Column vector with the constants c_i in
% c_MMA = 1.0*10^4*ones(m,1);   % Column vector with the constants c_i in  % update '19.11.09
d     = zeros(m,1);       % Column vector with the constants d_i in the terms 0.5*d_i*(y_i)^2.

%---  Ph.D. kang's code (IGA top)
% 
%       n=size(xval,1);
%       
%       m=1; a0=1; c=10000*ones(m,1); d=zeros(m,1); a=zeros(m,1);
%       epsimin = 0.0000001;
%       
%       xold1=xval; xold2=xval;
%       
%       low   = xmin;
%       upp   = xmax;
