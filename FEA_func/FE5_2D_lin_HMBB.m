
%%%%%%%%%% FE-ANALYSIS %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function [U]=FE5_2D_lin_HMBB(nelx,nely,KG,Fnew,theta,u_old,solver_type)
% [K0] = lk3(elem_size);
% 
% K = sparse(2*(nelx+1)*(nely+1), 2*(nelx+1)*(nely+1));
F = zeros(2*(nely+1)*(nelx+1),1); 
U = zeros(2*(nely+1)*(nelx+1),1);
% for elx = 1:nelx
%   for ely = 1:nely
%     n1 = (nely+1)*(elx-1)+ely; 
%     n2 = (nely+1)* elx   +ely;
%     edof = [2*n1-1; 2*n1; 2*n2-1; 2*n2; 2*n2+1; 2*n2+2; 2*n1+1; 2*n1+2];
%     K(edof,edof) = K(edof,edof) + (E_min+x(ely,elx)^penal*(E-E_min))*K0;
%   end
% end
% DEFINE LOADS AND SUPPORTS (HALF MBB-BEAM)
% F(2,1) = -1;
% fixeddofs   = union([1:2:2*(nely+1)],[2*(nelx+1)*(nely+1)]);
% alldofs     = [1:2*(nely+1)*(nelx+1)];
% freedofs    = setdiff(alldofs,fixeddofs);
% DEFINE LOADS AND SUPPORTS (Cantilever)

%theta=0; % 0 deg = x direction, + deg = upward
F(1,1)=Fnew*sin(theta*pi()/180); % x direction
F(2,1)=Fnew*cos(theta*pi()/180); % y direction
%F(2*(nelx+1)*(nely+1)-nely+1,1)=-0.1*Fnew;


%F(2*(nelx+1)*(nely+1),1)=-1.0;
fixeddofs   = union([1:2:2*(nely+1)],[2*(nelx+1)*(nely+1)]);
alldofs     = [1:2*(nely+1)*(nelx+1)];
freedofs    = setdiff(alldofs,fixeddofs);

%% SOLVING
[U_sparse]=FEA_solver_r1(KG,F,freedofs,u_old,solver_type);

U(freedofs,:)=full(U_sparse);
U(fixeddofs,:)= 0;

