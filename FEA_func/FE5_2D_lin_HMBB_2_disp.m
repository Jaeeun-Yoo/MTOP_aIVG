
%%%%%%%%%% FE-ANALYSIS %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function [U,F]=FE5_2D_lin_HMBB_2_disp(nelx,nely,KG,Fnew,u_old,solver_type)
% [K0] = lk3(elem_size);
% 
% K = sparse(2*(nelx+1)*(nely+1), 2*(nelx+1)*(nely+1));
F_size = 3; % Load set = 1, dummy variable 2 (for displacement)
F = zeros(2*(nely+1)*(nelx+1),F_size); 
U = zeros(2*(nely+1)*(nelx+1),F_size);
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
F(2*(nely+1),1)=Fnew; % x direction, load#1
% F(2*(nely+1)*round(nelx*16/30)+2*(nely+1),1)=Fnew; % y direction, load#1
F(2*(nely+1)*round(nelx*10/30)+2*(nely+1),1)=Fnew; % y direction, load#1

F(2*(nely+1),2)=1; % y direction,  dummy load, #1
% F(2*(nely+1)*round(nelx*16/30)+2*(nely+1),3)=1; % y direction,  dummy load, #2
F(2*(nely+1)*round(nelx*10/30)+2*(nely+1),3)=1; % y direction,  dummy load, #2

%F(2*(nelx+1)*(nely+1)-nely+1,1)=-0.1*Fnew;


%F(2*(nelx+1)*(nely+1),1)=-1.0;
fixeddofs   = union([1:2:2*(nely+1)],[2*(nelx+1)*(nely+1)]);
alldofs     = [1:2*(nely+1)*(nelx+1)];
freedofs    = setdiff(alldofs,fixeddofs);

%% SOLVING
[U_sparse]=FEA_solver_r2(KG,F,freedofs,u_old,solver_type);

U(freedofs,:)=full(U_sparse);
U(fixeddofs,:)= 0;

