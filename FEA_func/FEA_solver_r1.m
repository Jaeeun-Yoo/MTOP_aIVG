
%%%%%%%%%% FE-solving %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% solver_type :
%  'Direct' : Matlab Direct method (sparse matrix solver)
%  'PCG' : Preconditioned Conjugate Gradient iterative solver
%
% update : '19.01.19

function [U_sparse]=FEA_solver_r1(KG,F,freedofs,u_old,solver_type)

%% SOLVING

switch solver_type
    case 'Direct'
        %% Direct solver
        
        a=toc;
        U_sparse = KG(freedofs,freedofs) \ F(freedofs,1);
        b=toc;
        disp(['Direct solver time = ',num2str(b-a)])
                
        
    case 'PCG'
        %% iterative solver, PCG
        
        ts=toc;
        tolit = 1e-5;
        maxit = 2000;
        
        alpha = 0.01; % '19.12.14 
        L = ichol(KG(freedofs,freedofs), struct('type','ict','droptol',1e-3,'diagcomp',alpha));
        
        [U_sparse,flag,relres,iter,rv]  = pcg(KG(freedofs,freedofs),F(freedofs,1),tolit,maxit,L,L',u_old(freedofs,1));

        te=toc;

        disp(['iterative solver(PCG, ichol diagcomp ict1e-3) time = ',num2str(te-ts),'(iter=',num2str(iter),') / alpha = ',num2str(alpha)])
        
end
