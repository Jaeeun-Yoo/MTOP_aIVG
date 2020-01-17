
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
        
        % F2=F;
        % freedofs2=freedofs;
        % U2=U_sparse;
        % save('KG_debug.mat','KG','freedofs2','F2','U2')
        
        
    case 'PCG'
        %% iterative solver, PCG
        
        ts=toc;
        tolit = 1e-5;
        maxit = 2000;
        
%         alpha = 0.002*(max(sum(abs(KG(freedofs,freedofs)),2)./diag(KG(freedofs,freedofs)))-2); % '19.03.30
%         alpha = 0.1; % '19.03.19
        alpha = 0.01; % '19.12.14 
        L = ichol(KG(freedofs,freedofs), struct('type','ict','droptol',1e-3,'diagcomp',alpha));
%         L = ichol(KG(freedofs,freedofs), struct('type','nofill','michol','on'));
%         L = ichol(KG(freedofs,freedofs), struct('michol','on'));
%         L = ichol(KG(freedofs,freedofs));
        
        [U_sparse,flag,relres,iter,rv]  = pcg(KG(freedofs,freedofs),F(freedofs,1),tolit,maxit,L,L',u_old(freedofs,1));
        %[U_sparse,flag,relres,iter,rv]  = pcg(KG(freedofs,freedofs),F(freedofs,1),tolit,maxit,L,L');
        te=toc;
        %disp(['iterative solver(PCG) time = ',num2str(b-a)])
        disp(['iterative solver(PCG, ichol diagcomp ict1e-3) time = ',num2str(te-ts),'(iter=',num2str(iter),') / alpha = ',num2str(alpha)])
%         disp(['iterative solver(PCG, ichol diagcomp ict1e-3) time = ',num2str(te-ts),'(iter=',num2str(iter),') / Modified ichol'])
        %save('KG_debug.mat','KG','freedofs','F','U_sparse','L','alpha','flag','relres','iter','tolit')
        
end
