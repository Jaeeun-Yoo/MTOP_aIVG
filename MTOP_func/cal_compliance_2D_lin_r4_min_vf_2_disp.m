%OBJECTIVE FUNCTION AND SENSITIVITY ANALYSIS
function [c_disp,dc]=cal_compliance_2D_lin_r4_min_vf_2_disp(x,nelx,nely,penal,U,E,E_min,DLTOP_ratio,I_mat,F_vec)


N_r=DLTOP_ratio;
N_s=DLTOP_ratio;

nele=nelx*nely;
KE_size = 8;
sub_elem_tot = N_r*N_s;

edofMat = zeros(nele,KE_size);
x2 = zeros(nele,sub_elem_tot);
% size_m = size(U,2); % multiple load size
eff_m = 2; % effective load size except dummy load

elem_num = 0;
for elx = 1:nelx
    for ely = 1:nely
        elem_num = elem_num + 1;
        
        n1 = 2*(nely+1)*(elx-1)+2*ely-1;
        n2 = 2*(nely+1)* elx   +2*ely-1;
        
        edof = [n1 n1+1 n2 n2+1 n2+2 n2+3 n1+2 n1+3];
        edofMat(elem_num,:) = edof;
        
        sub_elem_num = 0;
        for i=1:N_r
            elx_d=DLTOP_ratio*(elx-1)+i;
            
            for j=1:N_s
                ely_d=DLTOP_ratio*(ely-1)+j;
                sub_elem_num=sub_elem_num+1;
                
                x2(elem_num,sub_elem_num) = x(ely_d,elx_d);
                
            end
        end
    end
end


I_mat2 = reshape(I_mat,KE_size,KE_size,sub_elem_tot);

ce1 = zeros(nele,sub_elem_tot);
ce2 = zeros(nele,sub_elem_tot);
for k=1:sub_elem_tot
    U_temp = U(:,1);
    Ud_temp = U(:,2);
    ce_temp = sum((U_temp(edofMat)*I_mat2(:,:,k)).*Ud_temp(edofMat),2);
    ce1(:,k) = ce_temp;
    
    U_temp = U(:,1);
    Ud_temp = U(:,3);
    ce_temp = sum((U_temp(edofMat)*I_mat2(:,:,k)).*Ud_temp(edofMat),2);
    ce2(:,k) = ce_temp;
end

% c_elem_temp = (E_min+x2(:).^penal*(E-E_min)).*ce(:);
% c_elem_temp_no_penal = (E_min+x2(:)*(E-E_min)).*ce(:);
% c = sum(c_elem_temp);
c_disp = zeros(eff_m,1);
dc = zeros(nely*N_s,nelx*N_r,eff_m);

for m=1:eff_m
    U_temp = U(:,1);
    c_disp(m,1) = F_vec(:,m+1)'*U_temp;
end

dc_temp1 = -penal*(E-E_min)*x2(:).^(penal-1).*ce1(:);
dc_temp2 = -penal*(E-E_min)*x2(:).^(penal-1).*ce2(:);
% dc_temp_no_penal = -1*(E-E_min).*ce(:);

dc(:,:,1) = reshape(permute(reshape(dc_temp1,nely,nelx,N_s,N_r),[3,1,4,2]),nely*N_s,nelx*N_r);
dc(:,:,2) = reshape(permute(reshape(dc_temp2,nely,nelx,N_s,N_r),[3,1,4,2]),nely*N_s,nelx*N_r);




