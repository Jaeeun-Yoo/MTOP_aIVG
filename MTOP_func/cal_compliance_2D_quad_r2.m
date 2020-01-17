%OBJECTIVE FUNCTION AND SENSITIVITY ANALYSIS
function [f,gf]=cal_compliance_2D_quad_r2(x,nelx,nely,penal,U,E,E_min,DLTOP_ratio,I_mat)


N_r=DLTOP_ratio;
N_s=DLTOP_ratio;

nele=nelx*nely;
KE_size = 18;
sub_elem_tot = N_r*N_s;

edofMat = zeros(nele,KE_size);
x2 = zeros(nele,sub_elem_tot);

elem_num = 0;
for elx = 1:nelx
    for ely = 1:nely
        elem_num = elem_num + 1;
        
        n1 = (2*nely+1)*(2* elx-2) + 2*ely-1;
        n2 = (2*nely+1)*(2* elx)   + 2*ely-1;
        n3 = (2*nely+1)*(2* elx-1) + 2*ely-1;
        edof = [2*n1-1 2*n1 2*n2-1 2*n2 2*n2+3 2*n2+4 2*n1+3 2*n1+4 2*n3-1 2*n3 2*n2+1 2*n2+2 2*n3+3 2*n3+4 2*n1+1 2*n1+2 2*n3+1 2*n3+2];
        
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

ce = zeros(nele,sub_elem_tot);
for k=1:sub_elem_tot
    ce_temp = sum((U(edofMat)*I_mat2(:,:,k)).*U(edofMat),2);
    ce(:,k) = ce_temp;
end

c = sum((E_min+x2(:).^penal*(E-E_min)).*ce(:));
dc_temp = -penal*(E-E_min)*x2(:).^(penal-1).*ce(:);

dc = reshape(permute(reshape(dc_temp,nely,nelx,N_s,N_r),[3,1,4,2]),nely*N_s,nelx*N_r);


f = c;
gf = dc;




