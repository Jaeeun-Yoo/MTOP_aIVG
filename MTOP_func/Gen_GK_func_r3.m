
%%%%%%%%%% ELEMENT STIFFNESS MATRIX %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% DLSIMP local stiffness matrix 
function [KG]=Gen_GK_func_r3(nelx,nely,nelz,x,penal,E,E_min,DLTOP_ratio,I_mat,elem_type)
%E = 210*10^9; % Elastic modulus, Pa, (210 GPa) 
%E_min= 10^(-3); % elastic modulus minimum, Pa,

% nelx_d=nelx*DLTOP_ratio;
% nely_d=nely*DLTOP_ratio;

switch (elem_type) 
    case 'QUAD4' 
        
        N_r=DLTOP_ratio;
        N_s=DLTOP_ratio;
                 
        nele=nelx*nely;
                
        KE_size = 8;      
        edofMat = zeros(nele,KE_size);
        x2 = zeros(nele,N_r*N_s);
        
        elem_num=0;
        for elx = 1:nelx
            for ely = 1:nely
                elem_num=elem_num+1;
                
                n1 = 2*(nely+1)*(elx-1)+2*ely-1;
                n2 = 2*(nely+1)* elx   +2*ely-1;
                
                edof = [n1 n1+1 n2 n2+1 n2+2 n2+3 n1+2 n1+3];
                edofMat(elem_num,:) = edof;                 
                
                sub_elem_num=0;
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
        
        I_mat2 = reshape(I_mat,KE_size*KE_size,N_r*N_s);
        sK = reshape(I_mat2*(E_min+x2'.^penal*(E-E_min)),KE_size*KE_size*nele,1);
        
        iK = reshape(kron(edofMat,ones(KE_size,1))',KE_size*KE_size*nele,1);
        jK = reshape(kron(edofMat,ones(1,KE_size))',KE_size*KE_size*nele,1);
        
        KG = sparse(iK,jK,sK);
        KG = (KG+KG')/2;
        
        
        
    case 'QUAD4_HT'
        N_r=DLTOP_ratio;
        N_s=DLTOP_ratio;
                 
        nele=nelx*nely;
        KE_size = 4;      
        edofMat = zeros(nele,KE_size);
        x2 = zeros(nele,N_r*N_s);
        
        elem_num=0;
        for elx = 1:nelx
            for ely = 1:nely
                elem_num=elem_num+1;
                
                n1 = (nely+1)*(elx-1)+ely;
                n2 = (nely+1)* elx   +ely;
                
                edof = [n1 n2 n2+1 n1+1];
                edofMat(elem_num,:) = edof; 
                
                sub_elem_num=0;                
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
        
        I_mat2 = reshape(I_mat,KE_size*KE_size,N_r*N_s);
        sK = reshape(I_mat2*(E_min+x2'.^penal*(E-E_min)),KE_size*KE_size*nele,1);
        
        iK = reshape(kron(edofMat,ones(KE_size,1))',KE_size*KE_size*nele,1);
        jK = reshape(kron(edofMat,ones(1,KE_size))',KE_size*KE_size*nele,1);
        
        KG = sparse(iK,jK,sK);
        KG = (KG+KG')/2;
        
        
    case 'QUAD8'   
                
        N_r=DLTOP_ratio;
        N_s=DLTOP_ratio;
                 
        nele=nelx*nely;
                
        KE_size = 18;      
        edofMat = zeros(nele,KE_size);
        x2 = zeros(nele,N_r*N_s);
        
        elem_num=0;
        for elx = 1:nelx
            for ely = 1:nely
                elem_num=elem_num+1;
                
                n1 = (2*nely+1)*(2* elx-2) + 2*ely-1;
                n2 = (2*nely+1)*(2* elx)   + 2*ely-1;
                n3 = (2*nely+1)*(2* elx-1) + 2*ely-1;
                edof = [2*n1-1 2*n1 2*n2-1 2*n2 2*n2+3 2*n2+4 2*n1+3 2*n1+4 2*n3-1 2*n3 2*n2+1 2*n2+2 2*n3+3 2*n3+4 2*n1+1 2*n1+2 2*n3+1 2*n3+2];
                
                edofMat(elem_num,:) = edof; 
                
                sub_elem_num=0;
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
        
        I_mat2 = reshape(I_mat,KE_size*KE_size,N_r*N_s);
        sK = reshape(I_mat2*(E_min+x2'.^penal*(E-E_min)),KE_size*KE_size*nele,1);
        
        iK = reshape(kron(edofMat,ones(KE_size,1))',KE_size*KE_size*nele,1);
        jK = reshape(kron(edofMat,ones(1,KE_size))',KE_size*KE_size*nele,1);
        
        KG = sparse(iK,jK,sK);
        KG = (KG+KG')/2;        
        
        
    case 'BRICK8'

        N_r=DLTOP_ratio;
        N_s=DLTOP_ratio;
        N_q=DLTOP_ratio;
                 
        nele=nelx*nely*nelz;
                
        KE_size = 24;      
        edofMat = zeros(nele,KE_size);
        x2 = zeros(nele,N_r*N_s*N_q);
        
        elem_num=0;
        for elz = 1:nelz
            for elx = 1:nelx
                for ely = 1:nely
                    
                    elem_num=elem_num+1;
                    
                    node_z_temp_1=(nely+1)*(nelx+1)*(elz-1);
                    node_z_temp_2=(nely+1)*(nelx+1)*(elz);
                    
                    node_x_temp_1=(nely+1)*(elx-1);
                    node_x_temp_2=(nely+1)* elx;
                    
                    n1 = node_z_temp_1+node_x_temp_1+ely;
                    n2 = node_z_temp_1+node_x_temp_2   +ely;
                    n3 = node_z_temp_2+node_x_temp_1 +ely;
                    n4 = node_z_temp_2+node_x_temp_2   +ely;
                    edof = [3*n1-2 3*n1-1 3*n1 3*n2-2 3*n2-1 3*n2 3*n2+1 3*n2+2 3*n2+3 3*n1+1 3*n1+2 3*n1+3 3*n3-2 3*n3-1 3*n3 3*n4-2 3*n4-1 3*n4 3*n4+1 3*n4+2 3*n4+3 3*n3+1 3*n3+2 3*n3+3];
                    
                    edofMat(elem_num,:) = edof;
                    
                    sub_elem_num=0;
                    for k=1:N_q
                        elz_d=DLTOP_ratio*(elz-1)+k;
                        for i=1:N_r
                            elx_d=DLTOP_ratio*(elx-1)+i;
                            for j=1:N_s
                                ely_d=DLTOP_ratio*(ely-1)+j;
                                
                                sub_elem_num=sub_elem_num+1;
                                
                                x2(elem_num,sub_elem_num) = x(ely_d,elx_d,elz_d);
                            end
                            
                        end
                    end
                end
            end
        end
        
        I_mat2 = reshape(I_mat,KE_size*KE_size,N_r*N_s*N_q);
        sK = reshape(I_mat2*(E_min+x2'.^penal*(E-E_min)),KE_size*KE_size*nele,1);
        
        iK = reshape(kron(edofMat,ones(KE_size,1))',KE_size*KE_size*nele,1);
        jK = reshape(kron(edofMat,ones(1,KE_size))',KE_size*KE_size*nele,1);
        
        KG = sparse(iK,jK,sK);
        KG = (KG+KG')/2;
        
        
    case 'BRICK8_HT'
        
        N_r=DLTOP_ratio;
        N_s=DLTOP_ratio;
        N_q=DLTOP_ratio;
                 
        nele=nelx*nely*nelz;
                
        KE_size = 8;      
        edofMat = zeros(nele,KE_size);
        x2 = zeros(nele,N_r*N_s*N_q);
        
        elem_num=0;
        for elz = 1:nelz
            for elx = 1:nelx
                for ely = 1:nely
                    
                    elem_num=elem_num+1;
                    
                    node_z_temp_1=(nely+1)*(nelx+1)*(elz-1);
                    node_z_temp_2=(nely+1)*(nelx+1)*elz;
                    
                    node_x_temp_1=(nely+1)*(elx-1);
                    node_x_temp_2=(nely+1)* elx;
                    
                    n1 = node_z_temp_1+node_x_temp_1+ely;
                    n2 = node_z_temp_1+node_x_temp_2   +ely;
                    n3 = node_z_temp_2+node_x_temp_1 +ely;
                    n4 = node_z_temp_2+node_x_temp_2   +ely;
                    edof = [n1 n2 n2+1 n1+1 n3 n4 n4+1 n3+1];
                    
                    edofMat(elem_num,:) = edof;
                    
                    sub_elem_num=0;
                    for k=1:N_q
                        elz_d=DLTOP_ratio*(elz-1)+k;
                        for i=1:N_r
                            elx_d=DLTOP_ratio*(elx-1)+i;
                            for j=1:N_s
                                ely_d=DLTOP_ratio*(ely-1)+j;
                                
                                sub_elem_num=sub_elem_num+1;
                                
                                x2(elem_num,sub_elem_num) = x(ely_d,elx_d,elz_d);
                            end
                            
                        end
                    end
                end
            end
        end
        
        I_mat2 = reshape(I_mat,KE_size*KE_size,N_r*N_s*N_q);
        sK = reshape(I_mat2*(E_min+x2'.^penal*(E-E_min)),KE_size*KE_size*nele,1);
        
        iK = reshape(kron(edofMat,ones(KE_size,1))',KE_size*KE_size*nele,1);
        jK = reshape(kron(edofMat,ones(1,KE_size))',KE_size*KE_size*nele,1);
        
        KG = sparse(iK,jK,sK);
        KG = (KG+KG')/2;
        
        
    case 'BRICK27'
        
        N_r=DLTOP_ratio;
        N_s=DLTOP_ratio;
        N_q=DLTOP_ratio;
                 
        nele=nelx*nely*nelz;
                
        KE_size = 81;      
        edofMat = zeros(nele,KE_size);
        x2 = zeros(nele,N_r*N_s*N_q);
        
        elem_num=0;
        for elz = 1:nelz
            for elx = 1:nelx
                for ely = 1:nely
                    
                    elem_num=elem_num+1;
                    
                    node_z_temp_1=(2*nely+1)*(2*nelx+1)*(2*elz-2);
                    node_z_temp_2=(2*nely+1)*(2*nelx+1)*(2*elz);
                    node_z_temp_3=(2*nely+1)*(2*nelx+1)*(2*elz-1);
                    
                    node_x_temp_1=(2*nely+1)*(2*elx-2);
                    node_x_temp_2=(2*nely+1)*(2*elx);
                    node_x_temp_3=(2*nely+1)*(2*elx-1);
                    
                    n1 = node_z_temp_1+node_x_temp_1+2*ely-1;
                    n2 = node_z_temp_1+node_x_temp_2+2*ely-1;
                    n3 = node_z_temp_1+node_x_temp_3+2*ely-1;
                    
                    n4 = node_z_temp_2+node_x_temp_1+2*ely-1;
                    n5 = node_z_temp_2+node_x_temp_2+2*ely-1;
                    n6 = node_z_temp_2+node_x_temp_3+2*ely-1;
                    
                    n7 = node_z_temp_3+node_x_temp_1+2*ely-1;
                    n8 = node_z_temp_3+node_x_temp_2+2*ely-1;
                    n9 = node_z_temp_3+node_x_temp_3+2*ely-1;
                    
                    edof = [3*n1-2; 3*n1-1; 3*n1; 3*n2-2; 3*n2-1; 3*n2; 3*n2+4; 3*n2+5; 3*n2+6; 3*n1+4; 3*n1+5; 3*n1+6 ;...
                        3*n4-2; 3*n4-1; 3*n4; 3*n5-2; 3*n5-1; 3*n5; 3*n5+4; 3*n5+5; 3*n5+6; 3*n4+4; 3*n4+5; 3*n4+6 ;...
                        3*n3-2; 3*n3-1; 3*n3; 3*n2+1; 3*n2+2; 3*n2+3; 3*n3+4; 3*n3+5; 3*n3+6; 3*n1+1; 3*n1+2; 3*n1+3 ;...
                        3*n7-2; 3*n7-1; 3*n7; 3*n8-2; 3*n8-1; 3*n8; 3*n8+4; 3*n8+5; 3*n8+6; 3*n7+4; 3*n7+5; 3*n7+6 ;...
                        3*n6-2; 3*n6-1; 3*n6; 3*n5+1; 3*n5+2; 3*n5+3; 3*n6+4; 3*n6+5; 3*n6+6; 3*n4+1; 3*n4+2; 3*n4+3 ;...
                        3*n3+1; 3*n3+2; 3*n3+3; 3*n9-2; 3*n9-1; 3*n9; 3*n8+1; 3*n8+2; 3*n8+3; 3*n9+4; 3*n9+5; 3*n9+6 ;3*n7+1; 3*n7+2; 3*n7+3 ;...
                        3*n6+1; 3*n6+2; 3*n6+3 ; 3*n9+1; 3*n9+2; 3*n9+3];
                    edofMat(elem_num,:) = edof;
                    
                    sub_elem_num=0;
                    for k=1:N_q
                        elz_d=DLTOP_ratio*(elz-1)+k;
                        for i=1:N_r
                            elx_d=DLTOP_ratio*(elx-1)+i;
                            for j=1:N_s
                                ely_d=DLTOP_ratio*(ely-1)+j;
                                
                                sub_elem_num=sub_elem_num+1;
                                
                                x2(elem_num,sub_elem_num) = x(ely_d,elx_d,elz_d);
                            end
                            
                        end
                    end
                end
            end
        end
        
        I_mat2 = reshape(I_mat,KE_size*KE_size,N_r*N_s*N_q);
        sK = reshape(I_mat2*(E_min+x2'.^penal*(E-E_min)),KE_size*KE_size*nele,1);
        
        iK = reshape(kron(edofMat,ones(KE_size,1))',KE_size*KE_size*nele,1);
        jK = reshape(kron(edofMat,ones(1,KE_size))',KE_size*KE_size*nele,1);
        
        KG = sparse(iK,jK,sK);
        KG = (KG+KG')/2;
        
        
end

