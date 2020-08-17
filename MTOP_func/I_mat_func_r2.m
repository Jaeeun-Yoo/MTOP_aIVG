
%%%%%%%%%% ELEMENT STIFFNESS MATRIX %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% DLSIMP local stiffness matrix 
function I_mat=I_mat_func_r2(elem_size,DLTOP_ratio,int_scale,nu,elem_type)

% a=elem_size/2;
% b=elem_size/2;
% c=elem_size/2;

a=elem_size/2;
b=elem_size/2;
c=elem_size/2;
        

switch (elem_type) 
    case 'QUAD4'
        
%         Ak = a*b/(DLTOP_ratio)^2/(DLTOP_ratio)^2;
%         Ak = (a*b)/4/(DLTOP_ratio)^2;
%         Ak = 1/(a*b);
        [wri,ri]=gausspt(int_scale);
        
        %nu = 0.3;
        %int_scale=1;
        
        D0=1/(1-nu^2)*[1 nu 0
            nu 1 0
            0 0 (1-nu)/2];
        
        % ri=[-1*sqrt(3)/3 +1*sqrt(3)/3]';
        % si=[-1*sqrt(3)/3 +1*sqrt(3)/3]';
        
        % N_r=size(ri,1);
        % N_s=size(si,1);
        
        N_r=DLTOP_ratio;
        N_s=DLTOP_ratio;
        %
        % ri=[-1+2/(N_r)/2:2/(N_r):1]';
        % si=[-1+2/(N_s)/2:2/(N_s):1]';
        
        ri_int=transpose(-1:2/N_r:1);
        si_int=transpose(-1:2/N_s:1);
        %qi_int=transpose(-1:2/N_q:1);
        
        
        % wri=2/N_r/int_scale;
        % wsi=2/N_s/int_scale;
        
        I_mat=zeros(8,8,N_r,N_s);
        %I_mat_temp=zeros(8,8,N_r,N_s);
        for j=1:N_s
            for i=1:N_r
                
                r_a=ri_int(i,1);
                r_b=ri_int(i+1,1);
                s_a=si_int(j,1);
                s_b=si_int(j+1,1);
                
                for j2=1:int_scale
                    si2=(s_b-s_a)/2*ri(j2,1)+(s_b+s_a)/2;
                    wsi2=(s_b-s_a)/2*wri(j2,1);
                    for i2=1:int_scale
                        ri2=(r_b-r_a)/2*ri(i2,1)+(r_b+r_a)/2;
                        wri2=(r_b-r_a)/2*wri(i2,1);
                        %I_mat=B_mat_func_r1(a,b,ri(j,1),si(i,1)); % 3x8 matrix result
                        I_mat(:,:,j,i)=I_mat(:,:,j,i)+wri2*wsi2*B_mat_func_r1(a,b,0,ri2,si2,0,elem_type)'*D0*B_mat_func_r1(a,b,0,ri2,si2,0,elem_type);
                    end
                end
                
                I_mat(:,:,j,i) = a*b*I_mat(:,:,j,i); 
                
            end
        end
%         I_mat = I_mat*Ak;
%         disp('test')
        
    case 'QUAD4_HT'
        
        [wri,ri]=gausspt(int_scale);
        
        %nu = 0.3;
        %int_scale=1;
        
        D0=[1 0
            0 1];
        
        % ri=[-1*sqrt(3)/3 +1*sqrt(3)/3]';
        % si=[-1*sqrt(3)/3 +1*sqrt(3)/3]';
        
        % N_r=size(ri,1);
        % N_s=size(si,1);
        
        N_r=DLTOP_ratio;
        N_s=DLTOP_ratio;
        %
        % ri=[-1+2/(N_r)/2:2/(N_r):1]';
        % si=[-1+2/(N_s)/2:2/(N_s):1]';
        
        ri_int=transpose(-1:2/N_r:1);
        si_int=transpose(-1:2/N_s:1);
        %qi_int=transpose(-1:2/N_q:1);
        
        
        % wri=2/N_r/int_scale;
        % wsi=2/N_s/int_scale;
        
        I_mat=zeros(4,4,N_r,N_s);
        %I_mat_temp=zeros(8,8,N_r,N_s);
        for j=1:N_s
            for i=1:N_r
                
                r_a=ri_int(i,1);
                r_b=ri_int(i+1,1);
                s_a=si_int(j,1);
                s_b=si_int(j+1,1);
                
                for j2=1:int_scale
                    si2=(s_b-s_a)/2*ri(j2,1)+(s_b+s_a)/2;
                    wsi2=(s_b-s_a)/2*wri(j2,1);
                    for i2=1:int_scale
                        ri2=(r_b-r_a)/2*ri(i2,1)+(r_b+r_a)/2;
                        wri2=(r_b-r_a)/2*wri(i2,1);
                        %I_mat=B_mat_func_r1(a,b,ri(j,1),si(i,1)); % 3x8 matrix result
                        I_mat(:,:,j,i)=I_mat(:,:,j,i)+wri2*wsi2*B_mat_func_r1(a,b,0,ri2,si2,0,elem_type)'*D0*B_mat_func_r1(a,b,0,ri2,si2,0,elem_type);
                    end
                end
                I_mat(:,:,j,i) = a*b*I_mat(:,:,j,i); 
                
            end
        end

    case 'QUAD8'
        
        
        [wri,ri]=gausspt(int_scale);
                
        %nu = 0.3;
        %int_scale=1;
        
        D0=1/(1-nu^2)*[1 nu 0
            nu 1 0
            0 0 (1-nu)/2];
        
        % ri=[-1*sqrt(3)/3 +1*sqrt(3)/3]';
        % si=[-1*sqrt(3)/3 +1*sqrt(3)/3]';
        
        % N_r=size(ri,1);
        % N_s=size(si,1);
        
        N_r=DLTOP_ratio;
        N_s=DLTOP_ratio;
        %
        % ri=[-1+2/(N_r)/2:2/(N_r):1]';
        % si=[-1+2/(N_s)/2:2/(N_s):1]';
        
        ri_int=transpose(-1:2/N_r:1);
        si_int=transpose(-1:2/N_s:1);
        %qi_int=transpose(-1:2/N_q:1);
        
        
        % wri=2/N_r/int_scale;
        % wsi=2/N_s/int_scale;
        
        I_mat=zeros(18,18,N_r,N_s);
        %I_mat_temp=zeros(18,18,N_r,N_s);
        for j=1:N_s
            for i=1:N_r
                
                r_a=ri_int(i,1);
                r_b=ri_int(i+1,1);
                s_a=si_int(j,1);
                s_b=si_int(j+1,1);
                
                for j2=1:int_scale
                    si2=(s_b-s_a)/2*ri(j2,1)+(s_b+s_a)/2;
                    wsi2=(s_b-s_a)/2*wri(j2,1);
                    for i2=1:int_scale
                        ri2=(r_b-r_a)/2*ri(i2,1)+(r_b+r_a)/2;
                        wri2=(r_b-r_a)/2*wri(i2,1);
                        %I_mat=B_mat_func_r1(a,b,ri(j,1),si(i,1)); % 3x8 matrix result
                        I_mat(:,:,j,i)=I_mat(:,:,j,i)+wri2*wsi2*B_mat_func_r1(a,b,0,ri2,si2,0,elem_type)'*D0*B_mat_func_r1(a,b,0,ri2,si2,0,elem_type);
                    end
                end
                I_mat(:,:,j,i) = a*b*I_mat(:,:,j,i); 
                
            end
        end
        
    case 'BRICK8'
        
        [wri,ri]=gausspt(int_scale);
                
        %nu = 0.3;
        %int_scale=1;
        
        D0=1/((1+nu)*(1-2*nu))*[1-nu nu nu 0 0 0
            nu 1-nu nu 0 0 0
            nu nu 1-nu 0 0 0
            0 0 0 (1-2*nu)/2 0 0
            0 0 0 0 (1-2*nu)/2 0
            0 0 0 0 0 (1-2*nu)/2];
        
        % ri=[-1*sqrt(3)/3 +1*sqrt(3)/3]';
        % si=[-1*sqrt(3)/3 +1*sqrt(3)/3]';
        
        % N_r=size(ri,1);
        % N_s=size(si,1);
        
        N_r=DLTOP_ratio;
        N_s=DLTOP_ratio;
        N_q=DLTOP_ratio;
        
        % ri=[-1+2/(N_r)/2:2/(N_r):1]';
        % si=[-1+2/(N_s)/2:2/(N_s):1]';
        % qi=[-1+2/(N_q)/2:2/(N_q):1]';
        
        ri_int=transpose(-1:2/N_r:1);
        si_int=transpose(-1:2/N_s:1);
        qi_int=transpose(-1:2/N_q:1);
        
        % wri=2/N_r/int_scale;
        % wsi=2/N_s/int_scale;
        % wqi=2/N_q/int_scale;
        
        
        I_mat=zeros(24,24,N_s,N_r,N_q);
        %I_mat_temp=zeros(24,24,24,N_s,N_r,N_q);
        %I_mat_temp=zeros(24,24,N_s,N_r,N_q);
        for j=1:N_s
            for i=1:N_r
                for k=1:N_q
                    
                    r_a=ri_int(i,1);
                    r_b=ri_int(i+1,1);
                    s_a=si_int(j,1);
                    s_b=si_int(j+1,1);
                    q_a=qi_int(k,1);
                    q_b=qi_int(k+1,1);
                    
                    for j2=1:int_scale
                        si2=(s_b-s_a)/2*ri(j2,1)+(s_b+s_a)/2;
                        wsi2=(s_b-s_a)/2*wri(j2,1);
                        
                        for i2=1:int_scale
                            ri2=(r_b-r_a)/2*ri(i2,1)+(r_b+r_a)/2;
                            wri2=(r_b-r_a)/2*wri(i2,1);
                            
                            for k2=1:int_scale
                                qi2=(q_b-q_a)/2*ri(k2,1)+(q_b+q_a)/2;
                                wqi2=(q_b-q_a)/2*wri(k2,1);
                                
                                %I_mat=B_mat_func_r1(a,b,ri(j,1),si(i,1)); % 3x8 matrix result
                                %I_mat_temp(:,:,j,i,k)=I_mat_temp(:,:,j,i,k)+B_mat_brick8_func_r1(a,b,c,ri_int(i_int,1),si_int(j_int,1),qi_int(q_int,1))'*D0*B_mat_brick8_func_r1(a,b,c,ri_int(i_int,1),si_int(j_int,1),qi_int(q_int,1));
                                I_mat(:,:,j,i,k)=I_mat(:,:,j,i,k)+wri2*wsi2*wqi2*B_mat_func_r1(a,b,c,ri2,si2,qi2,elem_type)'*D0*B_mat_func_r1(a,b,c,ri2,si2,qi2,elem_type);
                            end
                        end
                    end
                    I_mat(:,:,j,i,k) = a*b*c*I_mat(:,:,j,i,k); 
                    
                end
            end
        end
        
    case 'BRICK8_HT' 
        [wri,ri]=gausspt(int_scale);
               
        %nu = 0.3;
        %int_scale=1;
        
        D0=[1 0 0
            0 1 0
            0 0 1];
                
        % ri=[-1*sqrt(3)/3 +1*sqrt(3)/3]';
        % si=[-1*sqrt(3)/3 +1*sqrt(3)/3]';
        
        % N_r=size(ri,1);
        % N_s=size(si,1);
        
        N_r=DLTOP_ratio;
        N_s=DLTOP_ratio;
        N_q=DLTOP_ratio;
        
        % ri=[-1+2/(N_r)/2:2/(N_r):1]';
        % si=[-1+2/(N_s)/2:2/(N_s):1]';
        % qi=[-1+2/(N_q)/2:2/(N_q):1]';
        
        ri_int=transpose(-1:2/N_r:1);
        si_int=transpose(-1:2/N_s:1);
        qi_int=transpose(-1:2/N_q:1);
        
        % wri=2/N_r/int_scale;
        % wsi=2/N_s/int_scale;
        % wqi=2/N_q/int_scale;
                
        I_mat=zeros(8,8,N_s,N_r,N_q);
        %I_mat_temp=zeros(24,24,24,N_s,N_r,N_q);
        %I_mat_temp=zeros(24,24,N_s,N_r,N_q);
        for j=1:N_s
            for i=1:N_r
                for k=1:N_q
                    
                    r_a=ri_int(i,1);
                    r_b=ri_int(i+1,1);
                    s_a=si_int(j,1);
                    s_b=si_int(j+1,1);
                    q_a=qi_int(k,1);
                    q_b=qi_int(k+1,1);
                    
                    for j2=1:int_scale
                        si2=(s_b-s_a)/2*ri(j2,1)+(s_b+s_a)/2;
                        wsi2=(s_b-s_a)/2*wri(j2,1);
                        
                        for i2=1:int_scale
                            ri2=(r_b-r_a)/2*ri(i2,1)+(r_b+r_a)/2;
                            wri2=(r_b-r_a)/2*wri(i2,1);
                            
                            for k2=1:int_scale
                                qi2=(q_b-q_a)/2*ri(k2,1)+(q_b+q_a)/2;
                                wqi2=(q_b-q_a)/2*wri(k2,1);
                                
                                %I_mat=B_mat_func_r1(a,b,ri(j,1),si(i,1)); % 3x8 matrix result
                                %I_mat_temp(:,:,j,i,k)=I_mat_temp(:,:,j,i,k)+B_mat_brick8_func_r1(a,b,c,ri_int(i_int,1),si_int(j_int,1),qi_int(q_int,1))'*D0*B_mat_brick8_func_r1(a,b,c,ri_int(i_int,1),si_int(j_int,1),qi_int(q_int,1));
                                I_mat(:,:,j,i,k)=I_mat(:,:,j,i,k)+wri2*wsi2*wqi2*B_mat_func_r1(a,b,c,ri2,si2,qi2,elem_type)'*D0*B_mat_func_r1(a,b,c,ri2,si2,qi2,elem_type);
                            end
                        end
                    end
                    I_mat(:,:,j,i,k) = a*b*c*I_mat(:,:,j,i,k); 
                    
                end
            end
        end  
        
    case 'BRICK27'
        
        [wri,ri]=gausspt(int_scale);
                
        %nu = 0.3;
        %int_scale=1;
        
        D0=1/((1+nu)*(1-2*nu))*[1-nu nu nu 0 0 0
            nu 1-nu nu 0 0 0
            nu nu 1-nu 0 0 0
            0 0 0 (1-2*nu)/2 0 0
            0 0 0 0 (1-2*nu)/2 0
            0 0 0 0 0 (1-2*nu)/2];
        
        % ri=[-1*sqrt(3)/3 +1*sqrt(3)/3]';
        % si=[-1*sqrt(3)/3 +1*sqrt(3)/3]';
        
        % N_r=size(ri,1);
        % N_s=size(si,1);
        
        N_r=DLTOP_ratio;
        N_s=DLTOP_ratio;
        N_q=DLTOP_ratio;
        
        % ri=[-1+2/(N_r)/2:2/(N_r):1]';
        % si=[-1+2/(N_s)/2:2/(N_s):1]';
        % qi=[-1+2/(N_q)/2:2/(N_q):1]';
        
        ri_int=transpose(-1:2/N_r:1);
        si_int=transpose(-1:2/N_s:1);
        qi_int=transpose(-1:2/N_q:1);
        
        I_mat=zeros(81,81,N_s,N_r,N_q);
        %I_mat_temp=zeros(24,24,24,N_s,N_r,N_q);
        %I_mat_temp=zeros(24,24,N_s,N_r,N_q);
        for j=1:N_s
            for i=1:N_r
                for k=1:N_q
                    
                    r_a=ri_int(i,1);
                    r_b=ri_int(i+1,1);
                    s_a=si_int(j,1);
                    s_b=si_int(j+1,1);
                    q_a=qi_int(k,1);
                    q_b=qi_int(k+1,1);
                    
                    
                    for j2=1:int_scale
                        %j_int=int_scale*(j-1)+j2;
                        si2=(s_b-s_a)/2*ri(j2,1)+(s_b+s_a)/2;
                        wsi2=(s_b-s_a)/2*wri(j2,1);
                        for i2=1:int_scale
                            %i_int=int_scale*(i-1)+i2;
                            ri2=(r_b-r_a)/2*ri(i2,1)+(r_b+r_a)/2;
                            wri2=(r_b-r_a)/2*wri(i2,1);
                            for k2=1:int_scale
                                %q_int=int_scale*(k-1)+k2;
                                qi2=(q_b-q_a)/2*ri(k2,1)+(q_b+q_a)/2;
                                wqi2=(q_b-q_a)/2*wri(k2,1);
                                
                                %I_mat=B_mat_func_r1(a,b,ri(j,1),si(i,1)); % 3x8 matrix result
                                %I_mat_temp(:,:,j,i,k)=I_mat_temp(:,:,j,i,k)+B_mat_brick8_func_r1(a,b,c,ri_int(i_int,1),si_int(j_int,1),qi_int(q_int,1))'*D0*B_mat_brick8_func_r1(a,b,c,ri_int(i_int,1),si_int(j_int,1),qi_int(q_int,1));
                                I_mat(:,:,j,i,k)=I_mat(:,:,j,i,k)+wri2*wsi2*wqi2*B_mat_func_r1(a,b,c,ri2,si2,qi2,elem_type)'*D0*B_mat_func_r1(a,b,c,ri2,si2,qi2,elem_type);
                            end
                        end
                    end
                    I_mat(:,:,j,i,k) = a*b*c*I_mat(:,:,j,i,k); 
                    
                end
            end
        end
        
end

