function [x,x_d,x_nd_s,x_nd_v,new_volfrac_ini]=Define_design_domain_r2(nelx_d,nely_d,nelz_d,x,volfrac_ini,elem_type,ex_type)

%% Design / non-design domain
%nele
%x(1:nely_d,1:nelx_d,1:nelz_d) = volfrac_ini; 

switch ex_type
    case 'bridge'
        switch elem_type
            case {'QUAD4','QUAD8','QUAD4_HT'}
                x_nd_s_idx=zeros(nely_d,nelx_d); % non-design solid, density=1
                x_nd_v_idx=zeros(nely_d,nelx_d); % non-design void, density=0
                
                x_nd_s_cnt=0; % number of non-design solid x element
                x_nd_v_cnt=0; % number of non-design void x element
                
                x_nd_s=find(x_nd_s_idx);
                x_nd_v=find(x_nd_v_idx);
                
                
                all_x = [1:nely_d*nelx_d]';
                x_d_temp = setdiff(all_x,x_nd_s);
                x_d = setdiff(x_d_temp,x_nd_v);
                
                % initial design varialbe compensation by non-design variables
                total_variable_num = nely_d*nelx_d;
                volfrac_w_nd = (volfrac_ini*(total_variable_num-size(x_nd_s,1)-size(x_nd_v,1))+size(x_nd_s,1))/total_variable_num;
                x_temp=x(:)*volfrac_ini/volfrac_w_nd;
                x_temp(x_nd_s)=1;
                x_temp(x_nd_v)=0;
                x=reshape(x_temp,nely_d,nelx_d);
                
                disp('NDx Dx allocation... (2D)')
                                
            case {'BRICK8','BRICK27','BRICK8_HT'}
                x_nd_s_idx=zeros(nely_d,nelx_d,nelz_d); % non-design solid, density=1
                x_nd_v_idx=zeros(nely_d,nelx_d,nelz_d); % non-design void, density=0
                
                x_nd_s_cnt=0; % number of non-design solid x element
                x_nd_v_cnt=0; % number of non-design void x element
                for k=1:nelz_d  % solid design area
                    for j=1:nely_d
                        for i=1:nelx_d
                            if j>=15/30*nely_d+1 && j<=16/30*nely_d
                                x_nd_s_cnt=x_nd_s_cnt+1;
                                x_nd_s_idx(j,i,k)=1;
                            end
                        end
                    end
                end
                
                for k=1:nelz_d % void design area
                    for j=1:nely_d
                        for i=1:nelx_d
                            %if j>=6/30*nely_d-(DLTOP_ratio-1) && j<=15/30*nely_d
                            if j>=1 && j<=15/30*nely_d
                                %if k>=3/10*nelz_d-(DLTOP_ratio-1) && k<=8/10*nelz_d
                                if k>=1/5*nelz_d+1 %&& k<=9/10*nelz_d
                                    x_nd_v_cnt=x_nd_v_cnt+1;
                                    x_nd_v_idx(j,i,k)=1;
                                end
                            end
                        end
                    end
                end
                
                %x_nd_s=zeros(x_nd_s_cnt,3);
                %[x_nd_s(:,1),x_nd_s(:,2),x_nd_s(:,3)]=find(x_nd_s_idx);
                %[x_nd_v(:,1),x_nd_v(:,2),x_nd_v(:,3)]=find(x_nd_v_idx);
                x_nd_s=find(x_nd_s_idx);
                x_nd_v=find(x_nd_v_idx);
                
                
                all_x = [1:nely_d*nelx_d*nelz_d]';
                x_d_temp = setdiff(all_x,x_nd_s);
                x_d = setdiff(x_d_temp,x_nd_v);
                
                % initial design varialbe compensation by non-design variables
                total_variable_num = nely_d*nelx_d*nelz_d;
                volfrac_w_nd = (volfrac_ini*(total_variable_num-size(x_nd_s,1)-size(x_nd_v,1))+size(x_nd_s,1))/total_variable_num;
                x_temp=x(:)*volfrac_ini/volfrac_w_nd;
                x_temp(x_nd_s)=1;
                x_temp(x_nd_v)=0;
                x=reshape(x_temp,nely_d,nelx_d,nelz_d);
                
                disp('NDx Dx allocation... (3D)')
                
        end
    
        
    case {'bridge_case2','bridge_case3'}
        switch elem_type
            case {'QUAD4','QUAD8','QUAD4_HT'}
                x_nd_s_idx=zeros(nely_d,nelx_d); % non-design solid, density=1
                x_nd_v_idx=zeros(nely_d,nelx_d); % non-design void, density=0
                
                x_nd_s_cnt=0; % number of non-design solid x element
                x_nd_v_cnt=0; % number of non-design void x element
                
                x_nd_s=find(x_nd_s_idx);
                x_nd_v=find(x_nd_v_idx);
                
                
                all_x = [1:nely_d*nelx_d]';
                x_d_temp = setdiff(all_x,x_nd_s);
                x_d = setdiff(x_d_temp,x_nd_v);
                
                % initial design varialbe compensation by non-design variables
                total_variable_num = nely_d*nelx_d;
                volfrac_w_nd = (volfrac_ini*(total_variable_num-size(x_nd_s,1)-size(x_nd_v,1))+size(x_nd_s,1))/total_variable_num;
                x_temp=x(:)*volfrac_ini/volfrac_w_nd;
                x_temp(x_nd_s)=1;
                x_temp(x_nd_v)=0;
                x=reshape(x_temp,nely_d,nelx_d);
                
                disp('NDx Dx allocation... (2D)')
                
            case {'BRICK8','BRICK27','BRICK8_HT'}
                x_nd_s_idx=zeros(nely_d,nelx_d,nelz_d); % non-design solid, density=1
                x_nd_v_idx=zeros(nely_d,nelx_d,nelz_d); % non-design void, density=0
                
                x_nd_s_cnt=0; % number of non-design solid x element
                x_nd_v_cnt=0; % number of non-design void x element
                for k=1:nelz_d  % solid design area
                    for j=1:nely_d
                        for i=1:nelx_d
%                             if j>=1 && j<=1/15*nely_d % Nguyen ex
%                             if j>=1 && j<=1/44*nely_d  % Zegard ex
                            if j>=1 && j<=1/22*nely_d  % Zegard ex, '20.08.09 modified
                                x_nd_s_cnt=x_nd_s_cnt+1;
                                x_nd_s_idx(j,i,k)=1;
                            end
                        end
                    end
                end
                
%                 for k=1:nelz_d % void design area
%                     for j=1:nely_d
%                         for i=1:nelx_d
%                             %if j>=6/30*nely_d-(DLTOP_ratio-1) && j<=15/30*nely_d
%                             if j>=1 && j<=15/30*nely_d
%                                 %if k>=3/10*nelz_d-(DLTOP_ratio-1) && k<=8/10*nelz_d
%                                 if k>=1/10*nelz_d+1 && k<=9/10*nelz_d
%                                     x_nd_v_cnt=x_nd_v_cnt+1;
%                                     x_nd_v_idx(j,i,k)=1;
%                                 end
%                             end
%                         end
%                     end
%                 end
%                 
                %x_nd_s=zeros(x_nd_s_cnt,3);
                %[x_nd_s(:,1),x_nd_s(:,2),x_nd_s(:,3)]=find(x_nd_s_idx);
                %[x_nd_v(:,1),x_nd_v(:,2),x_nd_v(:,3)]=find(x_nd_v_idx);
                x_nd_s=find(x_nd_s_idx);
                x_nd_v=find(x_nd_v_idx);
                
                
                all_x = [1:nely_d*nelx_d*nelz_d]';
                x_d_temp = setdiff(all_x,x_nd_s);
                x_d = setdiff(x_d_temp,x_nd_v);
                
                % initial design varialbe compensation by non-design variables
                total_variable_num = nely_d*nelx_d*nelz_d;
                volfrac_w_nd = (volfrac_ini*(total_variable_num-size(x_nd_s,1)-size(x_nd_v,1))+size(x_nd_s,1))/total_variable_num;
                x_temp=x(:)*volfrac_ini/volfrac_w_nd;
                x_temp(x_nd_s)=1;
                x_temp(x_nd_v)=0;
                x=reshape(x_temp,nely_d,nelx_d,nelz_d);
                
                disp('NDx Dx allocation... (3D)')
                
        end
        
    case 'michell'
        switch elem_type
            case {'QUAD4','QUAD8','QUAD4_HT'}
                x_nd_s_idx=zeros(nely_d,nelx_d); % non-design solid, density=1
                x_nd_v_idx=zeros(nely_d,nelx_d); % non-design void, density=0
                
                x_nd_s_cnt=0; % number of non-design solid x element
                x_nd_v_cnt=0; % number of non-design void x element
                
                cr = (nelx_d*30/300);
                cx = (nelx_d*30/300)+0.5;
                cy = (nely_d/2)+0.5;
                
                % solid area
                for j=1:nely_d
                    for i=1:nelx_d
                        if (i-cx)^2+(j-cy)^2<=cr^2
                            x_nd_s_cnt=x_nd_s_cnt+1;
                            x_nd_s_idx(j,i)=1;
                        end
                    end
                end
                
                x_nd_s=find(x_nd_s_idx);
                x_nd_v=find(x_nd_v_idx);
                
                all_x = [1:nely_d*nelx_d]';
                x_d_temp = setdiff(all_x,x_nd_s);
                x_d = setdiff(x_d_temp,x_nd_v);
                
                % initial design varialbe compensation by non-design variables
                total_variable_num = nely_d*nelx_d;
                volfrac_w_nd = (volfrac_ini*(total_variable_num-size(x_nd_s,1)-size(x_nd_v,1))+size(x_nd_s,1))/total_variable_num;
                x_temp=x(:)*volfrac_ini/volfrac_w_nd;
                x_temp(x_nd_s)=1;
                x_temp(x_nd_v)=0;
                x=reshape(x_temp,nely_d,nelx_d);
                
                disp('NDx Dx allocation... (2D)')
                                
            case {'BRICK8','BRICK27','BRICK8_HT'}
                x_nd_s_idx=zeros(nely_d,nelx_d,nelz_d); % non-design solid, density=1
                x_nd_v_idx=zeros(nely_d,nelx_d,nelz_d); % non-design void, density=0
                
                x_nd_s_cnt=0; % number of non-design solid x element
                x_nd_v_cnt=0; % number of non-design void x element
%                 for k=1:nelz_d  % solid design area
%                     for j=1:nely_d
%                         for i=1:nelx_d
%                             if j>=15/30*nely_d+1 && j<=16/30*nely_d
%                                 x_nd_s_cnt=x_nd_s_cnt+1;
%                                 x_nd_s_idx(j,i,k)=1;
%                             end
%                         end
%                     end
%                 end
%                 
%                 for k=1:nelz_d % void design area
%                     for j=1:nely_d
%                         for i=1:nelx_d
%                             %if j>=6/30*nely_d-(DLTOP_ratio-1) && j<=15/30*nely_d
%                             if j>=1 && j<=15/30*nely_d
%                                 %if k>=3/10*nelz_d-(DLTOP_ratio-1) && k<=8/10*nelz_d
%                                 if k>=1/5*nelz_d+1 %&& k<=9/10*nelz_d
%                                     x_nd_v_cnt=x_nd_v_cnt+1;
%                                     x_nd_v_idx(j,i,k)=1;
%                                 end
%                             end
%                         end
%                     end
%                 end
                
                %x_nd_s=zeros(x_nd_s_cnt,3);
                %[x_nd_s(:,1),x_nd_s(:,2),x_nd_s(:,3)]=find(x_nd_s_idx);
                %[x_nd_v(:,1),x_nd_v(:,2),x_nd_v(:,3)]=find(x_nd_v_idx);
                x_nd_s=find(x_nd_s_idx);
                x_nd_v=find(x_nd_v_idx);
                
                
                all_x = [1:nely_d*nelx_d*nelz_d]';
                x_d_temp = setdiff(all_x,x_nd_s);
                x_d = setdiff(x_d_temp,x_nd_v);
                
                % initial design varialbe compensation by non-design variables
                total_variable_num = nely_d*nelx_d*nelz_d;
                volfrac_w_nd = (volfrac_ini*(total_variable_num-size(x_nd_s,1)-size(x_nd_v,1))+size(x_nd_s,1))/total_variable_num;
                x_temp=x(:)*volfrac_ini/volfrac_w_nd;
                x_temp(x_nd_s)=1;
                x_temp(x_nd_v)=0;
                x=reshape(x_temp,nely_d,nelx_d,nelz_d);
                
                disp('NDx Dx allocation... (3D)')
                
        end
    
    case 'optics_mount2'
        switch elem_type
            case {'QUAD4','QUAD8','QUAD4_HT'}
                x_nd_s_idx=zeros(nely_d,nelx_d); % non-design solid, density=1
                x_nd_v_idx=zeros(nely_d,nelx_d); % non-design void, density=0
                
                x_nd_s_cnt=0; % number of non-design solid x element
                x_nd_v_cnt=0; % number of non-design void x element
                
                x_nd_s=find(x_nd_s_idx);
                x_nd_v=find(x_nd_v_idx);
                
                
                all_x = [1:nely_d*nelx_d]';
                x_d_temp = setdiff(all_x,x_nd_s);
                x_d = setdiff(x_d_temp,x_nd_v);
                
                % initial design varialbe compensation by non-design variables
                total_variable_num = nely_d*nelx_d;
                volfrac_w_nd = (volfrac_ini*(total_variable_num-size(x_nd_s,1)-size(x_nd_v,1))+size(x_nd_s,1))/total_variable_num;
                x_temp=x(:)*volfrac_ini/volfrac_w_nd;
                x_temp(x_nd_s)=1;
                x_temp(x_nd_v)=0;
                x=reshape(x_temp,nely_d,nelx_d);
                
                disp('NDx Dx allocation... (2D)')
                
            case {'BRICK8','BRICK27','BRICK8_HT'}
                x_nd_s_idx=zeros(nely_d,nelx_d,nelz_d); % non-design solid, density=1
                x_nd_v_idx=zeros(nely_d,nelx_d,nelz_d); % non-design void, density=0
                
                x_nd_s_cnt=0; % number of non-design solid x element
                x_nd_v_cnt=0; % number of non-design void x element
                for k=1:nelz_d  % solid non-design area
                    for j=1:nely_d
                        for i=1:nelx_d
                            % #1, Bonding area
                            if i>=21/66*nelx_d+1 && i<=45/66*nelx_d
                                if j<=16/36*nely_d
                                    if k<=2/12*nelz_d
                                        x_nd_s_cnt=x_nd_s_cnt+1;
                                        x_nd_s_idx(j,i,k)=1;
                                    end
                                end
                            end
                            
                            % #2, Bolting area, Left bottom
                            if i<=16/66*nelx_d
                                if j>=34/36*nely_d+1
                                    if k<=nelz_d
                                        x_nd_s_cnt=x_nd_s_cnt+1;
                                        x_nd_s_idx(j,i,k)=1;
                                    end
                                end
                            end
                            
                            % #3, Bolting area, Right bottom
                            if i>=50/66*nelx_d+1
                                if j>=34/36*nely_d+1
                                    if k<=nelz_d
                                        x_nd_s_cnt=x_nd_s_cnt+1;
                                        x_nd_s_idx(j,i,k)=1;
                                    end
                                end
                            end
                            
                        end
                    end
                end
                
                for i=1:nelx_d/2 % void non-design area, left
                    for j=1:nely_d
                        for k=1:nelz_d
                            if k>((i-20/66*nelx_d)*tan(60*pi()/180))
                                if j<=23/36*nely_d
                                    x_nd_v_cnt=x_nd_v_cnt+2;
                                    x_nd_v_idx(j,i,k)=1;
                                    x_nd_v_idx(j,nelx_d-i+1,k)=1;
                                end
                            end
                            %
                            %                             if k>(-(i-48/66*nelx_d)*tan(60*pi()/180))
                            %                                 if j<=23/36*nely_d
                            %                                     x_nd_v_cnt=x_nd_v_cnt+1;
                            %                                     x_nd_v_idx(j,i,k)=1;
                            %                                 end
                            %                             end
                        end
                    end
                end
                
                %x_nd_s=zeros(x_nd_s_cnt,3);
                %[x_nd_s(:,1),x_nd_s(:,2),x_nd_s(:,3)]=find(x_nd_s_idx);
                %[x_nd_v(:,1),x_nd_v(:,2),x_nd_v(:,3)]=find(x_nd_v_idx);
                x_nd_s=find(x_nd_s_idx);
                x_nd_v=find(x_nd_v_idx);
                
                
                all_x = [1:nely_d*nelx_d*nelz_d]';
                x_d_temp = setdiff(all_x,x_nd_s);
                x_d = setdiff(x_d_temp,x_nd_v);
                
                % initial design varialbe compensation by non-design variables
                total_variable_num = nely_d*nelx_d*nelz_d;
                volfrac_w_nd = (volfrac_ini*(total_variable_num-size(x_nd_s,1)-size(x_nd_v,1))+size(x_nd_s,1))/total_variable_num;
                x_temp=x(:)*volfrac_ini/volfrac_w_nd;
                x_temp(x_nd_s)=1;
                x_temp(x_nd_v)=0;
                x=reshape(x_temp,nely_d,nelx_d,nelz_d);
                
                disp('NDx Dx allocation... (3D)')
                
        end
    
        
        
    otherwise
        switch elem_type
            case {'QUAD4','QUAD8','QUAD4_HT'}
                x_nd_s_idx=zeros(nely_d,nelx_d); % non-design solid, density=1
                x_nd_v_idx=zeros(nely_d,nelx_d); % non-design void, density=0
                
                x_nd_s_cnt=0; % number of non-design solid x element
                x_nd_v_cnt=0; % number of non-design void x element
                x_nd_s=find(x_nd_s_idx);
                x_nd_v=find(x_nd_v_idx);
                
                
                all_x = [1:nely_d*nelx_d]';
                x_d_temp = setdiff(all_x,x_nd_s);
                x_d = setdiff(x_d_temp,x_nd_v);
                
                % initial design varialbe compensation by non-design variables
                total_variable_num = nely_d*nelx_d;
                volfrac_w_nd = (volfrac_ini*(total_variable_num-size(x_nd_s,1)-size(x_nd_v,1))+size(x_nd_s,1))/total_variable_num;
                x_temp=x(:)*volfrac_ini/volfrac_w_nd;
                x_temp(x_nd_s)=1;
                x_temp(x_nd_v)=0;
                x=reshape(x_temp,nely_d,nelx_d);
                
                disp('NDx Dx allocation... (2D)')
                
            case {'BRICK8','BRICK27','BRICK8_HT'}
                
                x_nd_s_idx=zeros(nely_d,nelx_d,nelz_d); % non-design solid, density=1
                x_nd_v_idx=zeros(nely_d,nelx_d,nelz_d); % non-design void, density=0
                
                x_nd_s_cnt=0; % number of non-design solid x element
                x_nd_v_cnt=0; % number of non-design void x element
                x_nd_s=find(x_nd_s_idx);
                x_nd_v=find(x_nd_v_idx);
                
                
                all_x = [1:nely_d*nelx_d*nelz_d]';
                x_d_temp = setdiff(all_x,x_nd_s);
                x_d = setdiff(x_d_temp,x_nd_v);
                
                % initial design varialbe compensation by non-design variables
                total_variable_num = nely_d*nelx_d*nelz_d;
                volfrac_w_nd = (volfrac_ini*(total_variable_num-size(x_nd_s,1)-size(x_nd_v,1))+size(x_nd_s,1))/total_variable_num;
                x_temp=x(:)*volfrac_ini/volfrac_w_nd;
                x_temp(x_nd_s)=1;
                x_temp(x_nd_v)=0;
                x=reshape(x_temp,nely_d,nelx_d,nelz_d);
                
                disp('NDx Dx allocation...(3D)')
                                
        end
                
end

new_volfrac_ini = volfrac_ini*volfrac_ini/volfrac_w_nd;

