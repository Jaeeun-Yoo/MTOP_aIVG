function sens_map_mat=sub_sens_map_r4_5(dfdx,data_N)

%intv_data_num=hist(dfdx(:),data_N);
org_N = size(dfdx(:),1);
intv_data_num = zeros(data_N,1);
intv_data_num(:,1) = fix(org_N/data_N);
intv_data_num(end,1) = org_N - (data_N-1)*intv_data_num(1,1);


org_N=size(dfdx,1);
sens_map_mat=zeros(org_N,1);
 
[Y,I] = sort(dfdx,1,'ascend');


%dfdx_resize_weight=dfdx_resize_weight*dfdx_min;
k=0;
for j=1:data_N
    for i=1:intv_data_num(j)
        %sens_map_mat(j,1)=round((dfdx(j,1)-x_min)/x_intv)+1;
        k=k+1;
        map_level=j;
        %map_level=200;
        sens_map_mat(I(k,1),1)=map_level;
        %dfdx_resize(map_level,1)=dfdx_resize(map_level,1)+Y(j,1);
        %dfdx_resize_weight(map_level,1)=dfdx_resize_weight(map_level,1)+1;
    end
    
end

% disp('debug line');


