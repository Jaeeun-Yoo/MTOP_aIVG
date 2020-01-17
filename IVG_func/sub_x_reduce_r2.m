function x_reduce=sub_x_reduce_r2(x,sens_map_mat,data_N)

org_N=size(x,1);
x_reduce=zeros(data_N,1);
x_weight=zeros(data_N,1);


for j=1:org_N
    x_reduce(sens_map_mat(j,1),1)=x_reduce(sens_map_mat(j,1),1)+x(j,1);
    x_weight(sens_map_mat(j,1),1)=x_weight(sens_map_mat(j,1),1)+1;
    %save('debug_sub.mat','x_reduce','x_weight','j','sens_map_mat');
end

for j=1:data_N
    if x_weight(j,1)==0
        x_reduce(j,1)=x_reduce(j,1);
    else
        x_reduce(j,1)=x_reduce(j,1)/x_weight(j,1);
    end
end

%disp('test');


