function x_recover=sub_x_recover_r1(x,sens_map_mat,x_reduce)

org_N=size(x,1);
x_recover=zeros(org_N,1);

for j=1:org_N
    x_recover(j,1)=x_reduce(sens_map_mat(j,1),1);
   
end



