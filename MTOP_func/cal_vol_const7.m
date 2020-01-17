%%%% const #1 : displacement on the load point %%%%
function [c,ceq,gc,gceq,volfrac_out]=cal_vol_const7(volfrac,nelx,nely,x)

nele=nelx*nely;
sum_eff=0;
%nele_eff=0;
for j=1:nely
    for i=1:nelx
        
         sum_eff=sum_eff+x(j,i);
%         nele_eff=nele_eff+1;
        
    end
end
            
%c=sum_eff-(volfrac)*nele_eff;
c=sum_eff-(volfrac)*nele;
volfrac_out=sum_eff/nele;

%c=sum(x(:).^3)-volfrac*(nelx*nely);
%c=sum(x(:))-(volfrac+x_min/1)*(nelx*nely);
%c=sum(x(:))/(nelx*nely) - volfrac;
%c=sum(x(:))/(volfrac*nelx*nely)-1;

%delc=abs(c-c_old);

ceq=[];

%gc=(nelx*nely)*ones(nely,nelx);
gc=ones(nely,nelx);
%gc=ones(nely,nelx)/(nelx*nely);
%gc2=zeros(nely,nelx);
gceq = [];

end
        
