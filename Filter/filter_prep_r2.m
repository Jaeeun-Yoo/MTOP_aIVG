function [H,Hs]=filter_prep_r2(id_2d,inp_filter_prep,x_d)

N_xd=size(x_d,1);

if id_2d==1

%     nele = inp_filter_prep.nele;
    rmin = inp_filter_prep.rmin;
    nelx_d = inp_filter_prep.nelx_d;
    nely_d = inp_filter_prep.nely_d;
    
    step = ceil(rmin)-1;
    iH = zeros(N_xd*(2*step+1)^2,1);
    jH = zeros(size(iH)); vH = zeros(size(iH));
    n = 0;
    for kk=1:N_xd
        el=x_d(kk,1);
        [j,i] = ind2sub([nely_d,nelx_d],el);
        [ispan,jspan] = meshgrid(max(1,i-step):min(nelx_d,i+step),max(1,j-step):min(nely_d,j+step));
        dist = max(0,rmin-sqrt((ispan-i).^2 + (jspan-j).^2 ));
        vH(n+(1:numel(dist))) = dist(:)./rmin;
        iH(n+(1:numel(dist))) = el;
        %jH(n+(1:numel(dist))) = sub2ind([nely_d nelx_d],ispan,jspan);
        jH(n+(1:numel(dist))) = sub2ind([nely_d nelx_d],jspan,ispan);
        n = n + numel(dist);
    end
    eff_idx = vH > 0;
    iH = iH(eff_idx);
    jH = jH(eff_idx);
    vH = vH(eff_idx);
    
%     iH(n+1:end)=[]; jH(n+1:end)=[]; vH(n+1:end)=[];
%     H = sparse(iH,jH,vH);
    H_temp = sparse(iH,jH,vH);
    H = sparse(H_temp(x_d,x_d));
    Hs = sum(H,2);
    clear iH
    clear jH
    clear vH
    clear H_temp
    clear eff_idx
        
else
    
%     nele = inp_filter_prep.nele;
    rmin = inp_filter_prep.rmin;
    nelx_d = inp_filter_prep.nelx_d;
    nely_d = inp_filter_prep.nely_d;
    nelz_d = inp_filter_prep.nelz_d;
    
    step = ceil(rmin)-1;
    iH = zeros(N_xd*(2*step+1)^3,1);
    jH = zeros(size(iH)); vH = zeros(size(iH));
    n = 0;
    for kk=1:N_xd
        el=x_d(kk,1);
        [j,i,k] = ind2sub([nely_d,nelx_d,nelz_d],el);
        [ispan,jspan,kspan] = meshgrid(max(1,i-step):min(nelx_d,i+step),max(1,j-step):min(nely_d,j+step),max(1,k-step):min(nelz_d,k+step));
        dist = max(0,rmin-sqrt((ispan-i).^2 + (jspan-j).^2 + (kspan-k).^2));
        vH(n+(1:numel(dist))) = dist(:)./rmin;
        iH(n+(1:numel(dist))) = el;
        jH(n+(1:numel(dist))) = sub2ind([nely_d nelx_d nelz_d],jspan,ispan,kspan);
        n = n + numel(dist);
    end
    eff_idx = vH > 0;
    iH = iH(eff_idx);
    jH = jH(eff_idx);
    vH = vH(eff_idx);    
    
    H_temp = sparse(iH,jH,vH);
    H = sparse(H_temp(x_d,x_d));
%     H = sparse(iH,jH,vH);
    Hs = sum(H,2);
    clear iH
    clear jH
    clear vH
    clear H_temp
    clear eff_idx
    
end



disp('filter preparation done2')


