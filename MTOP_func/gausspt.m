% Gauss-Legendre quadrature point / weight
% Number of point : 2,3,4
% '16.11.12
% Jaeeun Yoo


function [wri,ri]=gausspt(int_scale)

if int_scale==2
    wri=[1;1];
    ri=[-sqrt(1/3);sqrt(1/3)];
   
elseif int_scale==3
    wri=[5/9;8/9;5/9];
    ri=[-sqrt(3/5);0;sqrt(3/5)];
    
elseif int_scale==4
    wri=[(18-sqrt(30))/36;(18+sqrt(30))/36;(18+sqrt(30))/36;(18-sqrt(30))/36];
    ri=[-sqrt(3/7+2/7*sqrt(6/5));-sqrt(3/7-2/7*sqrt(6/5));sqrt(3/7-2/7*sqrt(6/5));sqrt(3/7+2/7*sqrt(6/5))];
    
elseif int_scale==5
    
    wri=[(322-13*sqrt(70))/900;(322+13*sqrt(70))/900;128/225;(322+13*sqrt(70))/900;(322-13*sqrt(70))/900];
    ri=[-1/3*sqrt(5+2*sqrt(10/7));-1/3*sqrt(5-2*sqrt(10/7));0;1/3*sqrt(5-2*sqrt(10/7));1/3*sqrt(5+2*sqrt(10/7))];
else
    disp('Gauss point ERROR / no valid number of points')
end
