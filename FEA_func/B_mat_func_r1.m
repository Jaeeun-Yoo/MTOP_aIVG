
%%%%%%%%%% ELEMENT STIFFNESS MATRIX %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function [B0]=B_mat_func_r1(a,b,c,r,s,q,elem_type)

  
switch (elem_type) 
    case 'QUAD4' 
        b=-b; % '18.12.16 update, y axis direction compensation
                
        % '18.12.16 original version
        B0 = [(s-1)/(4*a) 0 (1-s)/(4*a) 0 (1+s)/(4*a) 0 -(1+s)/(4*a) 0
            0 (r-1)/(4*b) 0 -(1+r)/(4*b) 0 (1+r)/(4*b) 0 (1-r)/(4*b)
            (r-1)/(4*b) (s-1)/(4*a) -(1+r)/(4*b) (1-s)/(4*a) (1+r)/(4*b) (1+s)/(4*a) (1-r)/(4*b) -(1+s)/(4*a)];
        % '19.08.25 original version
%         B0 = [(s-1)/(4) 0 (1-s)/(4) 0 (1+s)/(4) 0 -(1+s)/(4) 0
%             0 (r-1)/(4) 0 -(1+r)/(4) 0 (1+r)/(4) 0 (1-r)/(4)
%             (r-1)/(4) (s-1)/(4) -(1+r)/(4) (1-s)/(4) (1+r)/(4) (1+s)/(4) (1-r)/(4) -(1+s)/(4)];

%           B0 = [(s-b)/(4*a*b) 0 (b-s)/(4*a*b) 0 (b+s)/(4*a*b) 0 -(b+s)/(4*a*b) 0
%             0 (r-a)/(4*a*b) 0 -(a+r)/(4*a*b) 0 (a+r)/(4*a*b) 0 (a-r)/(4*a*b)
%             (r-a)/(4*a*b) (s-b)/(4*a*b) -(a+r)/(4*a*b) (b-s)/(4*a*b) (a+r)/(4*a*b) (b+s)/(4*a*b) (a-r)/(4*a*b) -(b+s)/(4*a*b)];
     
    case 'QUAD4_HT'
        b=-b; % '18.12.16 update, y axis direction compensation
        
        % '19.01.14 heat transfer 2D
        B0 = [(s-1)/(4*a) (1-s)/(4*a) (1+s)/(4*a) -(1+s)/(4*a)
            (r-1)/(4*b) -(1+r)/(4*b) (1+r)/(4*b) (1-r)/(4*b)];

    case 'QUAD8'
        b=-b; % '18.12.16 update, y axis direction compensation
        
        dN1dx=1/(4*a)*(2*r-1)*(s^2-s);
        dN2dx=1/(4*a)*(2*r+1)*(s^2-s);
        dN3dx=1/(4*a)*(2*r+1)*(s^2+s);
        dN4dx=1/(4*a)*(2*r-1)*(s^2+s);
        dN5dx=-1/(a)*r*(s^2-s);
        dN6dx=1/(2*a)*(2*r+1)*(1-s^2);
        dN7dx=-1/(a)*r*(s^2+s);
        dN8dx=1/(2*a)*(2*r-1)*(1-s^2);
        dN9dx=-2/(a)*r*(1-s^2);
        
        dN1dy=1/(4*b)*(r^2-r)*(2*s-1);
        dN2dy=1/(4*b)*(r^2+r)*(2*s-1);
        dN3dy=1/(4*b)*(r^2+r)*(2*s+1);
        dN4dy=1/(4*b)*(r^2-r)*(2*s+1);
        dN5dy=1/(2*b)*(1-r^2)*(2*s-1);
        dN6dy=-1/(b)*(r^2+r)*s;
        dN7dy=1/(2*b)*(1-r^2)*(2*s+1);
        dN8dy=-1/(b)*(r^2-r)*s;
        dN9dy=-2/(b)*(1-r^2)*s;
                
        
        B0 = [dN1dx 0 dN2dx 0 dN3dx 0 dN4dx 0 dN5dx 0 dN6dx 0 dN7dx 0 dN8dx 0 dN9dx 0
            0 dN1dy 0 dN2dy 0 dN3dy 0 dN4dy 0 dN5dy 0 dN6dy 0 dN7dy 0 dN8dy 0 dN9dy
            dN1dy dN1dx dN2dy dN2dx dN3dy dN3dx dN4dy dN4dx dN5dy dN5dx dN6dy dN6dx dN7dy dN7dx dN8dy dN8dx dN9dy dN9dx];
        
    case 'BRICK8'
        b=-b; % '18.12.16 update, y axis direction compensation
                
        B0 = [(s-1)*(1-q)/(8*a) 0 0 (1-s)*(1-q)/(8*a) 0 0 (1+s)*(1-q)/(8*a) 0 0 -(1+s)*(1-q)/(8*a) 0 0 (s-1)*(1+q)/(8*a) 0 0 (1-s)*(1+q)/(8*a) 0 0 (1+s)*(1+q)/(8*a) 0 0 -(1+s)*(1+q)/(8*a) 0 0
            0 (r-1)*(1-q)/(8*b) 0 0 -(1+r)*(1-q)/(8*b) 0 0 (1+r)*(1-q)/(8*b) 0 0 (1-r)*(1-q)/(8*b) 0 0 (r-1)*(1+q)/(8*b) 0 0 -(1+r)*(1+q)/(8*b) 0 0 (1+r)*(1+q)/(8*b) 0 0 (1-r)*(1+q)/(8*b) 0
            0 0 -(1-r)*(1-s)/(8*c) 0 0 -(1+r)*(1-s)/(8*c) 0 0 -(1+r)*(1+s)/(8*c) 0 0 -(1-r)*(1+s)/(8*c) 0 0 (1-r)*(1-s)/(8*c) 0 0 (1+r)*(1-s)/(8*c) 0 0 (1+r)*(1+s)/(8*c) 0 0 (1-r)*(1+s)/(8*c)
            (r-1)*(1-q)/(8*b) (s-1)*(1-q)/(8*a) 0 -(1+r)*(1-q)/(8*b) (1-s)*(1-q)/(8*a) 0 (1+r)*(1-q)/(8*b) (1+s)*(1-q)/(8*a) 0 (1-r)*(1-q)/(8*b) -(1+s)*(1-q)/(8*a) 0 (r-1)*(1+q)/(8*b) (s-1)*(1+q)/(8*a) 0 -(1+r)*(1+q)/(8*b) (1-s)*(1+q)/(8*a) 0 (1+r)*(1+q)/(8*b) (1+s)*(1+q)/(8*a) 0 (1-r)*(1+q)/(8*b) -(1+s)*(1+q)/(8*a) 0
            -(1-r)*(1-s)/(8*c) 0 -(1-s)*(1-q)/(8*a) -(1+r)*(1-s)/(8*c) 0 (1-s)*(1-q)/(8*a) -(1+r)*(1+s)/(8*c) 0 (1+s)*(1-q)/(8*a) -(1-r)*(1+s)/(8*c) 0 -(1+s)*(1-q)/(8*a) (1-r)*(1-s)/(8*c) 0 -(1-s)*(1+q)/(8*a) (1+r)*(1-s)/(8*c) 0 (1-s)*(1+q)/(8*a) (1+r)*(1+s)/(8*c) 0 (1+s)*(1+q)/(8*a) (1-r)*(1+s)/(8*c) 0 -(1+s)*(1+q)/(8*a)
            0 -(1-r)*(1-s)/(8*c) -(1-r)*(1-q)/(8*b) 0 -(1+r)*(1-s)/(8*c) -(1+r)*(1-q)/(8*b) 0 -(1+r)*(1+s)/(8*c) (1+r)*(1-q)/(8*b) 0 -(1-r)*(1+s)/(8*c) (1-r)*(1-q)/(8*b) 0 (1-r)*(1-s)/(8*c) -(1-r)*(1+q)/(8*b) 0 (1+r)*(1-s)/(8*c) -(1+r)*(1+q)/(8*b) 0 (1+r)*(1+s)/(8*c) (1+r)*(1+q)/(8*b) 0 (1-r)*(1+s)/(8*c) (1-r)*(1+q)/(8*b)];
        
    case 'BRICK8_HT'
        b=-b; % '18.12.16 update, y axis direction compensation
        
        % '19.01.14 heat transfer 2D
        B0 = [(s-1)*(1-q)/(8*a) (1-s)*(1-q)/(8*a) (1+s)*(1-q)/(8*a) -(1+s)*(1-q)/(8*a) (s-1)*(1+q)/(8*a) (1-s)*(1+q)/(8*a) (1+s)*(1+q)/(8*a) -(1+s)*(1+q)/(8*a)
            (r-1)*(1-q)/(8*b) -(1+r)*(1-q)/(8*b) (1+r)*(1-q)/(8*b) (1-r)*(1-q)/(8*b) (r-1)*(1+q)/(8*b) -(1+r)*(1+q)/(8*b) (1+r)*(1+q)/(8*b) (1-r)*(1+q)/(8*b)
            -(1-r)*(1-s)/(8*c) -(1+r)*(1-s)/(8*c) -(1+r)*(1+s)/(8*c) -(1-r)*(1+s)/(8*c) (1-r)*(1-s)/(8*c) (1+r)*(1-s)/(8*c) (1+r)*(1+s)/(8*c) (1-r)*(1+s)/(8*c)];
        
        
    case 'BRICK27'
        b=-b; % '18.12.16 update, y axis direction compensation
        
        % basis function
        %
        % N1=1/(8)*(r^2-r)*(s^2-s)*(q^2-q);
        % N2=1/(8)*(r^2+r)*(s^2-s)*(q^2-q);
        % N3=1/(8)*(r^2+r)*(s^2+s)*(q^2-q);
        % N4=1/(8)*(r^2-r)*(s^2+s)*(q^2-q);
        %
        % N5=1/(8)*(r^2-r)*(s^2-s)*(q^2+q);
        % N6=1/(8)*(r^2+r)*(s^2-s)*(q^2+q);
        % N7=1/(8)*(r^2+r)*(s^2+s)*(q^2+q);
        % N8=1/(8)*(r^2-r)*(s^2+s)*(q^2+q);
        %
        % N9=-1/(4)*(r^2-1)*(s^2-s)*(q^2-q);
        % N10=-1/(4)*(r^2+r)*(s^2-1)*(q^2-q);
        % N11=-1/(4)*(r^2-1)*(s^2+s)*(q^2-q);
        % N12=-1/(4)*(r^2-r)*(s^2-1)*(q^2-q);
        %
        % N13=-1/(4)*(r^2-r)*(s^2-s)*(q^2-1);
        % N14=-1/(4)*(r^2+r)*(s^2-s)*(q^2-1);
        % N15=-1/(4)*(r^2+r)*(s^2+s)*(q^2-1);
        % N16=-1/(4)*(r^2-r)*(s^2+s)*(q^2-1);
        %
        % N17=-1/(4)*(r^2-1)*(s^2-s)*(q^2+q);
        % N18=-1/(4)*(r^2+r)*(s^2-1)*(q^2+q);
        % N19=-1/(4)*(r^2-1)*(s^2+s)*(q^2+q);
        % N20=-1/(4)*(r^2-r)*(s^2-1)*(q^2+q);
        %
        % N21=1/(2)*(r^2-1)*(s^2-1)*(q^2-q);
        % N22=1/(2)*(r^2-1)*(s^2-s)*(q^2-1);
        % N23=1/(2)*(r^2+r)*(s^2-1)*(q^2-1);
        % N24=1/(2)*(r^2-1)*(s^2+s)*(q^2-1);
        % N25=1/(2)*(r^2-r)*(s^2-1)*(q^2-1);
        % N26=1/(2)*(r^2-1)*(s^2-1)*(q^2+q);
        %
        % N27=-1/(1)*(r^2-1)*(s^2-1)*(q^2-1);
        
        % r derivative        
        dN1dr=1/(8*a)*(2*r-1)*(s^2-s)*(q^2-q);
        dN2dr=1/(8*a)*(2*r+1)*(s^2-s)*(q^2-q);
        dN3dr=1/(8*a)*(2*r+1)*(s^2+s)*(q^2-q);
        dN4dr=1/(8*a)*(2*r-1)*(s^2+s)*(q^2-q);
        
        dN5dr=1/(8*a)*(2*r-1)*(s^2-s)*(q^2+q);
        dN6dr=1/(8*a)*(2*r+1)*(s^2-s)*(q^2+q);
        dN7dr=1/(8*a)*(2*r+1)*(s^2+s)*(q^2+q);
        dN8dr=1/(8*a)*(2*r-1)*(s^2+s)*(q^2+q);
        
        dN9dr=-1/(4*a)*(2*r)*(s^2-s)*(q^2-q);
        dN10dr=-1/(4*a)*(2*r+1)*(s^2-1)*(q^2-q);
        dN11dr=-1/(4*a)*(2*r)*(s^2+s)*(q^2-q);
        dN12dr=-1/(4*a)*(2*r-1)*(s^2-1)*(q^2-q);
        
        dN13dr=-1/(4*a)*(2*r-1)*(s^2-s)*(q^2-1);
        dN14dr=-1/(4*a)*(2*r+1)*(s^2-s)*(q^2-1);
        dN15dr=-1/(4*a)*(2*r+1)*(s^2+s)*(q^2-1);
        dN16dr=-1/(4*a)*(2*r-1)*(s^2+s)*(q^2-1);
        
        dN17dr=-1/(4*a)*(2*r)*(s^2-s)*(q^2+q);
        dN18dr=-1/(4*a)*(2*r+1)*(s^2-1)*(q^2+q);
        dN19dr=-1/(4*a)*(2*r)*(s^2+s)*(q^2+q);
        dN20dr=-1/(4*a)*(2*r-1)*(s^2-1)*(q^2+q);
        
        dN21dr=1/(2*a)*(2*r)*(s^2-1)*(q^2-q);
        dN22dr=1/(2*a)*(2*r)*(s^2-s)*(q^2-1);
        dN23dr=1/(2*a)*(2*r+1)*(s^2-1)*(q^2-1);
        dN24dr=1/(2*a)*(2*r)*(s^2+s)*(q^2-1);
        dN25dr=1/(2*a)*(2*r-1)*(s^2-1)*(q^2-1);
        dN26dr=1/(2*a)*(2*r)*(s^2-1)*(q^2+q);
        
        dN27dr=-1/(1*a)*(2*r)*(s^2-1)*(q^2-1);
        
        % s derivative        
        dN1ds=1/(8*b)*(r^2-r)*(2*s-1)*(q^2-q);
        dN2ds=1/(8*b)*(r^2+r)*(2*s-1)*(q^2-q);
        dN3ds=1/(8*b)*(r^2+r)*(2*s+1)*(q^2-q);
        dN4ds=1/(8*b)*(r^2-r)*(2*s+1)*(q^2-q);
        
        dN5ds=1/(8*b)*(r^2-r)*(2*s-1)*(q^2+q);
        dN6ds=1/(8*b)*(r^2+r)*(2*s-1)*(q^2+q);
        dN7ds=1/(8*b)*(r^2+r)*(2*s+1)*(q^2+q);
        dN8ds=1/(8*b)*(r^2-r)*(2*s+1)*(q^2+q);
        
        dN9ds=-1/(4*b)*(r^2-1)*(2*s-1)*(q^2-q);
        dN10ds=-1/(4*b)*(r^2+r)*(2*s)*(q^2-q);
        dN11ds=-1/(4*b)*(r^2-1)*(2*s+1)*(q^2-q);
        dN12ds=-1/(4*b)*(r^2-r)*(2*s)*(q^2-q);
        
        dN13ds=-1/(4*b)*(r^2-r)*(2*s-1)*(q^2-1);
        dN14ds=-1/(4*b)*(r^2+r)*(2*s-1)*(q^2-1);
        dN15ds=-1/(4*b)*(r^2+r)*(2*s+1)*(q^2-1);
        dN16ds=-1/(4*b)*(r^2-r)*(2*s+1)*(q^2-1);
        
        dN17ds=-1/(4*b)*(r^2-1)*(2*s-1)*(q^2+q);
        dN18ds=-1/(4*b)*(r^2+r)*(2*s)*(q^2+q);
        dN19ds=-1/(4*b)*(r^2-1)*(2*s+1)*(q^2+q);
        dN20ds=-1/(4*b)*(r^2-r)*(2*s)*(q^2+q);
        
        dN21ds=1/(2*b)*(r^2-1)*(2*s)*(q^2-q);
        dN22ds=1/(2*b)*(r^2-1)*(2*s-1)*(q^2-1);
        dN23ds=1/(2*b)*(r^2+r)*(2*s)*(q^2-1);
        dN24ds=1/(2*b)*(r^2-1)*(2*s+1)*(q^2-1);
        dN25ds=1/(2*b)*(r^2-r)*(2*s)*(q^2-1);
        dN26ds=1/(2*b)*(r^2-1)*(2*s)*(q^2+q);
        
        dN27ds=-1/(1*b)*(r^2-1)*(2*s)*(q^2-1);
        
        % q derivative        
        dN1dq=1/(8*c)*(r^2-r)*(s^2-s)*(2*q-1);
        dN2dq=1/(8*c)*(r^2+r)*(s^2-s)*(2*q-1);
        dN3dq=1/(8*c)*(r^2+r)*(s^2+s)*(2*q-1);
        dN4dq=1/(8*c)*(r^2-r)*(s^2+s)*(2*q-1);
        
        dN5dq=1/(8*c)*(r^2-r)*(s^2-s)*(2*q+1);
        dN6dq=1/(8*c)*(r^2+r)*(s^2-s)*(2*q+1);
        dN7dq=1/(8*c)*(r^2+r)*(s^2+s)*(2*q+1);
        dN8dq=1/(8*c)*(r^2-r)*(s^2+s)*(2*q+1);
        
        dN9dq=-1/(4*c)*(r^2-1)*(s^2-s)*(2*q-1);
        dN10dq=-1/(4*c)*(r^2+r)*(s^2-1)*(2*q-1);
        dN11dq=-1/(4*c)*(r^2-1)*(s^2+s)*(2*q-1);
        dN12dq=-1/(4*c)*(r^2-r)*(s^2-1)*(2*q-1);
        
        dN13dq=-1/(4*c)*(r^2-r)*(s^2-s)*(2*q);
        dN14dq=-1/(4*c)*(r^2+r)*(s^2-s)*(2*q);
        dN15dq=-1/(4*c)*(r^2+r)*(s^2+s)*(2*q);
        dN16dq=-1/(4*c)*(r^2-r)*(s^2+s)*(2*q);
        
        dN17dq=-1/(4*c)*(r^2-1)*(s^2-s)*(2*q+1);
        dN18dq=-1/(4*c)*(r^2+r)*(s^2-1)*(2*q+1);
        dN19dq=-1/(4*c)*(r^2-1)*(s^2+s)*(2*q+1);
        dN20dq=-1/(4*c)*(r^2-r)*(s^2-1)*(2*q+1);
        
        dN21dq=1/(2*c)*(r^2-1)*(s^2-1)*(2*q-1);
        dN22dq=1/(2*c)*(r^2-1)*(s^2-s)*(2*q);
        dN23dq=1/(2*c)*(r^2+r)*(s^2-1)*(2*q);
        dN24dq=1/(2*c)*(r^2-1)*(s^2+s)*(2*q);
        dN25dq=1/(2*c)*(r^2-r)*(s^2-1)*(2*q);
        dN26dq=1/(2*c)*(r^2-1)*(s^2-1)*(2*q+1);
        
        dN27dq=-1/(1*c)*(r^2-1)*(s^2-1)*(2*q);
        
        
        % generation of displacement-strain matrix
        
        B0 = [dN1dr 0 0 dN2dr 0 0 dN3dr 0 0 dN4dr 0 0 dN5dr 0 0 dN6dr 0 0 dN7dr 0 0 dN8dr 0 0 dN9dr 0 0 dN10dr 0 0 dN11dr 0 0 dN12dr 0 0 dN13dr 0 0 dN14dr 0 0 dN15dr 0 0 dN16dr 0 0 dN17dr 0 0 dN18dr 0 0 dN19dr 0 0 dN20dr 0 0 dN21dr 0 0 dN22dr 0 0 dN23dr 0 0 dN24dr 0 0 dN25dr 0 0 dN26dr 0 0 dN27dr 0 0
            0 dN1ds 0 0 dN2ds 0 0 dN3ds 0 0 dN4ds 0 0 dN5ds 0 0 dN6ds 0 0 dN7ds 0 0 dN8ds 0 0 dN9ds 0 0 dN10ds 0 0 dN11ds 0 0 dN12ds 0 0 dN13ds 0 0 dN14ds 0 0 dN15ds 0 0 dN16ds 0 0 dN17ds 0 0 dN18ds 0 0 dN19ds 0 0 dN20ds 0 0 dN21ds 0 0 dN22ds 0 0 dN23ds 0 0 dN24ds 0 0 dN25ds 0 0 dN26ds 0 0 dN27ds 0
            0 0 dN1dq 0 0 dN2dq 0 0 dN3dq 0 0 dN4dq 0 0 dN5dq 0 0 dN6dq 0 0 dN7dq 0 0 dN8dq 0 0 dN9dq 0 0 dN10dq 0 0 dN11dq 0 0 dN12dq 0 0 dN13dq 0 0 dN14dq 0 0 dN15dq 0 0 dN16dq 0 0 dN17dq 0 0 dN18dq 0 0 dN19dq 0 0 dN20dq 0 0 dN21dq 0 0 dN22dq 0 0 dN23dq 0 0 dN24dq 0 0 dN25dq 0 0 dN26dq 0 0 dN27dq
            dN1ds dN1dr 0 dN2ds dN2dr 0 dN3ds dN3dr 0 dN4ds dN4dr 0 dN5ds dN5dr 0 dN6ds dN6dr 0 dN7ds dN7dr 0 dN8ds dN8dr 0 dN9ds dN9dr 0 dN10ds dN10dr 0 dN11ds dN11dr 0 dN12ds dN12dr 0 dN13ds dN13dr 0 dN14ds dN14dr 0 dN15ds dN15dr 0 dN16ds dN16dr 0 dN17ds dN17dr 0 dN18ds dN18dr 0 dN19ds dN19dr 0 dN20ds dN20dr 0 dN21ds dN21dr 0 dN22ds dN22dr 0 dN23ds dN23dr 0 dN24ds dN24dr 0 dN25ds dN25dr 0 dN26ds dN26dr 0 dN27ds dN27dr 0
            dN1dq 0 dN1dr dN2dq 0 dN2dr dN3dq 0 dN3dr dN4dq 0 dN4dr dN5dq 0 dN5dr dN6dq 0 dN6dr dN7dq 0 dN7dr dN8dq 0 dN8dr dN9dq 0 dN9dr dN10dq 0 dN10dr dN11dq 0 dN11dr dN12dq 0 dN12dr dN13dq 0 dN13dr dN14dq 0 dN14dr dN15dq 0 dN15dr dN16dq 0 dN16dr dN17dq 0 dN17dr dN18dq 0 dN18dr dN19dq 0 dN19dr dN20dq 0 dN20dr dN21dq 0 dN21dr dN22dq 0 dN22dr dN23dq 0 dN23dr dN24dq 0 dN24dr dN25dq 0 dN25dr dN26dq 0 dN26dr dN27dq 0 dN27dr
            0 dN1dq dN1ds 0 dN2dq dN2ds 0 dN3dq dN3ds 0 dN4dq dN4ds 0 dN5dq dN5ds 0 dN6dq dN6ds 0 dN7dq dN7ds 0 dN8dq dN8ds 0 dN9dq dN9ds 0 dN10dq dN10ds 0 dN11dq dN11ds 0 dN12dq dN12ds 0 dN13dq dN13ds 0 dN14dq dN14ds 0 dN15dq dN15ds 0 dN16dq dN16ds 0 dN17dq dN17ds 0 dN18dq dN18ds 0 dN19dq dN19ds 0 dN20dq dN20ds 0 dN21dq dN21ds 0 dN22dq dN22ds 0 dN23dq dN23ds 0 dN24dq dN24ds 0 dN25dq dN25ds 0 dN26dq dN26ds 0 dN27dq dN27ds];

        
end