clear
clc
datx = csvread("ITS.csv"); %Data with ~10000 generated paths
%datx = csvread("DataO.csv"); %Data original
n=80; %pack
iter=2;

for l=1:iter
disp(l)
for i=2:n+1
    un_val = 2;
    %un_val = 1;
    up_val = 9269; %degradation paths
    %up_val = 14;  %degradation paths
    r=randi([un_val up_val]);
    for j=2:11
        pack(j,i-1)=datx(j,r);
    end
%-------------------------------------------------------------------    
%Temperature generation & acceleration factor
%sigma = 5; TG = 10;
for i=1:n
    p = rand();
    mu = 25;
    sigma = 5;
    Tx(1,i) = norminv(p,mu,sigma);
    Tx(11,i)=Tx(1,i)+10;
    for j=2:10
        Tx(j,i)=Tx(j-1,i)+((Tx(11,i)-Tx(1,i))/10);
    end
end
for i=1:n
    for j=1:11
        AF(j,i)=exp((-1547.5/8.314)*((1/298)-(1/(Tx(j,i)+273))));
    end
end
%-------------------------------------------------------------------
    
tim=[0 17 42 62 102 122 202 442 722 782 1102]; 
    
end
%Apply AF
for j=1:11
    for i=1:n
        pack(j,i)=pack(j,i).*AF(j,i);
        dat=pack;
        nn=11;
        mm=n;
        
        x1 = [1.0e-3 1.5 1.8]';

        for i=2:mm
        tim=[0 17 42 62 102 122 202 442 722 782 1102];
        tm = tim';
        ym = dat(:,i);
    
        opt = optimset ( "TolX",1.e-20, "MaxFunEvals", 100000,"MaxIter", 100000,"DerivativeCheck",'on', "display","iter");
        [x, fval, info, output] = fminunc(@(c)myfun1(c,tm,ym,nn),x1,opt);
        %% PARAMETERs VALUES
        parameter(i-1,:)= x;
        x1=x;
        disp(x)
        disp(i)
        
        end
        
        
        
    end
end
%ax{l}=a;

end

function y = myfun1(x, tm, xm,n)
   a = x(1);
   b = x(2);
   c = x(3);
   ym = zeros(n-1,1);
   for j = 1:n-1
       dT = (tm(j+1))^c - (tm(j))^c;
       dX = (xm(j+1))-(xm(j));
       ym(j)=gampdf(dX,a*dT,b); 
   end
   y = -sum(log(ym));
end 