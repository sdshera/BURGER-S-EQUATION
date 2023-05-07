%% Srijan Dasgupta Assignment 3
clear all 
close all 
clc 
%% Problem Definition 
tolerance=1e-6;                                     %tolerance 
N=60;                                              %max modes number  
k = 1:1:N;                                          %range of modes 
Re=40;                                              %reynold's number 
v=1/Re;                                             %kinematic viscosity  
dt = 0.025*(Re/N^2);                                %timestep (CFL)
Uk = 1./k;                                          %dns initial fourier coefficient 
Uk_les=1./k;                                        %les initial fourier coefficient 

%% LES viscosity correction 
Ck=0.4523;                                          %kolmogorov constant for 1D burger                          
m=2;                                                %slope of energy spectrum 
vt_inf=0.31*((5-m)/(m+1))*sqrt(3-m)*Ck^(-3/2);      

%% DNS Solver 
error=1;                                            %error to initiate problem 
iter=0;                                             %iteration counter 
t=0;                                                %time counter 
while error>tolerance && iter<100000
    Uk(1)=1;                                        %energy input 
    Uk_0=Uk;                                        %original value 
    for j=k                                         %modes loop 
        C=conv(j,N,Uk_0);                           %convective term 
        D=((j^2*Uk_0(j))*v);                       %diffusive term
        Uk(j) = Uk_0(j) - dt*(C+D);                 %solver 
    end 
    Uk(1)=1;
    error=max(abs((Uk-Uk_0)/dt));                   %error 
    t=t+dt;
    iter=iter+1;
    
end 

E_DNS = Uk.*conj(Uk);                               %DNS energy spectrum 

%% LES solver 
error=1;                                            %error to initiate problem 
iter=0;                                             %iteration counter
t=0;                                                %time counter 
while error>tolerance && iter<100000
    Uk_les(1)=1;                                    %energy input 
    Uk_0=Uk_les;                                    %original value
    for j=k                                         %modes loop
        C=conv(j,N,Uk_0);                           %convective term
        
        vt_star=1+34.5*exp((-3.03*k(end))/j);       %non-dimensional eddy viscosity 
        E_loop=conj(Uk_0(j))*Uk_0(j);               
        vt=vt_inf*sqrt(E_loop/N)*vt_star;           %eddy viscosity 
        v_eff=(v)+(vt);                             %effective LES viscosity 
        
        D=((j^2*Uk_0(j))*v_eff);                    %DNS diffusive term
        Uk_les(j) = Uk_0(j) - dt*(C+D);
    end 
    Uk_les(1)=1;
    error=max(abs((Uk_les-Uk_0)/dt));               %error 
    t=t+dt;
    iter=iter+1;
end 

E_LES = Uk_les.*conj(Uk_les);                       %LES energy spectrum 

%% Plot DNS Solution 
figure (1)
loglog(k,E_DNS)                                 
title ('Energy Spetrum vs Wave number (DNS)')
xlabel('Wave number (k)')
ylabel('Turbulence energy spectrum (E)')
hold on 
%% Plot LES Solution  
loglog(k,E_LES)
title ('Energy Spetrum vs Wave number')
xlabel('Wave number (k)')
ylabel('Turbulence energy spectrum (E)')
DNS= sprintf('DNS N=%.1f Re=%.1f',N,Re);
LES = sprintf('LES N=%.1f Re=%.1f',N,Re);
legend (DNS,LES)
hold off 
%% Convective term function  
function sum=conv(k,N,U_f)                          % inputs 
sum=0;                                              %initiate sum 
for j=k                                             %for k mode 
    for p=(-N+j):N                                  %fourier series truncated for abs(k)=<N
        q=j-p;                                      
    
        if q~= 0 & p~= 0                            
            if q < 0                                %conditions for q          
            u_q = conj(U_f(abs(q)));
            else                        
            u_q = U_f(q);
            end
        
            if p < 0                                %conditions for p 
            u_p = conj(U_f(abs(p)));
            else
            u_p = U_f(p);
            end 
            sum = sum +  u_p * q * 1i * u_q;        %total convective term sum 
        end 
        
    end 
end 
end 
