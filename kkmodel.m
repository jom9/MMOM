% n is the input index of refraction vector
% k is the input complex index of refraction vector
% lamda is the input wavelengths
% d is the thickness

R=zeros(size(n)); % R is the predicted reflectance
T=zeros(size(n)); % T is the predicted transmitance
d=90*10^-6;
%n=n';
%lamda=lamda';
%k=k';
for i=1:1:size(n)
    R(i)=((n(i)-1)^2+k(i)^2)/( (n(i)+1)^2 +k(i)^2 );
    T(i)=exp(4*pi*k(i)*d/lamda(i))/(1-R(i))^2;
   
    
end
P1=polyfit(lamda,R,1);
P2=polyfit(lamda,R,2);
P3=polyfit(lamda,R,3);

l1=1*10^-6;
l2=10*10^-6;
step= .1*10^-6;


nout1 = zeros(l2-l1/step);
kout1 = zeros(l2-l1/step);
nout2 = zeros(l2-l1/step);
kout2 = zeros(l2-l1/step);
nout3 = zeros(l2-l1/step);
kout3 = zeros(l2-l1/step);


wavelength= zeros(l2-l1/step);
for i=1:1: (l2-l1)/step
    R1 = @(x) P1(1) + P1(2)*x;
    R2 = @(x) P2(1) + P2(2)*x+P2(3)*x.^2;
    R3 = @(x) P3(1) + P3(2)*x.^2+P3(3)*x.^3;
    wavelength(i) = ((i-1)*step+l1);
    fun1 = @(x) log(P1(1) + P1(2)*x)./(x.^2-wavelength(i)^2);
    fun2 = @(x) log(P2(1) + P2(2)*x+P2(3)*x.^2)./(x.^2-wavelength(i)^2);
    fun3 = @(x) log(P3(1) + P3(2)*x.^2+P3(3)*x.^3)./(x.^2-wavelength(i)^2);
    
    phi1 = wavelength(i)/pi * integral(fun1,0,Inf) +pi ;
    phi2 = wavelength(i)/pi * integral(fun2,0,Inf) + pi ;
    phi3 = wavelength(i)/pi * integral(fun3,0,Inf) + pi ;
    
    
    
    nout1(i)= ( 1-R1(wavelength(i)) )/ (1+R1(wavelength(i)) +2*(R1(wavelength(i))*cos(phi1))^.5);
    nout2(i)= ( 1-R2(wavelength(i)) )/ (1+R2(wavelength(i)) +2*(R2(wavelength(i))*cos(phi2))^.5);
    nout3(i)= ( 1-R3(wavelength(i)) )/ (1+R3(wavelength(i)) +2*(R3(wavelength(i))*cos(phi3))^.5);
  
   
    
    kout1(i) =  -2*(R1(wavelength(i))*sin(phi1))^.5/(1+R1(wavelength(i)+2*(R1(wavelength(i))*cos(phi1))^.5));
    kout2(i) =  -2*(R2(wavelength(i))*sin(phi2))^.5/(1+R2(wavelength(i)+2*(R2(wavelength(i))*cos(phi2))^.5));
    kout3(i) =  -2*(R3(wavelength(i))*sin(phi3))^.5/(1+R3(wavelength(i)+2*(R3(wavelength(i))*cos(phi3))^.5));
    
end

f1 = figure;
f2 = figure;
f3= figure;

figure(f1)
hold on
plot(wavelength,nout1)
plot(wavelength,nout2)
plot(wavelength,nout3)

figure(f2)
hold on
plot(wavelength,kout1)
plot(wavelength,kout2)
plot(wavelength,kout3)

Rout1 = ((nout1-1).^2+kout1.^2)./((nout1+1).^2+kout1.^2);
Rout2 = ((nout2-1).^2+kout2.^2)./((nout2+1).^2+kout2.^2);
Rout3 = ((nout3-1).^2+kout3.^2)./((nout3+1).^2+kout3.^2);

figure(f3)
hold on
plot(wavelength,Rout1)
plot(wavelength,Rout2)
plot(wavelength,Rout3)
