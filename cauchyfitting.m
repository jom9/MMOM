

% Modeling index of refraction using Cauchy relation
% n^2 = sum(ai/l^2i)
% Note this will only get you valid values outside absorption band
% n and l are the imported values for index of refraction and wavelength
l=lamda;

fit = createFit1(l,n); % creates fit
coeff = coeffvalues(fit); % stores fit coeffeiecents to vector

l1=.5; 
l2=20;
step=.001;
% This is the range of lamda values that will be fitted to

lout = transpose([l1:step:l2]);
nfitted=transpose(zeros(size(lout)));
for i=1:1:size(lout)
    fprintf("%d",coeff(1)+coeff(2)/lout(i)^2 + coeff(3)/lout(i)^4);
    nfitted(i)= coeff(1)+coeff(2)/lout(i)^2 + coeff(3)/lout(i)^4;
end
plot(l,n)
hold on
plot(lout,nfitted)
    