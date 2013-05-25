%% Self Consistent Kohn Sham resonance calculation

% Based on previous model_interacting.m and fgh1dcs.m tests

%% Definitions

clear all

tic;

i = sqrt(-1);
N = 299; % odd, position
a = 7; % boundary of mesh
dx = (2*a)/(N-1); % position spacing
n = (N-1)/2; % even, momentum

xtab = -a:dx:a;

theta = 0.403; % Scaling Angle

% Define Gaussian External Potential 
potextab = zeros(1,N);
pottab = zeros(1,N);
aa = 0.75;  % height
bb = 5;
b = 0.05;  % width of the gaussian
c = 4;  % sharpness of the step
d = 5;  % boundary of the step
pot = @(x) aa./(1+exp(-2.*c.*(x.*exp(i*theta)+d)))...
    - aa./(1+exp(-2.*c.*(x.*exp(i*theta)-d)))...
    - bb.*exp(-((x.^2).*exp(i*theta*2))./b);
pottab = pot(xtab);
potextab = pottab;

% Interaction
lambda = 1;
potee = @(x,y) lambda/sqrt(1+(x*exp(i*theta)-y*exp(i*theta)).^2);
poteetab = zeros(N,N);
for q = 1:N;
    for j = 1:N;
        poteetab(q,j) = potee(xtab(q), xtab(j));
    end
end



%% Initial guess for Kohn Sham orbital

H = zeros(N,N);

for q = 1:N;
   for j = 1:N;
       term1 = 0;
       for l = 1:n;
           term1 = term1 + cos((l*2*pi*(q-j))/N)*2*(((pi*l)/(N*dx))^2);
       end
       H(q,j) =  ((2*exp(-2*i*theta))/N)*term1;
       if (q==j)
          H(q,j) = H(q,j) + pottab(q); 
       end
   end
end

[vecs,vals] = eig(H);

% sort by size of imag part of the energy
for q = 1:N;
    if (imag(vals(q,q))<0)
        vals(q,q)
    end
    if (imag(vals(q,q))<0 & real(vals(q,q)<4))
        norm = 1/(sum(vecs(:,q).^2).*dx);
        wfn = sqrt(norm).*vecs(:,q);
        E = vals(q,q);
    end
end


%% SCF



for m = 1:5;

% Compute Kohn Sham potential:

den = wfn.*wfn;

% Hartree and exchange for KS potential
pothxtab = zeros(1,N);
pothxkernel = zeros(N,N);
for j=1:N;
    for q=1:N;
        pothxkernel(j,q) = den(q).*poteetab(j,q);
    end
end
for j=1:N;
    pothxtab(j) = sum(pothxkernel(j,:)).*dx;
end

% Total
kspot = zeros(1,N);
kspot = potextab + pothxtab -(1/2).*pothxtab;
pottab = kspot;

% Solve KS equation
H = zeros(N,N);

for q = 1:N;
   for j = 1:N;
       term1 = 0;
       for l = 1:n;
           term1 = term1 + cos((l*2*pi*(q-j))/N)*2*(((pi*l)/(N*dx))^2);
       end
       H(q,j) =  ((2*exp(-2*i*theta))/N)*term1;
       if (q==j)
          H(q,j) = H(q,j) + pottab(q); 
       end
   end
end

[vecs,vals] = eig(H);

for q = 1:N;
    if (real(vals(q,q))<0)
        vals(q,q)
    end
    if (real(vals(q,q))<0)
        norm = 1/(sum(vecs(:,q).^2).*dx);
        wfn = sqrt(norm).*vecs(:,q);
        E = vals(q,q);
    end
end

end

%% Calculate energy and lifetime

% Construct the density

den = 2.*(wfn.*wfn);

dendd = gradient(den,dx).*gradient(den,dx);
tkernel = (dendd./den);

vkernel = zeros(1,N);

for k = 1:N;
   vkernel(k) = den(k)*potextab(k); 
end
%vkernel = (den).*(pottab');

t = (1/8)*exp(-2*i*theta)*integ_trap_discrete(xtab,tkernel);
%t = (1/8)*exp(-2*i*theta)*sum(tkernel).*dx;
v = integ_trap_discrete(xtab,vkernel);
%v = sum(vkernel).*dx;

hkernel = zeros(N,N);
for q = 1:N;
    for j = 1:N;
        hkernel(q,j) = den(q)*den(j)*poteetab(q,j);
    end
end

hartree = sum(sum(hkernel)).*dx.*dx;

energy = t + v + hartree - (1/2)*hartree

toc;



%% Kohn Sham


% Hartree and exchange
pothxtab = zeros(1,N);
pothxkernel = zeros(N,N);
for j=1:N;
    for q=1:N;
        pothxkernel(j,q) = den(q).*poteetab(j,q);
    end
end
for j=1:N;
    pothxtab(j) = sum(pothxkernel(j,:)).*dx;
end

% Total
kspot = zeros(1,N);
kspot = potextab + pothxtab -(1/2).*pothxtab;
pottab = kspot;












