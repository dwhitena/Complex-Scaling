%% 1d - Complex Scaling Fourier Grid Hamiltonian Method

% Specifically, 1D, 1 particle scattering from a delta like gaussian well
% with exponential steps

% This is also solvable with the multistep potential method for comparison


%% Discretization

clear all

i = sqrt(-1);
N = 1999; % odd, position
a = 12; % boundary of mesh
dx = (2*a)/(N-1); % position spacing
n = (N-1)/2; % even, momentum

xtab = -a:dx:a;

theta = -0.35; % Scaling Angle

% Define Gaussian Potential 
aa = 0;  % height
bb = 9;
b = 0.5;  % width of the gaussian
c = 0;  % sharpness of the step
d = 0;  % boundary of the step
pot = @(x) aa./(1+exp(-2.*c.*(x.*exp(i*theta)+d)))...
    - aa./(1+exp(-2.*c.*(x.*exp(i*theta)-d)))...
    - bb.*exp(-((x.^2).*exp(i*theta*2))./b);
pottab = pot(xtab);


%% Calculate Discrete Hamiltonian

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
    if (real(vals(q,q))<-6)
        norm = 1/(sum(vecs(:,q).^2).*dx);
        wfn = sqrt(norm).*vecs(:,q);
        E = vals(q,q)
    end
end


%% Calculate van Weizsacker Energy

clear wfn

norm = 1/(sum(vecs(:,84).^2).*dx);
wfn = sqrt(norm).*vecs(:,84);

clear tkernel vkernel dendd den t v

den = wfn.^2;

dendd = gradient(den,dx).*gradient(den,dx);
tkernel = (dendd./den);

vkernel = zeros(1,N);

for k = 1:N;
   vkernel(k) = den(k)*pottab(k); 
end
%vkernel = (den).*(pottab');

t = (1/8)*exp(-2*i*theta)*integ_trap_discrete(xtab,tkernel);
%t = (1/8)*exp(-2*i*theta)*sum(tkernel).*dx;
v = integ_trap_discrete(xtab,vkernel);
%v = sum(vkernel).*dx;

energy = t + v


%% Calculate width from asymptotic behavior

eres = vals(240,240);

ek = sqrt(2*eres);

afunc = @(x) exp(-i.*ek.*x.*exp(i*theta));
afunctab = afunc(xtab);

kernel = sqrt(ek).*wfn.*afunctab';
gamma_loc = conj(kernel).*kernel;


%% Convergance

% a = 15; % boundary of mesh

% aa = 4;  % height
% b = 0.5;  % width of the gaussian
% c = 4;  % sharpness of the step
% d = 2;  % boundary of the step

% theta = -0.35;








