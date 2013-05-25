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


ntab = [99 109 149 175 199 299 399 499 599 699 799 899 999 1099 1199 1299];

retabvw = [1.596005833351029 1.601277459516470 1.613467163841009 1.617522608407151 1.619983892164033 1.624725605719650 ...
    1.626394039290449 1.627166502295570 1.627585891290257 ...
    1.627838609486579 1.628002536327466 1.628114865547914 ...
    1.628195177747528 1.628254576502269 1.628299738807182 ...
    1.628334875284116];

imtabvw = [-0.000429579075649 0.001496051983559 0.002212205468666 0.001411173206091 0.000678813913672 -0.001256891191870 ...
    -0.00210259579038 -0.002523374913724 -0.002759607892665 ...
    -0.00290460951986 -0.002999733014362 -0.003065400933010 ...
    -0.003112593848763 -0.003147627629109 -0.003174338650379 ...
    -0.003195164167098];

retab = [1.628514450298723 1.628518136817797 1.628539053170347 1.628537103027341 1.628536916614403 1.628536984690307 ...
    1.628536984587789 1.628536984569074 1.628536984555788 ...
    1.628536984547428 1.628536984539298 1.628536984534369 ...
    1.628536984528817 1.628536984523802 1.628536984523630 ...
    1.628536984519818];

imtab = [-0.003363022998107 -0.003291821247477 -0.003315012136485 -0.003316080111981 -0.003315665170612 -0.003315707068789 ...
    -0.003315707102066 -0.003315707076052 -0.003315707060193 ...
    -0.003315707047948 -0.003315707040241 -0.003315707034637 ...
    -0.003315707028926 -0.003315707019963 -0.003315707027171 ...
    -0.003315707015061];


etabvw = [1.6010 1.60280 1.61001 1.61331 1.6155 1.6203 1.6221 1.6230 1.6235 1.6237 1.6239...
    1.6241, 1.62416 1.62423 1.62428 1.62432];

etab = [1.6246 1.6246 1.6246 1.6246 1.6246 1.6246 1.6246 1.6246 1.6246 1.6246 1.6246...
    1.6246 1.6246 1.6246 1.6246 1.6246];


errre = abs(retabvw - retab)./abs(retab);
errim = abs(imtabvw - imtab)./abs(imtab);
err = abs(etabvw - etab)./abs(etab);





