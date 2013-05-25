%% Program to compute Tranmission vs. Energy using the Multistep
%% Potential Approximation (see J. Appl. Phys. 61 (4), 1987)

% Variables:

%   x - position
%   eng - energy

% Functions used in this program:

%   pot - potential well function
%   integ_trap - trapazoidal integration
%   bisection_noshow - rootfinding


%% Form Discrete Potential Steps

clear all

% Define Potential

i = sqrt(-1);
l = 0.1;
j = 0;
pot = @(x) (0.5*x^2 - j)*exp(-l*x^2) + j; %Moiseyev
xl = -8;
xr = 8;

% Form discrete quantum well potential:

M = 25;  % number of step potentials
vtab = zeros(1,M);  % values of steps

h = abs(xl - xr)/M;
xtab = zeros(1,M+1);  % discrete spacial coordinate

xtab(1) = xl;
for q = 2:(M+1)
    xtab(q) = xtab(q-1) + h; 
end

vtab(1) = j;
vtab(M) = j;
for q = 2:M-1
   vtab(q) = pot((xtab(q) + xtab(q+1))/2);
end


%% Calculate needed coefficients and Reflection

k = @(eng,n) sqrt(2*(eng-vtab(n)));

r = @(eng,m) (k(eng,m-1)-k(eng,m))/(k(eng,m-1)+k(eng,m));

ro = @(eng) r(eng,M);

for q = M:-1:3
    ro = @(eng) (r(eng,q-1) + ro(eng)*exp(2*i*k(eng,q-1)*xtab(q)))/...
        (1 + r(eng,q-1)*ro(eng)*exp(2*i*k(eng,q-1)*xtab(q)));
end

reflection = @(eng) abs(ro(eng))^2;
transmission = @(eng) 1 - reflection(eng);

plotfunc(transmission,j,4,200)

%% Discrete Reflection / Transmission

N = 100;
Rtab = zeros(1,N);
Ttab = zeros(1,N);
engtab = zeros(1,N);

eh = (10 - j)/N;
engtab = j:eh:10;

for q = 1:N
    Rtab(q) = reflection(engtab(q));
end
Ttab = 1 - Rtab;

figure
plot(engtab,Rtab)



%%

for q = 1:N;
   Rtab(q) = reflection(engtab(q));
   Ttab(q) = 1 - Rtab(q);
end





