%% 2D Finite Difference Method - Gaussian

%% Parameters

clear all

tic;

i = sqrt(-1);

l = 12;
N = 1999;

h = (2*l)/(N-1);

xtab = -l:h:l;
ytab = -l:h:l;

psi = zeros(N,N);
psi(:,:) = 1;

% Boundary conditions
psi(:,1) = 0;
psi(:,N) = 0;
psi(1,:) = 0;
psi(N,:) = 0;

E = 3.25 - 1*i;

% Scaling Parameter
theta = -0.35;

% Discretize Potential 
aa = 0;  % height
bb = 9; % depth of the gaussian
b = 0.5;  % width of the gaussian
c = 0;  % sharpness of the step
d = 0;  % boundary of the step
lambda = 1;
pot = @(x,y) aa/(1+exp(-2*c*(x*exp(i*theta)+d)))...
    - aa/(1+exp(-2*c*(x*exp(i*theta)-d)))...
    - bb*exp(-((x.^2)*exp(i*theta*2))/b)...
    + aa/(1+exp(-2*c*(y*exp(i*theta)+d)))...
    - aa/(1+exp(-2*c*(y*exp(i*theta)-d)))...
    - bb*exp(-((y.^2)*exp(i*theta*2))/b)...
    + lambda/sqrt(1+(x*exp(i*theta)-y*exp(i*theta)).^2);
pottab = zeros(N,N);
for q = 1:N;
    for j = 1:N;
        pottab(q,j) = pot(xtab(q), ytab(j));
    end
end

potnon = @(x,y) aa/(1+exp(-2*c*(x*exp(i*theta)+d)))...
    - aa/(1+exp(-2*c*(x*exp(i*theta)-d)))...
    - bb*exp(-((x^2)*exp(i*theta*2))/b)...
    + aa/(1+exp(-2*c*(y*exp(i*theta)+d)))...
    - aa/(1+exp(-2*c*(y*exp(i*theta)-d)))...
    - bb*exp(-((y.^2)*exp(i*theta*2))/b);
potnontab = zeros(N,N);
for q = 1:N;
    for j = 1:N;
        potnontab(q,j) = potnon(xtab(q), ytab(j));
    end
end



%% Compute better initial wavefunction

clear psi


n = (N-1)/2; % even, momentum

H = zeros(N,N);

for q = 1:N;
   for j = 1:N;
       term1 = 0;
       for l = 1:n;
           term1 = term1 + cos((l*2*pi*(q-j))/N)*2*(((pi*l)/(N*h))^2);
       end
       H(q,j) =  ((2*exp(-2*i*theta))/N)*term1;
       if (q==j)
          H(q,j) = H(q,j) + potnontab(q); 
       end
   end
end

[vecs,vals] = eig(H);

%%

for q = 1:N;
    if (real(vals(q,q))<-6)
        vals(q,q)
    end
    if (real(vals(q,q))<-6)
        norm = 1/(sum(vecs(:,q).^2).*h);
        wfn = sqrt(norm).*vecs(:,q);
        E = 2*vals(q,q);
    end
end


%% Construct 2D wavefunction

psi = zeros(N,N);

for j = 1:N;
    for k = 1:N;
        psi(j,k) = wfn(j)*wfn(k);
    end
end


%% Apply Correction Formula

error = 100;

while (error > 0.00000005)


for q = 2:(N-1);
    for j = 2:(N-1);
        if (q==(N/2) & j==(N/2))
            Enew = ((exp(-2*i*theta)*(4*psi(q,j) - psi(q-1,j) - psi(q+1,j) - psi(q,j-1)...
                - psi(q,j+1)))/(2*(h^2)*psi(q,j))) + pottab(q,j);
        elseif (q==((N-1)/2) & j==((N-1)/2))
            Enew = ((exp(-2*i*theta)*(4*psi(q,j) - psi(q-1,j) - psi(q+1,j) - psi(q,j-1)...
                - psi(q,j+1)))/(2*(h^2)*psi(q,j))) + pottab(q,j);
        else
        psi(q,j) = -(psi(q-1,j) + psi(q+1,j) + psi(q,j-1) + psi(q,j+1))/...
            (-4 + 2*exp(2*i*theta)*(h^2)*(E - pottab(q,j)));
        end
    end
end

error = (abs(E - Enew)/abs(Enew))*100;

%if (real(Enew)<real(E))
    E = Enew
%end

end

E


%% Functional


clear tkernel vkernel dendd den t v

norm = 2/(sum(sum(psi.*psi)).*h.*h);
psi = sqrt(norm).*psi;

den = sum(psi.*psi).*h;

%%

dendd = gradient(den,h).*gradient(den,h);
tkernel = (dendd./den);

vkernel = zeros(1,N);

for k = 1:N;
   vkernel(k) = den(k)*potnontab(k); 
end
%vkernel = (den).*(pottab');

t = (1/8)*exp(-2*i*theta)*integ_trap_discrete(xtab,tkernel);
%t = (1/8)*exp(-2*i*theta)*sum(tkernel).*dx;
v = integ_trap_discrete(xtab,vkernel);
%v = sum(vkernel).*dx;

potee = @(x,y) lambda/sqrt(1+(x*exp(i*theta)-y*exp(i*theta)).^2);
poteetab = zeros(N,N);
for q = 1:N;
    for j = 1:N;
        poteetab(q,j) = potee(xtab(q), ytab(j));
    end
end

hkernel = zeros(N,N);
for q = 1:N;
    for j = 1:N;
        hkernel(q,j) = den(q)*den(j)*poteetab(q,j);
    end
end

hartree = sum(sum(hkernel)).*h.*h;

energy = t + v + hartree - (1/2)*hartree

toc;

%% Exchange only correction


clear lambda potee poteetab

lambda = 1;
potee = @(x,y) lambda/sqrt(1+(x*exp(i*theta)-y*exp(i*theta)).^2);
poteetab = zeros(N,N);
for q = 1:N;
    for j = 1:N;
        poteetab(q,j) = potee(xtab(q), ytab(j));
    end
end

excor = sum(sum(psi.*poteetab.*psi)).*h.*h;

energyex = energy + excor


%% Local Correlation (contact interacting, Phys Rev A 70, 032508 (2004))

para = (-1/24);
parb = -0.00436143;
pard = 0.252758;
pare = 0.0174457;

kerntop = (para.*(den.^3) + parb.*(den.^2));
kernbot = ((den.^2) + pard.*den + pare);

energyc = sum(kerntop./kernbot).*h

energywithc = energy + energyc

%% Compute Kohn-Sham potential

ks = zeros(1,N);
sqrtdendd = zeros(1,N);

sqrtdendd = gradient(gradient(sqrt(den),h),h);

ks = exp(-2*i*theta).*(sqrtdendd./(2.*sqrt(den))) + (3.256527226273844 - 0.006844587724273i);
%ks = exp(-2*i*theta).*(sqrtdendd./(2.*sqrt(den)));

%%

for q = 1:N
   if (real(ks(q)) < 0)
       ks(q) = 0;
       q
   end
end

%% Solve Kohn-Shame system

clear H q j l wfn pottab

pottab = kstest;

H = zeros(N,N);

n = (N-1)/2; % even, momentum

for q = 1:N;
   for j = 1:N;
       term1 = 0;
       for l = 1:n;
           term1 = term1 + cos((l*2*pi*(q-j))/N)*2*(((pi*l)/(N*h))^2);
       end
       H(q,j) =  ((2*exp(-2*i*theta))/N)*term1;
       if (q==j)
          H(q,j) = H(q,j) + pottab(q); 
       end
   end
end

[vecs,vals] = eig(H);

for q = 1:N;
    if (imag(vals(q,q))<0 & real(vals(q,q)<4))
        norm = 1/(sum(vecs(:,q).^2).*h);
        wfn = sqrt(norm).*vecs(:,q);
        E = vals(q,q);
    end
end


% Construct the density

den = 2.*(wfn.*wfn);







