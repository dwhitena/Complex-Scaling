%% 2D Finite Difference Method - Gaussian

%% Parameters

clear all

i = sqrt(-1);

l = 5;
N = 599;

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
aa = 4;  % height
b = 0.5;  % width of the gaussian
c = 4;  % sharpness of the step
d = 2;  % boundary of the step
pot = @(x,y) aa/(1+exp(-2*c*(x*exp(i*theta)+d)))...
    - aa/(1+exp(-2*c*(x*exp(i*theta)-d)))...
    - aa*exp(-((x^2)*exp(i*theta*2))/b)...
    + aa/(1+exp(-2*c*(y*exp(i*theta)+d)))...
    - aa/(1+exp(-2*c*(y*exp(i*theta)-d)))...
    - aa*exp(-((y.^2)*exp(i*theta*2))/b);
pottab = zeros(N,N);
for q = 1:N;
    for j = 1:N;
        pottab(q,j) = pot(xtab(q), ytab(j));
    end
end


%% Insert KS potential is using

load('kstestwhole');

pottab = zeros(N,N);

for q = 1:N;
    for j = 1:N;
        pottab(q,j) = ks(q) + ks(j);
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
          H(q,j) = H(q,j) + pottab(q); 
       end
   end
end

[vecs,vals] = eig(H);

for q = 1:N;
    if (imag(vals(q,q))<0 & real(vals(q,q)<2))
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

while (error > 0.000001)


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


%%  Von Weizsacker

clear tkernel vkernel dendd den t v

norm = 2/(sum(sum(psi.*psi)).*h.*h);
psi = sqrt(norm).*psi;

den = sum(psi.*psi).*h;

dendd = gradient(den,h).*gradient(den,h);
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









