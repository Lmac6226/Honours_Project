% Setup of Grid and Initial Wavepacket

%Defining variables
h = 6.64e-34;
nm = 1e-9;
eV = 1.60218e-19;
hbar = 1.055e-34;
m_e = 9.11e-31; % rest mass of the electron
m_p = 1.67262192e-27;


% Defining spatal and temporal steps as well as grid
Nx = 2001;
Nt = 150000;
L = 150e-9; % (m)
dx = L/Nx;
C1 = 0.005;
dt = ((2*m_e*dx^2)/hbar)*C1;
C2 = eV*dt/hbar;

x = dx*((-Nx/2):1:(Nx/2)-1); % vector of grid points




% Specifying Initial Guassian Wave Packet
s = L/150; % Width of gaussian
E = 20; % Energy of electron determines wavelength below
wL = sqrt(h^2/(2*m_e*(eV)*E));
xc = 1e-8;        % dx*(-Nx/4), initial position of wavepacket.

yR = [];  % Real Wavefunction
yI = [];  % Imaginary Wavefunction

for nx = 1 : Nx
 yR(nx) = exp(-0.5*((x(nx)-xc)/s)^2)*cos((2*pi*(x(nx)-xc))/wL);
 
 
 yI(nx) = exp(-0.5*((x(nx)-xc)/s)^2)*sin((2*pi*(x(nx)-xc))/wL);
end

Y = yR + 1i*yI;
Y_conjugate = conj(Y);
YY_initial = Y.*Y_conjugate;


% Normalising Intitial Gaussian Wavefunction

Y = Y./sum(Y.*Y_conjugate);
%Normalisation_Test = simpsons(Y,0,L);  % Correctly Normalised if == 1



% Defining the Vibrating Potential Array


ts = 1:Nt;
U = zeros(length(ts),2001);


for t = 1:Nt

U_upper = 10*eV;
U_lower = 0;
qWell_wide = 2*(nm); % Must be an even number
CentreWellSep = 1.2*(nm);%1.8*(nm); % Must be even too



Centre_well_1 = 0; % + 0.5*nm*sin(2*pi*0.00001*t); Oscillating Barrier
Centre_well_2 = Centre_well_1 + CentreWellSep;
Centre_well_3 =  Centre_well_2 + CentreWellSep;
Centre_well_4 =  Centre_well_3 + CentreWellSep;
Centre_well_5 =  Centre_well_4 + CentreWellSep;
Centre_well_6 =  Centre_well_5 + CentreWellSep;
Centre_well_7 =  Centre_well_6 + CentreWellSep;
Centre_well_8 =  Centre_well_7 + CentreWellSep; 
Centre_well_9 =  Centre_well_8 + CentreWellSep;
Centre_well_10 =  Centre_well_9 + CentreWellSep;

    
RW1 = Centre_well_1+qWell_wide*0.5;
RW2 = Centre_well_2+qWell_wide*0.5;
RW3 = Centre_well_3+qWell_wide*0.5;
RW4 = Centre_well_4+qWell_wide*0.5;
RW5 = Centre_well_5+qWell_wide*0.5;
RW6 = Centre_well_6+qWell_wide*0.5;
RW7 = Centre_well_7+qWell_wide*0.5;
RW8 = Centre_well_8+qWell_wide*0.5; 
RW9 = Centre_well_9+qWell_wide*0.5;
RW10 = Centre_well_10+qWell_wide*0.5;

LW1 = Centre_well_1-qWell_wide*0.5;
LW2 = Centre_well_2-qWell_wide*0.5;
LW3 = Centre_well_3-qWell_wide*0.5;
LW4 = Centre_well_4-qWell_wide*0.5;
LW5 = Centre_well_5-qWell_wide*0.5;
LW6 = Centre_well_6-qWell_wide*0.5;
LW7 = Centre_well_7-qWell_wide*0.5;
LW8 = Centre_well_8-qWell_wide*0.5;
LW9 = Centre_well_9-qWell_wide*0.5;
LW10 = Centre_well_10-qWell_wide*0.5;


U(t,x<LW1) = U_lower;
U(t,x>=LW1 & x<=RW1) = U_upper;
U(t,x>RW1) = U_lower;
%U(t,x>RW1 & x<LW2) = U_lower;
%U(t,x>=LW2 & x<=RW2) = U_upper;
%U(t,x>RW2 & x<LW3) = U_lower;
%U(t,x>=LW3 & x<=RW3) = U_upper;
%U(t,x>RW3 & x<LW4) = U_lower;
%U(t,x>=LW4 & x<=RW4) = U_upper;
%U(t,x>RW4 & x<LW5) = U_lower;
%U(t,x>=LW5 & x<=RW5) = U_upper;
%U(t,x>RW5 & x<LW6) = U_lower;
%U(t,x>=LW6 & x<=RW6) = U_upper;
%U(t,x>RW6 & x<LW7) = U_lower;
%U(t,x>=LW7 & x<=RW7) = U_upper;
%U(t,x>RW7 & x<LW8) = U_lower;
%U(t,x>=LW8 & x<=RW8) = U_upper;
%U(t,x>RW8 & x<LW9) = U_lower;
%U(t,x>=LW9 & x<=RW9) = U_upper;
%U(t,x>RW9 & x<LW10) = U_lower;
%U(t,x>=LW10 & x<=RW10) = U_upper;
%U(t,x>RW10) = U_lower;


end


% Implementing Crank-Nicolsen Matrix Method



% Constructing Matrices

alpha = -1*1i*hbar*dt/(dx^2*m_e);
beta = 2*1i*dt/hbar;

diagonals = [2*(1+alpha)*ones(Nx,1),-alpha*ones(Nx,2)];
A = spdiags(diagonals,[0 -1 1],Nx,Nx);

I = speye(Nx); % Implement Boundary Conditions
A([1 Nx],:) = I([1 Nx],:);

B = [];
B = ones(Nx,1);
B(1) = 0;
B(end) = 0;


% Updating Wavefunction and Plotting

for m = 1:Nt
    
    for j = 2:Nx-1
    B(j) = 2*(1-alpha+(beta/2)*U(1,j))*Y(j)+alpha*Y(j+1)+alpha*Y(j-1);
    end
   
   
    Y = A\B;
    Y_conjugated = conj(Y);

  
    

    if rem(m,10)==0
        plot(x,Y_normalized)
        drawnow
    end

end