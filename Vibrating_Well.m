
E_vals = [];
o = 0;
T_vals = [];


%for Ampl = 0:0.1:100
    
    
o = o + 1;





% Creating Array of Frequencies
W = [];
k = 0;
%for r = 1:30
%    for z = 0.25:0.25:1
%        k = k + 1;
 %       W(k) = z*10^r;
  %  end
%end


%for w = 1:length(W)

%i = i + 1;

%Defining variables
h = 6.64e-34;
nm = 1e-9;
eV = 1.60218e-19;
hbar = 1.055e-34;
m = 9.11e-31; % rest mass of the electron


% Defining spatal and temporal steps as well as grid
Nx = 2001;
Nt = 15000;
L = 150e-9; % (m)
dx = L/Nx;
C1 = 0.0005;
dt = ((2*m*dx^2)/hbar)*C1;

%limit_dt = hbar/(((2*hbar^2)/(m*dx^2))+max(abs(U_const))*(eV));


%Defining Grid Parameters
s = L/150;
x = dx*((-Nx/2):1:(Nx/2)-1);




% Specifying Initial Wave Packet

E = 8; % Energy of electron determines wavelength below
wL = sqrt(h^2/(2*m*(eV)*E));
xc = -2e-8;        % dx*(-Nx/4), initial position of wavepacket.


yR = [];  % Real Wavefunction
yI = [];  % Imaginary Wavefunction

for nx = 1 : Nx
 yR(nx) = exp(-0.5*((x(nx)-xc)/s)^2)*cos((2*pi*(x(nx)-xc))/wL);
 
 
 yI(nx) = exp(-0.5*((x(nx)-xc)/s)^2)*sin((2*pi*(x(nx)-xc))/wL);
 
end



% Normalising Intitial Gaussian Wavefunction

y2 = yR.^2 + yI.^2;
A = simpsons(y2,0,L);

yR = yR./sqrt(A); % Normalised Initial Real Wavefunction
yI = yI./sqrt(A); % Normalised Initial Imaginary Wavefunction
yR0 = yR;
yI0 = yI;

y3 = yR.^2+yI.^2;
Normalisation_Test = simpsons(y3,0,L);  % Correctly Normalised if == 1

%plot(x,yR)
%hold on
%plot(x,yI)

% Defining the Vibrating Potential Well Array

time = 700000;

ts = 1:time;
U = zeros(length(ts),2001);


for t = 1:time

U_upper = 10;
U_lower = 0;
qWell_wide = 0.2*(nm); % Must be an even number
CentreWellSep = 2.6*(nm);%1.8*(nm); % Must be even too



Centre_well_1 = -1e-8 + 0.5*sin(2*pi*0.00001*t)*CentreWellSep; %was 4*pi*10^12
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
U(t,x>RW2) = U_lower;
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


% UPDATE EQUATIONS

C2 = eV*dt/hbar;


% Boundary Conditions Imposed
yR(1) = 0;
yR(end) = 0;
yI(1) = 0;
yI(end) = 0;
yR0(1) = 0;
yR0(end) = 0;
yI0(1) = 0;
yI0(end) = 0;



for nt = 1 : time  % run-time of simulation (will cause issues if packet reaches boundaries).
    
 for nx = 2 : Nx - 1
     yR(nx) = yR(nx) - C1*(yI(nx+1)-2*yI(nx)+yI(nx-1)) + C2*U(nt,nx)*yI(nx);
 end
 
 for nx = 2 : Nx-1
     yI(nx) = yI(nx) + C1*(yR(nx+1)-2*yR(nx)+yR(nx-1)) - C2*U(nt,nx)*yR(nx);
 end

 if rem(nt,1000) == 0
     plot(x,yR+3e4)
     plot(x,yR.^2+yI.^2)
     xlabel('x position (10e-8)')
     ylabel('Probability Density')
     hold on
     plot(x,8e3*U(nt,:))
     plot(x,1e8*U(nt,:))
     hold off
     drawnow 
 end
end



Prob  = yR.^2 + yI.^2;
Norm = sum(Prob);

yR = yR./sqrt(Norm); % Normalised Real Wavefunction
yI = yI./sqrt(Norm); % Normalised Imaginary Wavefunction


Prob = yR.^2 + yI.^2;

Transmission = 0;

Transmission = sum(Prob(1000:2001));


%T_vals(o) = Transmission;

%Ampl_vals(o) = Ampl;




%plot(Ampl_vals,T_vals)
%title('Transmission for Multi Barrier System Oscillating at Constant Frequency For Increasing Oscillation Amplitude (5eV Energy Wavepacket, 10^8 Hz )')
%xlabel('Increasing Amplitude of Oscillation')
%ylabel('Transmission of Wavepacket through Array')
%hold on



