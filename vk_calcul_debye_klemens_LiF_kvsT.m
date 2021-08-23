% This is code V.1. for numerical calculation of thermal conductivity of 
% SHI irradiated non-metals
clear 

% =========================== parameters ==================================
kB = 1.38064852e-23;            % Boltzman constant in J/K
hr = 6.62607015e-34 / (2 * pi); % hbar, reduced Plack constant in J*s
TD = 735;                       % Debye temperature of LiF crystal taken 
                                % from "Atomic vibrations" chapter of 
                                % "Quantum theory of the solid state" by Kantorobich
M = 25.939;                     % Molecular mass LiF in g/mol
MF = 18.998;
MLi = 6.941;

g = 1.5;                        % Gruneisen constant of LiF, a measure of anharmonicity
a = 4.03e-10;                   % a^3 atomic volume in m
d = 1;                          % dimensions of the sample in m 

T = [1 : 100 110 : 10 : 1000];              % temperature
V0 = 4/3 * pi * (0.5 * a)^3;                % atomic volume
C11=11.12e11;C12=4.2e11;C44=6.49e11;        % elasticity constants at 300K
                                            % based on Briscoe and Squire, 1957
vl = 6602;                                  % in m/s based on Briscoe and Squire, 1957
vt = 4938;
%v=(1./(1./vl.^3./3+2./vt.^3./3))^(1/3);     % summation of velocity
v=(1./(1./vl.^2./2+1./vt.^2./2))^(1/2);    

wD=(6*pi.^2./V0)^(1/3).* v;         % Debye frequency
wm = wD;
ww=(1:50)/50;                       % setting frequency range
dw=wm/50;w=wm*ww;      

% 
% wm=wD;                                % maximum frequency set equal to Debye
% ww=(1:50)/50;                         % setting frequency range
% dwm=wm/50;w=wm*ww;      



% load DFT simulation data
k_num = load('numdata.txt');
k_num_4A = k_num(:,1);
k_num_7A = k_num(:,2);
k_num_8A = k_num(:,3);


% load Berman experimental simulation data
k_Ber = load('LiF_k_Berman.txt');
T_Ber = k_Ber(:,1);
k_Ber = k_Ber(:,2)* 1e+2;

% load Pohl experimental simulation data
k_Pohl = load('LiF_k_Pohl.txt');
T_Pohl = k_Pohl(:,1);
k_Pohl = k_Pohl(:,2)* 1e+2;

%T = T_Ber';
% =========================== calculation =================================
% kL - thermal conductivity of the lattice
c1 = kB / (2 * pi^2 * v);       % constant coefficients 
c2 = (kB / hr)^3;

x = hr * w' ./ (kB * T);         % dimensionless parameter


% =========================================================================
% FITTING PARAMETERS

% estimated
%dM = 1.8*M * 1e-1;                % temporary value
dM = 0.9*MF;

c = 1.55e-6;                       % the point defect concentration (per atom)

% coefficients
% v.1.  (Tritt 2004)
%a0 = 1.1e+28;                      
a0 = 1.175e+28;                      
%alpha = 3.2;
alpha = 3.2;
%==========================================================================




% ------------------phonon scattering relaxation time ----------------- 

% Umklapp scattering
% v.1. according to Tritt
tU = (a0 * hr * g^2 * w'.^2 .* T .* exp(-TD ./ (alpha*T)) / (M * v^2 * TD)).^(-1);
% ***
% % v.2. according to Klemens
% tU = (a0 * 8 * g^2 * kB * TD * v / M * T .* exp(- TD ./ ( alpha * T ))).^(-1);
disp('t Umklapp')
display(max(max(tU)))

% scattering due to point defect
A = c * (dM / M)^2 * 4 * pi^3 * a^3 / v^3;
tPD = (A * w'.^4).^(-1);     % t^(-1) = v * P (probability by Klemens)
disp('t point defect')
display(max(max(tPD)))

% boundary scattering
tB = (v / d)^(-1);
disp('t boundary')
display(max(max(tB)))

% summation of phonon relaxation time
tq = (tU.^(-1) + tPD.^(-1) + tB.^(-1)).^(-1);
%tq = (tU.^(-1) + tB.^(-1)).^(-1);
disp('t sum')
display(max(max(tq)))

% phonon-phonon normal scattering the relaxation rate
a1 = 1e-11;  % constant
tN = 1 ./ (a1 * w' .* T.^3);   % for LiF and diamond
disp('t phonon-phonon')
display(max(max(tN)))

% combined relaxation time
tc = (tq.^(-1) + tN.^(-1)).^(-1);    
disp('t combined')
display(max(max(tc)))

% dx = zeros(length(w), length(T));
% I1 = zeros(length(w), length(T));
% I2 = zeros(length(w), length(T));
% I3 = zeros(length(w), length(T));
%for i = 1 : length(T)
%     dx(:, i) = hr * dw ./ (kB * T(i));
% 
%     I1(:, i) = tc .* x.^4 .* exp(x) ./ (exp(x) - 1).^2 * dx(:, i);
%     I2(:, i) = tc ./ tN .* x.^4 .* exp(x) ./ (exp(x) - 1).^2 * dx(:, i);
%     I3(:, i) = tc ./ (tN .* tq) .* x.^4 .* exp(x) ./ (exp(x) - 1).^2 * dx(:, i);
%end
    
dx = hr * dw ./ (kB * T);

I1 = trapz(tc .* x.^4 .* exp(x) ./ (exp(x) - 1).^2, 1) .* dx;
I2 = trapz(tc ./ tN .* x.^4 .* exp(x) ./ (exp(x) - 1).^2, 1) .* dx;
I3 = trapz(tc ./ (tN .* tq) .* x.^4 .* exp(x) ./ (exp(x) - 1).^2, 1) .* dx;

    
k1 = c1 * c2 * T.^3   .*   I1;
k2 = c1 * c2 * T.^3   .*   I2.^2 ./ I3;
k = k1 + k2;

% I = tq .* x.^4 .* exp(x) ./ (exp(x) - 1).^2 * dx;
% kL = c1 * c2 * T.^3 .* I;

% display the temperature @ peak
%display(find(kL == max(kL)))
% display (max(k_num_8A))
% display (max(kL))

figure
lw = 1.3;
plot(  T, k_num_8A, 'r--', T, k, 'b',  'LineWidth', lw)
hold on
plot(T_Ber, k_Ber,'-s','MarkerSize',5, 'color',[0.0 0.5 0.5],  'LineWidth', lw)
plot(T_Pohl, k_Pohl,'s','MarkerSize',5, 'color',[0.2 0.0 0.7],  'LineWidth', lw)

title('LiF T dependent thermal conductivity')
xlabel('temperature, K')
ylabel('thermal conductivity, W/mK')
legend('k_L_i_F DFT', 'k_L_i_F Debye-Klemens', 'experimental Berman 1956', ...
    'experimental Pohl 1960') 
hold off
display(find(k == max(k)))
display(max(k))
% display(k(120))
%plot(T, k, 'k',  'LineWidth', lw)
%title('LiF T dependent thermal conductivity Debye-Klemens')
%xlabel('temperature, K')
%ylabel('thermal conductivity, W/mK')
% axis([0 200 0 3000])

%__________________________________________________________________________
% Fit accuracy
% mean_dat = mean(k);
% SS_res = sum((k' - k_num_8A).^2);
% SS_tot = sum((k' - mean_dat ).^2);
% 
%ss = sqrt(sum((k' - k_Ber).^2));
% aberr = ss / mean(k);
% display(acc)
%__________________________________________________________________________

