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

g = 1.6;                        % Gruneisen constant of LiF, a measure of anharmonicity
a = 4.03e-10;                   % a^3 atomic volume in m
d = 1;                          % dimensions of the sample in m 

%T = [1 : 100 110 : 10 : 1000];              % temperature
T = 300;
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

NA = 6.02214076*10^23;              % Avogadro constant

% 
% wm=wD;                                % maximum frequency set equal to Debye
% ww=(1:50)/50;                         % setting frequency range
% dwm=wm/50;w=wm*ww;      

% ============================ load data ==================================

% % load DFT simulation data
% k_num = load('numdata.txt');
% k_num_4A = k_num(:,1);
% k_num_7A = k_num(:,2);
% k_num_8A = k_num(:,3);
% 
% 
% % load Berman experimental simulation data
% k_Ber = load('LiF_k_Berman.txt');
% T_Ber = k_Ber(:,1);
% k_Ber = k_Ber(:,2)* 1e+2;
% 
% % load Pohl experimental simulation data
% k_Pohl = load('LiF_k_Pohl.txt');
% T_Pohl = k_Pohl(:,1);
% k_Pohl = k_Pohl(:,2)* 1e+2;

% load defects concentration data
ndef = load('lif_defects.txt');
ff = ndef(:,1);
ll = [250 317	377	445	450	518	 530 540];  % spectral positions of defects
nF = ndef(:,2:9)*1e+18;
nF(4,1) = nF(4,1)-0.02e+19;
nF(:,2) = zeros(length(ff), 1);

fff = logspace(0,14,201)';
% parameters of the fit for 250 F
% A = [2.89e+19 7.0e+18 5.2e+18 1.09e+19 2.9e+19 2.8e+18 4e+17 2e+18];
% A1 = [1.0e+11 1.5e+10 2.0e+11 1.65e+11 1.6e+11 2e+11 0.1e+11 0.25e+11];
%
% A = [7.85e+19 7.85e+18 1.51e+19 3.025e+19 8.025e+19 7.75e+18 1.11e+18 5.5e+18];
% A1 = [3.5e+11 1.25e+10 4.0e+11 3.6e+11 3.6e+11 3.6e+11 3.6e+11 3.6e+11];
% A2 = [4.6e+11 4.6e+11 4.6e+11 4.6e+11 4.6e+11 4.6e+11 4.6e+11 4.6e+11];
% %               I          II          III         IV            V        VI
% A =         [39.6e+18   0.75e+19    0.675e+19    1.8e+19     3.65e+19   5.5e+18    1.11e+18    5.5e+18];
% A1 =        [1.e+11    2.0e+10     1.8e+11      2.5e+11     2e+11      3.6e+11     3.6e+11     3.6e+11];
% A2 = A1 .*  [10         100           15.0         4.6         15         4.6         4.6         4.6];
%               I          II          III         IV            V        VI
A =         [39.6e+18   0.75e+19    0.675e+19    1.8e+19     3.65e+19   5.5e+18    1.11e+18    5.5e+18];
A1 =        [1.e+11    2.0e+10     1.8e+11      2.5e+11     2e+11      3.6e+11     3.6e+11     3.6e+11];
A2 = A1 .*  [10         100           15.0         4.6         15         4.6         4.6         4.6];


nnF = A .* (1 - exp(-fff ./ A1))./ (1 + fff ./ A2);
%nnF = A .* (1 - exp(-fff ./ A1));

% special fit for F3
A(2) = 0.0e+19;
A1(2) = 2.0e+10;
A2(2) = 5e+10;%0.15e+11;
%nnF(:,2) = A(2) .* (1 - exp(-fff ./ A1(2)));
nnF(:,2) = A(2) .* exp(-(fff-A2(2)).^2 ./ A1(2)^2);
%nnF(:,2) = flip(nnF(:,2));
%nnF(:,2) = A(2) .* (1 - exp(-fff ./ A1(2)))./ (1 + fff ./ A2(2));


% figure
% axes_handle_1 = axes;
% axes_position = get(axes_handle_1, 'Position');
% axes_handle_2 = axes('Position', axes_position);
% % semilogx(axes_handle_1, fff(1:ind1), nnF(1:ind1, :), 'Linewidth', 1.3)
% % semilogx(axes_handle_2, fff(ind2:end),nnF(ind2:end, :), 'LineWidth', 1.3);
% q = semilogx(ff, nF(:,:), 'o', 'Linewidth', 1.3);
% q(1).Color = [0 0.4470 0.7410];
% q(2).Color = [0.8500 0.3250 0.0980];
% q(3).Color = [0.9290 0.6940 0.1250];
% q(4).Color = [0.4940 0.1840 0.5560];
% q(5).Color = [0.4660 0.6740 0.1880];
% q(6).Color = [0.3010 0.7450 0.9330];
% q(7).Color = [0.6350 0.0780 0.1840];
% q(8).Color = [0 1 0];
% 
% 
% 
% hold on
% q2 = semilogx(fff, nnF(:,:), '--', 'Linewidth', 1.3);
% q2(1).Color = [0 0.4470 0.7410];
% q2(2).Color = [0.8500 0.3250 0.0980];
% q2(3).Color = [0.9290 0.6940 0.1250];
% q2(4).Color = [0.4940 0.1840 0.5560];
% q2(5).Color = [0.4660 0.6740 0.1880];
% q2(6).Color = [0.3010 0.7450 0.9330];
% q2(7).Color = [0.6350 0.0780 0.1840];
% q2(8).Color = [0 1 0];
% %semilogx(fff, nnF(), '--', 'color', [0.5 0.5 0.5], 'Linewidth', 1.3)
% hold off
% 
% legend([q(1) q(2) q(3) q(4) q(5) q(6) q(7) q(8)], 'F (250 nm)', 'F_3_(_I_) (317 nm)',	'F_3_(_I_I_) (377 nm)',	...
%     'F_2 (445 nm)',	'Colloids (450 nm)',	'F_4_(_I_) (518 nm)',	...
%     'Colloids (530nm)', 'F_4_(_I_I_) (540 nm)')
% legend('Location','northwest')
% 
% %Link the y axis limits and fontsize property of the axes objects
% linkaxes([axes_handle_1 axes_handle_2], 'y');
% linkprop([axes_handle_1 axes_handle_2], 'FontSize');
% %Set the x range limits and tick mark positions of the first axes object
% set(axes_handle_1, 'XLim', [1 100], ...
%       'XTick', (1))
% %Set the x range limits and tick mark positions for the second axes object.
% %Also set the background color to 'none', which makes the background
% %transparent.
% set(axes_handle_2, 'Color', 'none', ...
%       'YTickLabel', [], ...
%       'XLim', [1e+9 2.5e+13], ...
%       'XTick', [ 1e+10 1e+11 1e+12 1e+13])
% %xlabel('fluence, cm^-^2')
% ylabel('Defects concentration, cm^-^3')
% set(gca,'FontSize',14) 
%axis([1 2.5e+13 0 3e+19])
%xticks([1 10^10 10^13])
% x0=10;
% y0=10;
% width=550;
% height=400;
% %set(gcf,'position',[x0,y0,width,height])
% set(gcf,'position',[width,height])

% TDTR measured values
temp = load('LiF_exp2.txt');
fl_exp = temp(:,1);
fl_exp(1) = 1;
k_exp = temp(:,2);
dk_exp = temp(:,3);

% MTR measured values
k_exp_mtr = [11 10.13 7.7 8.69 6.84];
fl_exp_mtr = fl_exp(1:end-1);
fl_exp_mtr(2) = 2e+10;

% =========================== calculation =================================
%T = T_Pohl';
% FITTING PARAMETERS

% estimated
dM = [MF 3*MF 3*MF 2*MF MLi 4*MF MLi 4*MF];      % difference in mass due to defects
%dM = MF;
% % c1 = 2.5e-2;
% % beta = 0.006;
% 
% % %c = linspace(1.55e-5, 7e-4, 201);   % the point defect concentration (per atom)
% % c1 = 4.5e-3;
% % beta = 0.06;
% % %beta = 1;
% % D0 = 2e-2;                         % Marat's model parameter
% % % stretched the tail up (smaller) or down (larger)
% 
% %c1 = 7.5e-2;
% c1 = 7.5e-2;
% 
% %beta = 0.015;
% beta = 0.025;
% 
% %D0 = 1.8e-2;                         % Marat's model parameter
% D0 = 0.25e-28;                         % Marat's model parameter


% parameters calculated by optimization program
% D0 = 9.6667e-29;
% beta = 0.005;
% c1 = 0.033;

% %D0 = 3.1e-29;
% D0 = 3.5354e-29;
% beta = 2.75000e-02;
% c1 = 0.0172;

D0 = 3.1000e-29;
beta = 2.75000e-02;
c1 =0.0272;
% ------------------phonon scattering relaxation time ----------------- 

% Umklapp scattering
% v.1. according to Tritt
% coefficients
%a0 = 1.1e+28;       % fixed by T dependent model                      
%a0 = 1.175e+28;       % fixed by T dependent model                      
a0 = 1.033e+28;       % fixed by T dependent model                      
alpha = 3.2;        % fixed by T dependent model

tU = (a0 * hr * g^2 * w'.^2 .* T .* exp(-TD ./ (alpha*T)) / (M * v^2 * TD)).^(-1);
% ***
% % v.2. according to Klemens
% tU = (a0 * 8 * g^2 * kB * TD * v / M * T .* exp(- TD ./ ( alpha * T ))).^(-1);
% disp('t Umklapp')
% display(max(max(tU)))

% scattering due to point defect
c = c1*nnF / NA;
%A0 = 3.1261e-44;   % in T-dependent case
A0 = sum( c .* (dM / M).^2 * 4 * pi^3 * a^3 / v^3, 2);
tPD = (A0 .* w.^4).^(-1);     % t^(-1) = v * P (probability by Klemens)
%tPD = 1;
% disp('t point defect')
% disp(max(max(tPD)))

% boundary scattering
tB = (v / d)^(-1);
% disp('t boundary')
% display(max(max(tB)))

% track-core scattering 
rF = 22e-7;     % radius of F center cylinder in cm
sigma= pi * rF^2 * beta;   % damaged cross-section in cm-2

%D=D0./(exp(- sigma*fff));            % v.1 according to Marat's model
D=D0.*(1 - exp(- sigma*fff));            % v.1 according to Marat's model
tT = (D.*w.^3).^(-1);       
%tT = 1;
% disp('t tracks')
% disp(max(max(tT)))

% Nd = c(:,1);                            % v.2 according to the Chernatynskiy et al. 
% tT = (Nd .* rF^4/v^2 * w.^3).^(-1);

% summation of photon relaxation times withour cylindrical tracks, point
% defects only
tqPD = (tU'.^(-1) + tPD.^(-1) + tB.^(-1)).^(-1);

% summation of photon relaxation times
tq = (tU'.^(-1) + tPD.^(-1) + tB.^(-1) + tT.^(-1)).^(-1);
% disp('t sum')
% display(max(max(tq)))

% phonon-phonon normal scattering the relaxation rate
BN = 1e-11;  % constant
tN = 1 ./ (BN * w .* T.^3);   % for LiF and diamond
% disp('t phonon-phonon')
% display(max(max(tN)))

% combined relaxation time
tc = (tq.^(-1) + tN.^(-1)).^(-1);   

% combined relaxation time with point defects only
tcPD = (tqPD.^(-1) + tN.^(-1)).^(-1); 

% disp('t combined')
% display(max(max(tc)))

   
x = hr * w ./ (kB * T);            % dimensionless parameter
dx = hr * dw ./ (kB * T);           % discritization 

% point defects only
I1PD = trapz(tcPD .* x.^4 .* exp(x) ./ (exp(x) - 1).^2, 2) .* dx;
% sum of point defects and tracks
I1 = trapz(tc .* x.^4 .* exp(x) ./ (exp(x) - 1).^2, 2) .* dx;

% point defects only
kPD = kB / (2 * pi^2 * v) * (kB / hr)^3 * T.^3   .*   I1PD ;

% all terms included
k = kB / (2 * pi^2 * v) * (kB / hr)^3 * T.^3   .*   I1 ;

% thermal conductivity of the pure crystal
kpure = k(1) * ones(length(k),1);

% % load thermal conductivity data with pure point defect effect
% kPD = load('kPD.txt');

figure
%Create two overlapping axes
axes_handle_1 = axes;
axes_position = get(axes_handle_1, 'Position');
axes_handle_2 = axes('Position', axes_position);
% 
%Plot the two sections of data on different axes objects
lw = 1.3;
fsz = 20;
% plot(fff, k, 'r', 'LineWidth', lw)
% xlabel('number of defects per atom')
% ylabel('k (Wm^-^1K^-^1)')
ind1 = 4;
ind2 = 141;
semilogx(axes_handle_1,fff(1:ind1),k(1:ind1), 'k--',  ...
    fl_exp(1),k_exp(1),'ro', fl_exp_mtr(1),k_exp_mtr(1),'bs', 'LineWidth', 1.2);
%semilogx(axes_handle_1,fff(1:ind1),k(1:ind1), 'k--', fl_exp_mtr(1),k_exp_mtr(1),'bs', 'LineWidth', 1.2)

% hold on
% errorbar(fl_exp(1), k_exp(1), dk_exp(1),'ro', 'LineWidth', lw)
% hold off

p = semilogx(axes_handle_2, fff(ind2:end),k(ind2:end), 'k--',...
    fff(ind2:end),kPD(ind2:end), '--', fff(ind2:end),kpure(ind2:end), ':',...
     fl_exp_mtr(2:end), k_exp_mtr(2:end), 'bs', 'LineWidth', 1.3);
p(2).Color = [0.5 0.5 0.5];
p(3).Color = [0.5 0.5 0.5];
%fl_exp(2:end),k_exp(2:end),'ro',
hold on
%errorbar(fl_exp, k_exp, dk_exp,'s', 'Parent',ax2, 'LineWidth', lw)
p2 = errorbar(fl_exp(2:end), k_exp(2:end), dk_exp(2:end),'ro', 'LineWidth', lw);
hold off
legend([p(3) p(2) p(1) p(4) p2], 'pristine', 'point defects', 'point defects + cylindrical tracks',...
    'MTR', 'TDTR')
legend('Location','southwest')
%Link the y axis limits and fontsize property of the axes objects
linkaxes([axes_handle_1 axes_handle_2], 'y');
linkprop([axes_handle_1 axes_handle_2], 'FontSize');
%Set the x range limits and tick mark positions of the first axes object
set(axes_handle_1, 'XLim', [1 100], ...
      'XTick', (1))
%Set the x range limits and tick mark positions for the second axes object.
%Also set the background color to 'none', which makes the background
%transparent.
set(axes_handle_2, 'Color', 'none', 'YLim', [4 13],...
      'YTickLabel', [], ...
      'XLim', [1e+9 2.5e+13], ...
      'XTick', [ 1e+10 1e+11 1e+12 1e+13])
set(gca,'FontSize',fsz)  
xlabel('fluence, cm^-^2')
ylh = ylabel('k, Wm^-^1K^-^1');
ylh.Position(1) = ylh.Position(1) - 2;
% figure
% lw = 1.3;
% fsz = 16;
% ind2 = 130;
% semilogx(fff(ind2:end),k(ind2:end), 'k--',fl_exp(2:end),k_exp(2:end),...
%     'ro', fl_exp(2:end-1), k_exp_mtr, 'bs', 'LineWidth', 1.3);
% set(gca,'FontSize',fsz)  