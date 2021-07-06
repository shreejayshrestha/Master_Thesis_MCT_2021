
%% Program to calculate the natural frequency of a thin plate and floating floor
% % Boundary condition: simply supported
% Created Feb-May 2021
% Author: shreejay
% shreejayshrestha@gmail.com

clc
clear variables
%% Thin plate
% E = 2.96e9; % youngs modulus unit Pa
%  E = 1.2e9; % 
%   E = 1.71e9; % youngs modulus unit Pa
% E = 1.68e9; % youngs modulus unit Pa PS: this value seems to be close enough with the f0 observed
%  E = 2.77e9; % very close : finally commited on 03:06:2021
E = 25e9; % for f_ef
% E = 2.78e9; % very close
v = 0.18; % Poisson's ratio
h = 0.14; % thickness of the plate/concrete slab
B = E*h^3/(12*(1-v^2)); % bending stiffness

%% Thin plate
a = 3.61;% length of the plate
b = 3.12;% width of the plate
% rho_conc = 2300;  % kg/m^3 density of concrete slab
rho_conc = 2450;  % kg/m^3 density of concrete slab

m = rho_conc * h; % mass per unit area of the concrete slab

nxmax = 10;
nzmax = 10;

nx = 1: nxmax;
nz = 1:nzmax;

floor_natfreq_analytical = zeros(9,1);
c = 1;
for i = 1:3
    for j = 1:3
        floor_natfreq_analytical(c,1) = (pi/2)*sqrt(B/m)*((nx(i)/a)^2 + (nz(j)/b)^2);
        c=c+1;
    end
end

floor_natfreq_analytical= sort(floor_natfreq_analytical);
floor_f0_analytical= floor_natfreq_analytical;
% plate_nat_freq_ansys = [52.09
% 118.45
% 141.4
% 203.67
% 227.82
% 286.75
% 308.21
% 344.88
% 376.57
% ];
%% Export
% save Final_Output_from_allscripts_23_04_2021\floor_analytical_natural_freq\floor_natfreq_analytical.mat floor_natfreq_analytical% cell1 is mf and cel2 is ff.
%  save Final_Output_from_allscripts_08_05_2021\acc\floor_f0_analytical.mat floor_f0_analytical

%% Floating floor
lplate = 1.25; % meter
bplate = 0.72; % meter
tplate = 0.05; % meter
mplate = rho_conc * tplate; % mass per unit area of the concrete slab in the floating floor
trockwool = 20; % thickness in mm
dsrockwool = 12.8; % dyn.stiff. MN/m^3 for 20x600x1200 mm trynnlydplate source refer anders buen email
dseps25 = 20; % MN/m^3 source Schutz Quadro Takk 25mm, canes webpage
dseps50 = 10; % MN/m^3 source Schutz Vari Takk Pro 50mm, canes webpage

dseff = (dsrockwool*dseps25*dseps50)/(2*dseps25*dseps50+ dsrockwool*dseps25+dsrockwool*dseps50);
dsair = 115/(trockwool) *1000000;
% meff = mplate + m;
meff = m;

% % F0_ff = 1/(2*pi) * sqrt((dseff+111/(2*trockwool)) *1000000 * (1/(mplate) + 1/(m)));
F0_ff = 1/(2*pi) * sqrt((dseff+111/(2*trockwool+0.075)) *1000000 * (1/(mplate) + 1/(m)));
F0_ff_1 = 1/(2*pi) * sqrt((dseff+115/(2*trockwool)) *1000000 * (1/(mplate) + 1/(m)));

F0_ff2 = 1/(2*pi) * sqrt((dseff+115/(2*trockwool) *1000000) * (1/(mplate+m)));



ke_calc = 1/12.8 + 1/12.8 + 1/10 + 1/20;
ke = 1/ke_calc ;
F0_ff3 = 1/(2*pi)* sqrt((ke+dsair)/meff);
Ln_improve = 1.5*F0_ff3;
Ln_imp_foff = 1.5*F0_ff;
%% Natural freq. of RC slab on elastic fondation
ssr = 6000 ; % kN/m^3 soil subrage reaction
nx = 1: nxmax;
nz = 1:nzmax;

floor_natfreq_analytical = zeros(9,1);
f_ef = zeros(9,1);
c = 1;
for i = 1:3
    for j = 1:3
        f_ef(c,1) =  (1/2*pi) * sqrt(pi^4*B*(1/m)*(((nx(i)/a)^2 + (nz(j)/b)^2)^2)+ssr/m );
        c=c+1;
    end
end
