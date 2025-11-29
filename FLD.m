%% this code is dedicated to compute estimates of forces and moments

Num_fin = 4; % Number of Fins
AR = ; % Aspect Ratio
Mach = ;% Mach Number
Lc = ; % Midchord Line Sweep Angle
Lr = ; % reference length
Lf = ; % fin mid-chord line
Af = ; % Fin Area
Ar = ; % Reference Area
T1 = ; % distance from nose tip to fin root chord leading edge
xt = ; % distance between fin root leading edge and fin tip leading edge parallel to body
ct = ; % fin tip chord
cr = ; % fin root chord
Lr = ; % reference length
lambda = ct/cr; % fin tamper ratio ct/cr
rt = ; % body radius at tail
Yt = ; % Spanwise Center of Pressure
s = ; % fin semispan
R = ; % radius of body
d = 2*R; % diameter
A_BN = ; % nose bsae area
AB = ; % base area
VB = ; % volume of body
L0 = ; % Total Body Length
LN = ; % length of nose
% beta = sqrt (M^2-1), supersonic; sqrt (1-M^2), subsonic;
if Mach > 1
    beta = sqrt (Mach^2-1);
else
    beta = sqrt (1-Mach^2);
end

%% Tail %%

% comments:
% cp of tail does not depend on number of fins
 
% Tail Normal Force Coefficient Derivative
CN_alpha1 = 2*pi/(2+sqrt(4+(beta*AR/cos(Lc))^2)); % subsonic for one-fin
CN_alpha2 = (Num_fin*pi*AR(Af/Ar))/(2+sqrt(4+(beta*AR/cos(Lc))^2)); % subsonic for multi-fin
CN_alpha3 = (Num_fin/2)*(CN_alpha1);  % Meant for 3-4 Fin Configuration
% Supersonic Condititons using Strip Theory (requires sims - OpenRocket,...
% AeroRocket, etc.)

% Tail Normal Force Coefficient Derivative (another form)
CN_F = 4*N*(1+R/(s+R))*((s/d)^2/(1+sqrt(1+(2*Lf/(cr+ct))^2)));

% Tail Center of Pressure
X_T = T1+(xt/3)*((cr+2*ct)/(cr+ct))+(1/6)*(cr+ct-cr*ct/(cr+ct));

% Tail Roll Forcing Moment Coefficient
C_Taill_subfoRoll = Num_fin*(CN_alpha1)*(Yt/Lr);          % Subsonic Forcing Moment coefficient

% Tail Roll Damping Moment Coefficient
C_Tail_subRoll = (-1)*(Num_fin*cr*s/(6*Lr^2))*(CN_alpha1)*((1+3*lambda)*s^2+4*(1+2*lambda)*s*rt+6*(1+lambda)*rt^2);% tail roll damping moment coefficient derivative - subsonic
C_Tail_supRoll = 1000*Num_fin*C_l; % Tail roll damping moment coefficient derivative- supersonic

%% Nosecone %%

% comments:
% the nose is overlooked often
% this part may be neg

% Nose Normal Force Coefficient Derivative (simplifed constant)
CN_nose = 2 ;

% Nose Center of Pressure
Nosetype = "cone";

if Nosetype == "cone"
    XN = (2/3)*LN;
elseif Nosetype == "ogive"
    XN = (0.466)*LN;
elseif Nosetype == "parabolic"
    XN = (0.5)*LN;
end
%% Body %%

% Body Normal Force Coefficient Derivative in Subsonic Flow
CN_alpha_body = 2* A_BN/Ar; % normal force coefficient derivative

% Body Center of Pressure
%X_B = Lr*CM_body/CN_alpha_body;
X_B = L0-VB/AB;

% Body Pitching Moment Coefficient
CM_body = X_B*CN_alpha_body/Lr;