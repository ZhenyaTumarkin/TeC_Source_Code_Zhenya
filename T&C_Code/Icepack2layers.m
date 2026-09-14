%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%   Subfunction  Compute Ice  Balance       %
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function[TsF,ICE,ICE_D,IP_wc,WR_IP,dQ12,Qfm,Im,Tdpice,Tice]=Icepack2layers(dt,Ts,Tstm1,ICEtm1,IP_wctm1,...
    Pr_liq,EICE,Rn,H,QE,G,Qv,Cwat,Ccrown,Cfol_H,Csno,Cicew,Ice_wc_sp,WR_SP,Tdpicetm1,Ticetm1,hSTL,min_SPD)
%%%INPUTS
%dt time step [s]
dth = dt/3600; %% [h]
%%% OUTPUTS
%%% PARAMETERS
ci = 2093; %%% [J/Kg K] specific heat ice
lan_ice=2.2; % W m^-1 K^-1
row = 1000; % water density [kg/m^3]
roi = 916.2; %% ice density [kg/m^3]
Lf= 1000*333.7; %% [J/Kg] %% Latent heat melt/freezing

ICE_D_tm1=0.001*ICEtm1*row/roi; % [m]   %% [m] last ice tghicnkess
%L= 1000*(2501.3 - 2.361*(Ta)); %%% Latent heat vaporization/condensaition [J/kg]
%Ls= Lf +L; %%% Latent Heat of sublimation  [J/kg]
%%% WR_SP Water Released from Snowpack [mm]
%%%%%%%%%%%%%
%%%%%%%%% Liquid Precipitation on Ice
Pr_liq= Pr_liq*dth; %% [mm]
%%%%%%%%%%%%%%%%%%%%%%%%%%
%%% UPDATE ICE amount
ICE = ICEtm1 - EICE*dth;  %%% [mm] Icepack After Sublimation and new
%%% for this script Csno = 0
Csno=0;

%%%%%%%%%%% ICEPACK ENERGETIC AND MASS BALANCE
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

G=0; %ASSUME at 2m the temp profile is stabilised (ie dT/dz ~= 0)
Gbot=G;

%%build discretisation
ms_sno=5; % length(Tdpsnowtm1);
ICE_D_eb=max(2,ICE_D_tm1);
Zs_sno= [0 exp([linspace(log(1000*hSTL),log(1000*ICE_D_eb),ms_sno)])];  %% [mm] ***


lan_ice =lan_ice*ones(1,ms_sno); %%% [W/m K ] Thermal conductivity debris
% specific heat capacity ice = 2093 J / kg K
cv_ice = roi*ci*ones(1,ms_sno); % [J m^-3 K^-1] Volumetric heat capcity snow
[G0,Tdpice]=Soil_Heat_Profile_Normal(Ts,dt,Tdpicetm1,ms_sno,Zs_sno,lan_ice,cv_ice,NaN,NaN,Gbot,3);
Tdpice(Tdpice>0)=0; %%% energy used for melting likely 


%Tsnow=sum(Tdpsnow.*diff(Zs_sno)/(1000*Dtm1));

%%% Surface budget computed on top min_SPD/2 (not as relevant as the temperature discretization)
D_eb= min_SPD/2; % %% [m]
ICEteb = D_eb*roi/row*1000; %% [mm]
ICEteb=min(ICEteb,ICE); 




    %%%% ENERGY BALANCE
    %Qfm = (Ts < 0)*(row*(IP_wc_L/1000)*Lf/dt) - (Ts >= 0)*(row*(IP_wc_F/1000)*Lf/dt); %%% [W/m^2]  %% Freezing/Melting Water in the icepack (not working)
% Qfm=0; %%%
% dQ = Rn + Qv + Qfm - H - QE -G;   % [W/m^2]  Variation Heat in the Icepack
% dTs = (1000*dQ*dt)/(ci*(ICEeb)*row); %% Variation temperature in the Icepack [�C]
% TsF = Tstm1 + dTs ; %%  Temperature Icepack [�C]
% Qcc =ci*(ICEeb)*TsF*row/1000; %%  Icepack Cold Content [J �C/m^2 K]
% %%%%% MASS BALANCE -- TEMPERATURE UPDATE
% if Qcc <= 0;
%     Im =0 ; % Ice melt [mm]
% else
%     TsF = 0;
%     Im = Qcc/Lf; %  Ice melt [Kg/m^2]
%     Im=  1000*(Im/row); %% Ice melt [mm]
%     Im = min(Im,ICE);
%     ICE =ICE - Im;
% end


%%%%%%%%%
dQ1 = Rn + Qv - H - QE -G0;   % [W/m^2]  Variation Heat in the Snowpack
dTs = (1000*dQ1*dt)/(ci*(ICEteb)*row); %% Variation temperature in the snowpack [�C]
TsF = Tstm1 + dTs ; %%  Temperature Snowpack [�C]
Qcc =ci*(ICEteb)*TsF*row/1000; %%  Snowpack Cold Content [J �C/m^2 K]
%%%%% MASS BALANCE -- TEMPERATURE UPDATE
if Qcc <= 0
    Im_top =0 ; % Snow melt [mm]
    dQres =0 ;
else
    TsF = 0;
    Im = Qcc/Lf; %  Snow melt [Kg/m^2]
    Im=  1000*(Im/row); %% Snow melt [mm]
    Im_top = min(Im,ICEteb);
    %%%%%%%%%%%%%%%%%%%%%%%%%
    dQres= Qcc/dt -  Lf*(Im_top)*row/(1000*dt) ; %% [W/m^2] Additional Energy to melt below top layer
    %Qcc =0; % Icepack Cold Content [J �C/m^2 K]
    dQ1=dQ1-dQres; 
end


%%%% Icepack water content 
f_frozen_to_liquid = (Tdpice >= -0.01).*(Tdpicetm1< -0.01).*diff(Zs_sno)/(1000*2); %%%% ASSUME ICE_D>2m always
f_liquid_to_frozen = (Tdpicetm1 >= -0.01).*(Tdpice < -0.01).*diff(Zs_sno)/(1000*2); 
IP_wc_L = sum(f_liquid_to_frozen)*IP_wctm1; %%% Snowpack Liquid to frozen in the snowpack [mm]
IP_wc_F = sum(f_frozen_to_liquid)*IP_wctm1; %%% Snowpack Frozen to liquid in the snowpack [mm]
%%%%%%%%%%



Qfm = (row*(IP_wc_L/1000)*Lf/dt) - (row*(IP_wc_F/1000)*Lf/dt); %%% [W/m^2]  %% Freezing/Melting Water in the snowpack
ICEeb= max(0.1,ICE-ICEteb); %%%[mm] Minumum snow-pack layer to compute energy budget temperature variations  
ICEeb= min(2000,ICEeb); %%%[mm]  Maximum Snow for computing energy budget (no more than 2m)

%%%% ENERGY BALANCE
dQ2 = G0 +  Qfm +  dQres -G;   % [W/m^2]  Variation Heat in the Snowpack
dTs = (1000*dQ2*dt)/(ci*(ICEeb)*row); %% Variation temperature in the snowpack [�C]
Tice = Ticetm1 + dTs ; %%  Temperature Snowpack [�C]
Qcc =ci*(ICEeb)*Tice*row/1000; %%  Snowpack Cold Content [J �C/m^2 K]
%%%%% MASS BALANCE -- TEMPERATURE UPDATE
if Qcc <= 0
    dQres =0 ; 
    Im=0;
else
    Tice = 0;
    Im = Qcc/Lf; %  Snow melt [Kg/m^2]
    Im=  1000*(Im/row); %% Snow melt [mm]

   
    %%%%%%%%%%%%%%%%%%%%%%%%%
    dQres= Qcc/dt -  Lf*(Im)*row/(1000*dt) ; %% [W/m^2] Energy Lost -- theroetically to increase TsF or G
    %Qcc =0; % Snowpack Cold Content [J �C/m^2 K]
end



dQ12 = dQ1+dQ2; 
ICE =ICE -Im - Im_top;
dTice=Tice-sum(Tdpice.*diff(Zs_sno)/(1000*2));  %% assumes ICE_D>=2m

Im=Im+Im_top;   %%% add melt together


%%%%%%%%%%
%%% ICE-DEPTH AND DENSITY
ICE_D= 0.001*ICE*(row/roi); %% Icepack [m]
%%%%%%%% LIQUID WATER SNOWPACK
IP_wc_max = Ice_wc_sp*ICE ; %% [mm] Ice water retention factor
KT= ICE/1000; %% [h] %% Vertical Glacier Flow Time Constant  V = 1.0 [m/h] 
IP_outflow = min(IP_wctm1,(IP_wctm1/KT));  %%% [mm/h] 
IP_wc =  IP_wctm1 + Im + WR_SP - IP_outflow*dth + Pr_liq*(1-Csno)*(1 - (1-Cicew)*Cwat - sum(Ccrown.*Cfol_H));%% Ice Pack Water Content [mm]
WR_IP = (IP_wc - IP_wc_max)*(IP_wc >= IP_wc_max); %% Water Released from Icepack filling [mm]
IP_wc =  IP_wc - WR_IP; %%[mm] %% Snow Pack Water Content [mm]
WR_IP = IP_outflow*dth + WR_IP; %% Total Water Released from Icepack [mm]
%%%%%%%%%%%%%%%%%%%%%%%%%
end




