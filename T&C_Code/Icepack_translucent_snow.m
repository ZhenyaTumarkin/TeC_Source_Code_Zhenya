%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%   Subfunction  Compute Ice Balance        %
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Energy balance on ice with snow surface properties,
% mass balance on snow. All melt is assigned to snow.
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function [TsF,ICE,ICE_D,IP_wc,WR_IP,dQ12,Qfm,Sm,Tdpice,Tice,...
          SWE,In_SWE,D,ros,SP_wc,WR_SP,U_SWE,NIn_SWE] = ...
Icepack_translucent_snow(dt,Ts,Tstm1,ICEtm1,IP_wctm1,...
    Pr_liq,Rn,H,QE,G,Qv,Cwat,Ccrown,Cfol_H,Csno,Cicew,...
    Ice_wc_sp,Tdpicetm1,Ticetm1,Tdew,Ws,t_sls,ESN_In,...
    In_max_SWE,In_SWEtm1,SWEtm1,ESN,Dtm1,rostm1,ros_max1,...
    ros_max2,hSTL,min_SPD,Ta,...
    dw_SNO,SP_wctm1,Pr_sno)

%% Parameters
dth = dt/3600;
ci = 2093; 
lan_ice=2.2 ;% W m^-1 K^-1
row = 1000; roi = 916.2; Lf = 333.7e3;
tiny = 1e-10;

%% Initialise states and diagnostics
ICE = max(0,ICEtm1); SWE = max(0,SWEtm1);
IP_wc = max(0,IP_wctm1); SP_wc = max(0,SP_wctm1);
Sm = 0; Sm_top = 0; Sm1 = 0; Sm2 = 0;
WR_SP = 0; WR_IP = 0; Qfm = 0; dQres = 0;

Pr_sno = max(0,Pr_sno*dth);
Pr_liq = max(0,Pr_liq*dth);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% Snow interception and SWE update
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

[U_SWE,In_SWE,NIn_SWE] = Snow_Interception(dth,Pr_sno,Tdew,Ws,...
    t_sls,ESN_In,In_max_SWE,In_SWEtm1);

SWE = SWE + Pr_sno*(1-(1-Csno)*Cwat) - NIn_SWE + U_SWE - ESN*dth;
SWE = max(0,SWE);

% Retained from original logic
Csno = 0;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% Ice temperature profile
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

G = 0; Gbot = G;
ICE_D_tm1 = 0.001*ICEtm1*row/roi;

ms_sno = 5;
ICE_D_eb = max(2,max(hSTL,ICE_D_tm1));
Zs_sno = [0 exp(linspace(log(1000*hSTL),log(1000*ICE_D_eb),ms_sno))];

lan_ice = lan_ice*ones(1,ms_sno);
cv_ice = roi*ci*ones(1,ms_sno); % [J m^-3 K^-1]

[G0,Tdpice] = Soil_Heat_Profile_Normal(min(Ts,0),dt,Tdpicetm1,...
    ms_sno,Zs_sno,lan_ice,cv_ice,NaN,NaN,Gbot,3);

Tdpice = min(Tdpice,0);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% Upper active layer energy balance
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

D_eb = max(min_SPD/2,0.001);
ICEteb = min(D_eb*roi/row*1000,SWE);

dQ1 = Rn + Qv - H - QE - G0;

if ICEteb > tiny
    TsF = Tstm1 + dQ1*dt/(ci*ICEteb);
    Qcc = ci*ICEteb*TsF;

    if Qcc > 0
        TsF = 0;
        Im_top = Qcc/Lf;
        Sm_top = min(Im_top,SWE);
        dQres = max(0,(Qcc - Lf*Sm_top)/dt);
    end
else
    TsF = min(Tstm1,0);
    dQres = max(0,dQ1);
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% Icepack liquid water / phase-change term
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

dz = diff(Zs_sno)/1000;

f_frozen_to_liquid = (Tdpice >= -0.01).*(Tdpicetm1 < -0.01).*dz/2;
f_liquid_to_frozen = (Tdpicetm1 >= -0.01).*(Tdpice < -0.01).*dz/2;

IP_wc_L = sum(f_liquid_to_frozen)*IP_wctm1;
IP_wc_F = sum(f_frozen_to_liquid)*IP_wctm1;

Qfm = Lf*(IP_wc_L - IP_wc_F)/dt;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% Lower icepack energy balance
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

ICEeb = min(2000,max(0.1,ICE-ICEteb));
dQ2 = G0 + Qfm + dQres - G;

Tice = Ticetm1 + dQ2*dt/(ci*ICEeb);
Qcc = ci*ICEeb*Tice;

if Qcc > 0
    Tice = 0;
    Im = Qcc/Lf;

    SWE_available = max(0,SWE-Sm_top);
    In_SWE = max(0,In_SWE);

    if SWE_available + In_SWE > tiny
        r = SWE_available/(SWE_available + In_SWE);
    else
        r = 0;
    end

    Sm1 = min(r*Im,SWE_available);
    Sm2 = min((1-r)*Im,In_SWE);
    In_SWE = In_SWE - Sm2;
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% Snow mass, depth, density and liquid water
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

Sm = Sm_top + Sm1 + Sm2;
SWE = max(0,SWE - Sm_top - Sm1);

dQ12 = dQ1 + dQ2;

[D,ros] = Snow_Depth(dt,Sm_top+Sm1,row,Ts,Ta,Pr_sno,SWE,...
    SWEtm1,Dtm1,rostm1,ros_max1,ros_max2);

if SWE > tiny && ros > 0
    SP_wc_max = SWE*(0.03*(ros >= 200) + ...
        (ros < 200)*(0.03 + 0.07*(200-ros)/200));
else
    SP_wc_max = 0;
end

Pr_liq_snow = Pr_liq*(1-(1-Csno)*Cwat - ...
    sum(Ccrown.*Cfol_H)*(1-dw_SNO));

SP_wc = max(0,SP_wc + Sm + max(0,Pr_liq_snow));
WR_SP = max(0,SP_wc - SP_wc_max);
SP_wc = SP_wc - WR_SP;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% Ice depth, icepack water storage and outflow
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

ICE_D = 0.001*ICE*row/roi;
IP_wc_max = max(0,Ice_wc_sp*ICE);

KT = max(ICE/1000,0.001);
IP_outflow = min(IP_wc,IP_wc/KT*dth);

Pr_liq_ice = Pr_liq*(1-Csno)*(1-(1-Cicew)*Cwat - sum(Ccrown.*Cfol_H));
IP_wc = max(0,IP_wc + WR_SP - IP_outflow + max(0,Pr_liq_ice));

WR_IP_fill = max(0,IP_wc - IP_wc_max);
IP_wc = IP_wc - WR_IP_fill;
WR_IP = IP_outflow + WR_IP_fill;
end



%%%%%% CANOPY SNOW INTERCEPTION
%%REFERENCES [Hedstrom and Pomeroy 1998; Pomeroy et al., 1998 --- Gelfan et al 2004 ---
%Pomeroy et al 2002 --- Liston and Elder 2006 ]
%%%% Parameter --- Hedstrom and Pomeroy 1998;
function[U_SWE,In_SWE,NIn_SWE]=Snow_Interception(dth,Pr_sno,Tdew,Ws,t_sls,ESN_In,In_max_SWE,In_SWEtm1)
Cp=1;  %%
%a=1.157e-006;%% [1/s]
%ah=a*3600;%% [1/h]
ah=4.1e-3;%% [1/h]
if In_max_SWE <= 0
    U_SWE = 0;
    In_SWE =0;
    NIn_SWE =0;
else
    NIn_SWE = 0.7*(In_max_SWE - In_SWEtm1)*(1- exp(-Cp*Pr_sno/In_max_SWE)); %%% New interepted snowpack  [mm]
    In_SWE = In_SWEtm1 + NIn_SWE -ESN_In*dth; %%% First updated Interception  [mm]
    %%%%%
    if (Tdew > 0 && Ws > 0.5 ) %%%  Gelfan et al 2004
        U_SWE = In_SWE; %% Unload Intercepted Snowpack [mm]
        In_SWE = 0 ; %% Intercepted snowpack [mm]
    else
        %U_SWE = (1-exp(-a*dth))*exp(-a*t_sls)*In_SWE; %%% Unload Intercepted Snowpack [mm]
        %U_SWE = -((1-dth*ah/2)/(1+ah*dth/2)-1)*In_SWE;
        U_SWE = (1-exp(-ah*dth))*In_SWE;
        In_SWE = In_SWE - U_SWE; % Intercepted snowpack [mm]
        if In_SWE <= 0.01  %% [mm]
            U_SWE = U_SWE + In_SWE ;
            In_SWE = 0 ;
        end
    end
end
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%% SNOW DEPTH AND DENSITY FUNCTION
% References (Bras, 1990)  [Essery et al., 1999]
%[Hedstrom and Pomeray 1998]  [Douville et al 1995]  [Belair et al 2003]
function[D,ros]=Snow_Depth(dt,Sm,row,Ts,Ta,Pr_sno,SWE,SWEtm1,Dtm1,rostm1,ros_max1,ros_max2)
%%% INPUTS
% dt = [s]
% Sm snow melting
%row % water density [kg/m^3]
%Ta  % air temperature [�C]
% Pr_sno [mm]
% SWEtm1 [mm]
% Dtm1 [m]
% rostm1 [kg/m^3] density of snow
%%% OUTPUTS
% D [m] snow depth
% ros [kg/m^3] snow density
ANS = 2;
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
switch ANS
    case 1 %%%%%%%%% Compaction due to snowfall and New Snow mixing
        %ros_max = 1000; %% [kg/m^3]
        ros_max = 450; %% [kg/m^3] 1000
        if (SWEtm1 > 0) || (Pr_sno>0)
            %%PARAMETERS
            %tau_f = 0.24; %%
            %tau_f = 0.0072; %%
            %tau_1 = 86400; %% [s]
            %ros_n = 67.92 + 51.25*exp(Ta/2.59); %%  [Hedstrom and Pomeray 1998]
            ros_n = 1000*(0.05 + (((9/5)*Ta +32)> 0).*(((9/5)*Ta+32)./100).^2); %% [kg/m^3] new snow density [from Bras 1990]
            if rostm1 == 0
                Dn = 0.001*SWE*(row/ros_n); %% New Snowpack [m]
            else
                Dn = 0.001*Pr_sno*(row/ros_n); %% New Snowpack [m]
            end
            if rostm1 > 0
                %%% Snow Compaction ??Temperature  %%%% Douville et al 1995
                %if Ts == 0
                %    ros =  ros_max -(ros_max-rostm1)*exp(-tau_f*dt/tau_1);
                %else
                ros =rostm1 ;
                %end
                %%% Snow Compaction ??Temperature
                D = 0.001*row*SWEtm1/ros; %% First Update Snow Depth [m]
                %%% Anderson and Crawford, (1964)
                del_D = 0.0254*( ((Pr_sno/25.4)*(1000*Dtm1/25.4))/(SWEtm1/25.4))*((1000*Dtm1/254)^0.35) ; %%% [m]  Snow compaction new snow
                %%% Approximation D = Dtm1 and SWEtm1 = SWE
                D_sm = 0.001*row*Sm/ros; %% [m] melting snow
                D = D - del_D + Dn - D_sm; %% [m] new snow depth
                ros = row*(SWE/(D*1000)); %% New density of snow [kg/m^3]
            else
                D= Dn ;% [m] new snow depth
                ros = ros_n; %% [kg/m^3]
            end
            %%%%%%%
        else
            D = 0;
            ros = 0;
        end
        %%% Final Check
        D(D<0) = 0;  ros(ros < 0) = 0;
        D(ros > ros_max) =  0.001*row*SWE/ros_max;
        ros(ros > ros_max)= ros_max;
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    case 2 %%% Belair 2003
        %%%%%%%%%%%%%%%
        %ros_max1=520; %600; %%% [kg/m^3]
        %ros_max2=320; %450; %%% [kg/m^3]
        if Sm >= 0.1 %%[mm]
            ros_max =( 1000/row)*(ros_max1- 20.47*(1-exp(-Dtm1/0.0673))/Dtm1); % [kg/m^3]
        else
            ros_max =( 1000/row)*(ros_max2- 20.47*(1-exp(-Dtm1/0.0673))/Dtm1);% [kg/m^3]
        end
        if (SWEtm1 > 0) || (Pr_sno>0)
            %%%PARAMETERS
            tau_f = 0.24; %%
            tau_1 = 86400; %% [s]
            %ros_n = 67.92 + 51.25*exp(Ta/2.59); %%  [Hedstrom and Pomeray 1998]
            ros_n = 1000*(0.05 + (((9/5)*Ta +32)> 0).*(((9/5)*Ta+32)./100).^2); %% [kg/m^3] new snow density [from Bras 1990]
            ros_n(ros_n > 158.5) = 158.5; % Added MM 04/04/2024
            %ros_n = (1000/row)*(109+ (6*Ta) +26*sqrt(ua)); %% Brun et al 1989
            if rostm1 == 0
                Dn = 0.001*SWE*(row/ros_n); %% New Snowpack [m]
            else
                Dn = 0.001*Pr_sno*(row/ros_n); %% New Snowpack [m]
            end
            if rostm1 > 0
                %%% Snow Compaction
                if rostm1 < ros_max
                    ros1 =  ros_max -(ros_max-rostm1)*exp(-tau_f*dt/tau_1);
                else
                    ros1 = rostm1 ;
                end
                %%% Snow mixing with new snow
                ros = (ros_n*Pr_sno + SWEtm1*ros1)/(Pr_sno+SWEtm1); %% [kg/m^3]
                D= 0.001*SWE*(row/ros); %% Snowpack [m]
                %%%%%%
            else
                D= Dn ;% [m] new snow depth
                ros = ros_n; %% [kg/m^3]
            end
            %%%%%%%
        else
            D = 0;
            ros = 0;
        end
        %%% Final Check
        D(D<0) = 0;  ros(ros < 0) = 0;
end
end

