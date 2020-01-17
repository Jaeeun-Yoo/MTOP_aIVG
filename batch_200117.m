clear
close all
clc


%% '20.01.17

% 2D HMBB, 67,500 elem, DLratio=5, R=5.0 elem, 1.0 FE elem
MTOP_aIVG_R10_2s(192,64,0,2.0,1,4.0,0.5,1,5,3,3,'QUAD4','HMBB','Direct',0) % 2D, 96x32 FE / 480x160 density, DLratio=5, aIVG
MTOP_aIVG_R10_2s(192,64,0,2.0,1,4.0,0.5,0,5,3,3,'QUAD4','HMBB','Direct',0) % 2D, 96x32 FE / 480x160 density, DLratio=5
MTOP_aIVG_R10_2s(192,64,0,0.4,1,4.0,0.5,0,1,3,3,'QUAD4','HMBB','Direct',0) % 2D, 480x160 FE / 480x160 density, DLratio=1
