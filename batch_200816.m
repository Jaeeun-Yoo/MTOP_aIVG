clear
close all
clc


%% '20.08.16

% 2D HMBB, 67,500 elem, DLratio=5, R=5.0 elem, 1.0 FE elem
MTOP_aIVG_R11_4s(192,64,0,1.0,1,2.0,0.5,1,5,3,3,'QUAD4','HMBB','Direct',0) % 2D, 192x64 FE / 960x320 density, DLratio=5, aIVG
MTOP_aIVG_R11_4s(192,64,0,1.0,1,2.0,0.5,0,5,3,3,'QUAD4','HMBB','Direct',0) % 2D, 192x64 FE / 960x320 density, DLratio=5
MTOP_aIVG_R11_4s(192,64,0,0.2,1,2.0,0.5,0,5,3,3,'QUAD4','HMBB','Direct',0) % 2D, 960x320 FE / 960x320 density, DLratio=1


