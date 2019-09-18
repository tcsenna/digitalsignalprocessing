%SETPATH: program to set path for SeismicLab and Demos 
 
% --------------------------------------------------------
% You need to modify Dir and DirD according to your system  
% --------------------------------------------------------
 
% modify YYY and NNN
% or ask Vicente...  vicente@ualberta.ca 

 Dir1='C:\Users\Kenji\Dropbox\Kenji-DROPBOX\MATLAB\SeismicMatlab\SeismicLab\codes\';
 Dir2='C:\Users\Kenji\Dropbox\Kenji-DROPBOX\MATLAB\SeismicMatlab\SeismicLab\';
 Dir3='C:\Users\Kenji\Dropbox\Kenji-DROPBOX\MATLAB\Kenji_Desconvolucao\Codigos';
 Dir4='C:\Users\Kenji\Dropbox\Kenji-DROPBOX\MATLAB\Kenji_Desconvolucao\gipsa-lab';
 
 path(path, strcat(Dir1,'bp_filter'));
 path(path, strcat(Dir1,'decon'));
 path(path, strcat(Dir1,'dephasing'));
 path(path, strcat(Dir1,'fx'));
 path(path, strcat(Dir1,'interpolation'));
 path(path, strcat(Dir1,'kl_transform'));
 path(path, strcat(Dir1,'radon_transforms'));
 path(path, strcat(Dir1,'scaling_tapering'));
 path(path, strcat(Dir1,'segy'));
 path(path, strcat(Dir1,'seismic_plots'));
 path(path, strcat(Dir1,'synthetics'));
 path(path, strcat(Dir1,'velan_nmo'));
 path(path, strcat(Dir1,'spectra'));

 path(path, Dir2);
 path(path, strcat(Dir2,'SeismicLab_demos'));
 path(path, strcat(Dir2,'SeismicLab_data'));
 
 path(path, Dir3);
 path(path, Dir4);
