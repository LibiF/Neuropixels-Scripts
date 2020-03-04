up = load('D:\Libi\Category stimuli\Diff_slope\Rise_to_pure_tone160.mat') ;
up = up.stim_up ;
% down = load('C:\Users\owner\Dropbox\Category stimuli\Teaching_stims_1sec\Down_base_12.mat') ;
% down = down.stim_down ;
audiowrite('puretone16.wav',up,500000)
% audiowrite('down12_6.wav',down,500000)