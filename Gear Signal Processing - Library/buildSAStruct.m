function sa_struct = buildSAStruct(vib_sig, rps_nom, z, Fs_original, Fs)
%UNTITLED Summary of this function goes here
%   Detailed explanation goes here

%% downsampling %%
if nargin == 5
    vib_sig_original = vib_sig ;
    vib_sig = resample(vib_sig_original, Fs, Fs_original, 10) ;
else
    Fs = Fs_original ;
end

%% angular resampling %%
t_vib = (1/Fs)*[0:(length(vib_sig)-1)]' ;
rps_t = rps_nom * ones(size(t_vib)) ;
d_cyc = rps_nom / Fs ;
RFs = 2^nextpow2(1/d_cyc) ;

phase_sig = cumsum(rps_t)/Fs ;
phase_sig = phase_sig - phase_sig(1) ;
phase_sig_resamp = [0:(1/RFs):(floor(phase_sig(end))-1/RFs)]' ;  % the phase signal in fine resolution.
resamp_t = interp1(phase_sig, t_vib, phase_sig_resamp,'linear') ; % obtain the new time vector (i.e., the non-normalized cycle vector) by interpolation.
resamp_data = interp1(t_vib, vib_sig, resamp_t, 'PCHIP') ;

%% sa calculations %%
cyc_vctr = [0:(1/RFs):(1-1/RFs)]' ;
all_cycles = vec2mat(resamp_data,RFs) ;
sa_sig = mean(all_cycles)' ;
% sa_sig = sa_sig - mean(sa_sig) ;
% res_sig = SA_filtering(sa_sig, RFs, z, [], [1:7]) ;
diff_sig = SA_filtering(sa_sig,RFs, z, 2, [1:7]) ;

%% spectrum calculations %%
max_ord = floor(Fs/2/rps_nom) ;
ord_vctr = [0:max_ord-1]' ;
ord_spec = 2*(abs(fft(sa_sig,RFs))/RFs) ;
ord_spec = ord_spec(1:max_ord) ;

%% output %%
[sa_struct.cyc_vctr, sa_struct.sa, sa_struct.diff, sa_struct.ord_vctr, sa_struct.fft_sa] = ...
    deal(cyc_vctr, sa_sig, diff_sig, ord_vctr, ord_spec) ;

end

