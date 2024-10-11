function filtered_SA = SA_filtering(SA_sig,RFs,GM,num_sidebands,particular_ords)
% DESCRIPTION:
% This function roughly filteres undesired orders from the SA signal by
% setting a zero value to these orders in the spectrum.
% There are mainly three uses for this function:
% 1.) Building the RESIDUAL signal (eliminating GM harmonics).
% 2.) Building the DIFFERNECE signal (eliminating GM harmonics + AM SB's).
% 3.) Eliminating particular orders, on demand (e.g., shaft harmonics).
% =====
% INPUTS:
% * SA_sig - The SA signal in the cycle domain.
% * RFs - Sampling rate after angular resampling, equals to the SA length.
% * GM - The fundamental gear mesh order. In case GM eliminating is not
% necessary put NaN (GM=[])!
% * num_sidebands - The number of SB's around each side of the GM harmonics
% to eliminate. usually used for the difference calculation to eliminate
% the close pair of AM sidebands. In case sidebands eliminating is not
% necessary, put NaN (num_sidebands=[])!
% * particular_ords - A row vector consists of the integer orders to 
% eliminate. for example: particular_ords=[2,3,6,7]. In case eliminating
% particular orders is uneccesary, put NaN (particular_ords=[]).
% =====
% OUTPUTS:
% * filtered_SA - The SA signal in the cycle domain after filtering.
% =====
% IN FUNCTION VARIABLES:
% * ords2delete - A vector consisting of the desired orders to delete.
% * max_order - The maximum order within theoretical bandwidth.
% * num_of_GM_harms - Number of GM harmonics available within bandwidth.
% * GM_harms_inds - A vector with the indices of GM harmonics orders.
% * GM_harm_ord - The order of the ii-th GM harmonic.
% * ords2delete_positive_inds / ords2delete_negative_inds - The indices 
% of the positive & negative orders to eliminate. recall that the FFT
% algorithm distributes the energy of each order symetrically between its
% positive and negative values. For more information, assist help on 'fft'.
% * ords2delete_inds - The total indices representing the orders to
% eliminate from the SA signal, sorted by index.
% * fft_SA - The Fourier Transform of the fft signal with length of RFs.
% =====
% Created by Ph.D. student Lior Bachar, Gears Team, 2021.
% PHM Laboratory, Ben-Gurion University of the Negev, Be'er Sheva, Israel.
% Email: liorbac@post.bgu.ac.il
%%

ords2delete = particular_ords ; % initialize ords2delete with the particular orders and then calculate the rest. 
max_order = floor(RFs/2) ; % the theoretical (!) bandwidth. 

if GM % in case eliminating the GM harmonics is desired
    num_of_GM_harms = floor(max_order./GM) ;
    if num_sidebands % in case eliminating the GM harmonics and its surrounding sidebands is desired
        for ii = 1:num_of_GM_harms
            GM_harm_ord = GM*ii ;
            ords2delete = [ords2delete,...
                (GM_harm_ord-num_sidebands):(GM_harm_ord+num_sidebands)] ; % each band includes the order of the ii-th GM and its sidebands from both sides.
        end % of for ii=1:num_of_GM_harms
    else % i.e., eliminating only the GM harmonics without sidebands.
        ords2delete = [1:num_of_GM_harms]*GM ;
    end % of if num_sidebands
end % of if GM

ords2delete(ords2delete>max_order) = [] ; % in case some SB's crossed the max order
ords2delete_positive_inds = ords2delete + 1 ; % convert from order to positive index
ords2delete_negative_inds = RFs - ords2delete_positive_inds + 2 ; % convert from order to negative index
ords2delete_inds = sort([ords2delete_positive_inds,ords2delete_negative_inds]) ;

fft_SA = fft(SA_sig,RFs) ; % fft with resolution of 1 order
fft_SA(ords2delete_inds) = 0 ; % eliminate the requested orders
filtered_SA = ifft(fft_SA,RFs) ; % inverse transform to the cycle domain.

end % of function