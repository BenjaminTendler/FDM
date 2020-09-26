function FDM = FDM_implementation(phase,mag,mask,TE)
%This function provides a practical implementation for FDM as described in:
%Tendler, Benjamin C., and Richard Bowtell. "Frequency difference mapping applied to the corpus callosum at 7T." Magnetic resonance in medicine 81.5 (2019): 3017-3031.

%This code has been written for 2D images (as used in the paper), but can be subsequently modified for 3D

%Written by Benjamin Tendler - contact benjamin.tendler@ndcn.ox.ac.uk

%%
%Output corresponds to the frequency difference map

%Inputs are as follows
%phase  Multi echo gradient echo phase
%mag    Multi echo gradient echo magnitude
%mask   Mask corresponding to each echo (used to correct the linear phase variation along read direction)
%TE     Echo times (in seconds)

%%
%Eq. [3]

for k=2:size(phase,3)
    div_phase(:,:,k-1)=angle(exp(1i*phase(:,:,k))./(exp(1i*phase(:,:,1))));
end
 
%%
%Eq. [4]

for k=1:size(div_phase,3)
    FD(:,:,k)=angle(exp(1i*div_phase(:,:,k))./(exp(1i*(k)*div_phase(:,:,1))));
end
 
%%
%Correction for linear phase variation in read direction - Eq. [8]

for k=1:size(FD,3)
    %Multiply magnitude by exponential of phase data and mask
    z=mask(:,:,k).*mag(:,:,k).*exp(1i*FD(:,:,k));
    %Sum along axis to get 1d vector which will display wrap (you may need to change this and subsequent code depending on axis)
    z_sum=sum(z,1)./sum(mask(:,:,k),1);
    %Remove nans
    z_sum(isnan(z_sum)==1)=0;
    %Unwrap 1d line
    pz=unwrap(angle(z_sum));
    %Generate matrix to calculate the gradient and offset of our unwrapped 1D line 
    a=[ones(size(z,2),1),(-size(z,2)/2:size(z,2)/2-1)'];
    %Generate magnitude vector
    mz=abs(z_sum);
    %Use lscov algorithm to obtain fitting variables, additionally attaching weights with the magnitude data
    ls_fit=lscov(a,pz',mz');
    %Generate line from ls_fit variables
    un_vec=ls_fit(2)*(-size(z,2)/2:size(z,2)/2-1)+ls_fit(1);
    %Replicate along y
    un_vec_mat=repmat(un_vec,[size(FD,2),1]);
    %Apply to matrix
    FD_corr(:,:,k)=angle(exp(1i*FD(:,:,k))./exp(1i*un_vec_mat));
end
 
%%
%Eq. [5]

for k=2:size(FD_corr,3)
    FDM(:,:,k)=FD_corr(:,:,k)/(2*pi*(k-1)*(TE(2)-TE(1)));    
end
