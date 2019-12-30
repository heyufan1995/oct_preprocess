function oct_vol_norm = normPeaks2(oct_vol,header,bds,peak_vals_out)
% normalize the mean value and standard deviation of the  retina image in
% each Bscan to constants
retina_mask = convertBoundariesToLabels(bds(:,:,[1 end]),size(oct_vol),true);
oct_vol_norm = oct_vol;

for ii = 1:size(oct_vol,3)
    bscan = oct_vol_norm(:,:,ii);
    bscan_mask = find(retina_mask(:,:,ii)==1);
    retina = bscan(bscan_mask);
    a=0.15/std(retina);
    b=0.37-a*mean(retina);
    oct_vol_norm(:,:,ii)=a*oct_vol_norm(:,:,ii)+b;
end
oct_vol_norm(oct_vol_norm<0) = 0;
oct_vol_norm(oct_vol_norm>1) = 1;