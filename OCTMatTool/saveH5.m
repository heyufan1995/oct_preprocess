function saveH5(savename,Imh5,dataname,inv)
if nargin<4
    inv=true;
end
delete savename;
if inv
    sizes = flip(size(Imh5));
    Imh5 = permute(Imh5,length(sizes):-1:1);
end
h5create(savename,['/',dataname],sizes);
h5write(savename,['/',dataname],Imh5);
sizes