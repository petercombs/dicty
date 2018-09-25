addpath mex util


t = readtable('test2/strain_genotypes.tsv', 'FileType', 'text');
a = readtable('test2/strain_abundances.tsv', 'FileType', 'text');

T = table2array(t);
T = T(2:end, 2:end);
[n_snps, n_strains] = size(T)

A = table2array(a);
A = A(2:end, 2:end);

D = T * A';

opt0 = opt_Integerfac_findvert('nsamples', 11, 'aggr', 'bestsingle', 'affine', true);

[That, Ahat, status] = Integerfac_findvert_cpp(D, n_strains, [0, 1], opt0);

Dhat = That * Ahat;
diff = D - Dhat;
disp('RMS')
disp(rms(diff(:)))
disp('avg amt')
disp(mean(D(:)))


ahatt = table(Ahat);
thatt = table(That);

writetable(ahatt);
writetable(thatt);
