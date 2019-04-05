library(diptest)
con = file('tmp/spore_freqs.txt', 'r')

out = list()
num_processed = 0
min_pval = 1
lines = readLines(con)
for (i in 1:length(lines)){
    if (i %% 1000 == 0)
    {
        print(c(i, min_pval))
    }

    line = lines[i]
    line = strsplit(line, '\t')[[1]];
    snp = line[1];
    data = as.numeric(line[2:length(line)]);
    dipres = diptest::dip.test(data);
    if (! startsWith(snp, 'DDB0169550') ){next;}
    if (dipres$p.value < min_pval && 0.4 < mean(data) && mean(data) < 0.6)
    {
        min_pval = dipres$p.value;
        best_pval = snp;
        print(snp)
    }
    #if (dipres$p.value < .005){
        #print(c(snp, dipres$p.value));
    #}
    out[snp] = dipres$p.value;
}

close(con)
