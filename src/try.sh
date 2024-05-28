
    ./../bin/vgan trailmix -t 60 -fq1 test.fq --iter 1000000 --burnin 10000 -k 1 --chains 5 --depth -1 \
                           -o downsampling_results/Q \
                           --deam3p ../share/damageProfiles/dmid3p.prof --deam5p ../share/damageProfiles/dmid5p.prof
