python ../scripts/pydeseq2_btest_permutation_microbe.py -x table_sim_microbes.txt \
                                             -y metadata_sim.txt \
                                             -z 'disease' \
                                             -k 2 \
                                             -p 2 \
                                             -G G_matrix_sim.txt \
                                             -pt 0.05 \
                                             -o test_output
