from MSTD_v1 import MSTD
import os



replicate=['simulate1','simulate2']

noises=['0.04','0.08','0.12','0.16','0.20']



for rep in replicate:
    print(rep)
    os.makedirs('DOMAINS/'+rep)
    for noise in noises:
        print(noise)
        MSTD('../simulate_data/'+rep+'/sim_ob_'+noise+'.chr5', 'DOMAINS/'+rep+'/'+rep+'_'+noise+'.chr5', MDHD=10,symmetry=1,window=10,visualization=0)


