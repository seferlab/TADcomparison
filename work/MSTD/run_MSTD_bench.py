from MSTD_v1 import MSTD
import os



chrs=list(range(1,101))
chrs = [str(x) for x in chrs]


rep='bench'
print(rep)
os.makedirs('DOMAINS/'+rep)
for chrom in chrs:
	print(chrom)
	MSTD('../benchmarks/maptest_'+chrom+'.txt', 'DOMAINS/'+rep+'/'+rep+'_'+chrom+'.txt', MDHD=10,symmetry=1,window=10,visualization=0)

