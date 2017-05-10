from old_code.util import run_cmd

dir = '/data1/shunfu1/SNPCalling/data_real_HeteSNPsCodingRegion/'
snp = '%s/chr15_snp_m_codingRegion_noHomo.txt'%dir
read_dir_n = 'reads_N10000000_L100_Err0.00'
out_fn = 'snp_read_cov.txt'

cmd = 'python sim_data_generator.py --snp_read_cov1 '+\
      '-s %s '%snp+\
      '-m %s/%s/intermediate/reads_m.bed '%(dir, read_dir_n)+\
      '-p %s/%s/intermediate/reads_p.bed '%(dir, read_dir_n)+\
      '-o %s/%s/intermediate/%s'%(dir, read_dir_n, out_fn)
run_cmd(cmd)