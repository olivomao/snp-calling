from old_code.final_caller19thjuly_m import process_count_line
import pdb

with open('tmp/sim_snps_err_reads_fp_cases.txt', 'r') as f:
    for line in f:
        tokens = line.split()
        if len(tokens)<8: continue
        print(line)
        print(process_count_line(tokens))
        pdb.set_trace()