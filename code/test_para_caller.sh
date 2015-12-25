#split
do_split=false
if $do_split
then
  count_dir=/home/sreeramkannan/singleCell/SNP_Calling_Summer15/data_0827_SNP1k_Reads10M/tophat_out/
  count_fn=count.txt
  count_split_dir=/count_split/

  cd $count_dir
  mkdir -p $count_split_dir

  split -l 167880 -d $count_fn
  mv x* $count_split_dir
fi

#para_caller
do_para_caller=true
if $do_para_caller
then
  parallel python final_caller19thjuly_m.py x{} caller_op_x{}.txt caller_op_exception_x{}.txt caller_op_snp_x{}.txt ::: 08 10 18 #00 01 02 03 04 05 06 07 08 09 10 11 12 13 14 15 16 17 18 19
#python final_caller19thjuly_m.py x08 caller_op_x08.txt caller_op_exception_x08.txt caller_op_snp_x08.txt
fi
