#num_snp 10000
#ispleiotropy nopleiotropy
#percentage_lab 0.8 (80% variants are selection targets for lab traits)
#set 1-100
#pop+environment p1-p5(fluc) p6-p10(cons)

echo -e "num_snp\tispleiotropy\tper_lab\tset\tenvir\tpop" > 0parameters

for num_snp in 10000
do
  for ispleiotropy in nopleio
  do
    for per_lab in 0.3333 ##approx. equal number of SNPs that constribtuing to the lab, 23, 2818 separately
    do
      for set in {1..100}
      do
      #fluc
      echo -e "${num_snp}\t${ispleiotropy}\t${per_lab}\t${set}\tfluc\tp1" >> 0parameters
      echo -e "${num_snp}\t${ispleiotropy}\t${per_lab}\t${set}\tfluc\tp2" >> 0parameters
      echo -e "${num_snp}\t${ispleiotropy}\t${per_lab}\t${set}\tfluc\tp3" >> 0parameters
      echo -e "${num_snp}\t${ispleiotropy}\t${per_lab}\t${set}\tfluc\tp4" >> 0parameters
      echo -e "${num_snp}\t${ispleiotropy}\t${per_lab}\t${set}\tfluc\tp5" >> 0parameters
      #cons
      echo -e "${num_snp}\t${ispleiotropy}\t${per_lab}\t${set}\tcons\tp6" >> 0parameters
      echo -e "${num_snp}\t${ispleiotropy}\t${per_lab}\t${set}\tcons\tp7" >> 0parameters
      echo -e "${num_snp}\t${ispleiotropy}\t${per_lab}\t${set}\tcons\tp8" >> 0parameters
      echo -e "${num_snp}\t${ispleiotropy}\t${per_lab}\t${set}\tcons\tp9" >> 0parameters
      echo -e "${num_snp}\t${ispleiotropy}\t${per_lab}\t${set}\tcons\tp10" >> 0parameters
      done
    done
  done
done

