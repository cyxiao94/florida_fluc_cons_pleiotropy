echo -e "filename\tseed" > 02seed_pleiotrpic_effects.log
for ispleiotropy in nopleio
do
  for num_snp in 10000
  do
    for per_lab in 0.3333 ##approx. equal number of SNPs that constribtuing to the lab, 23, 2818 separately
    do
      for set in {1..100}
      do
      echo "slim -t -m -d \"ispleiotropy='${ispleiotropy}'\" -d \"num_snp='${num_snp}'\" -d \"set='${set}'\" -d \"per_lab='${per_lab}'\" 02_add_pleiotropic_effect.slim"
      done
    done
 done
done
