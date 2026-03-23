echo -e "filename\tseed" > 02seed_pleiotrpic_effects.log

for ispleiotropy in nopleio
do
  for num_snp in 10000
  do
    for per_lab in 0.001 ##10 SNPs contributing to the lab, 4995 SNPs constributing to the 23 and 2818 separately
    do
      for set in {1..100}
      do
      echo "slim -t -m -d \"ispleiotropy='${ispleiotropy}'\" -d \"num_snp='${num_snp}'\" -d \"set='${set}'\" -d \"per_lab='${per_lab}'\" 02_add_pleiotropic_effect.slim"
      done
    done
 done
done
