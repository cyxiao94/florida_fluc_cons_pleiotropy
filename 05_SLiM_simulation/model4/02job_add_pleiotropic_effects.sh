echo -e "filename\tseed" > 02seed_pleiotrpic_effects.log

for ispleiotropy in pleio
do
  for num_snp in 10000
  do
    for set in {1..100}
    do
    echo "slim -t -m -d \"ispleiotropy='${ispleiotropy}'\" -d \"num_snp='${num_snp}'\" -d \"set='${set}'\" 02_add_pleiotropic_effect.slim"
    done
 done
done
