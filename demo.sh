#!/bin/bash

echo "cybase"

echo "tt..."
./tt_BBEd_floats cybasebr_s-sr032_f-cc.floats > cybasebr_s-sr032_f-cc_floats_BBEd_levenshtein.tt

echo "pivots..."
python ba_pivots_klaesa1_BBEd.py cybase cybase

for k in 1 5 10 20
do

echo $k

echo "kreal..."
python ba_kreal_BBEd_k.py $k cybase cybase &> ba_kreal_BBEd_k.log.txt

echo "kaesa..."
python ba_kaesa_BBEd_k.py $k cybase cybase &> ba_kaesa_BBEd_k.log.txt

echo "laesa..."
python ba_klaesa1_pivots_BBEd_k.py $k cybase cybase &> ba_klaesa1_pivots_BBEd_k.log.txt

echo "laesaeabbedext..."
python ba_kblaesa1BBEdExt_pivots_BBEd_k.py $k cybase cybase &> ba_kblaesa1BBEdExt_pivots_BBEd_k.log.txt

echo "laesaeabubu..."
python ba_kblaesa1BBEdExt1bubu_pivots_BBEd_k.py $k cybase cybase &> ba_kblaesa1BBEdExt1bubu_pivots_BBEd_k.log.txt

done

echo "bye."

