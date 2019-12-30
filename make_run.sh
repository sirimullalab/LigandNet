for csv in ./pharos_database/actives_fingerprints/*.csv
do
	echo "python train.py --actives $csv"
done >> run.sh
