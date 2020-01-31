for sdf in ../../../ZINC_sdf/*/*/*.sdf
do
	echo "python predict.py --sdf $sdf"
done >> run.sh
