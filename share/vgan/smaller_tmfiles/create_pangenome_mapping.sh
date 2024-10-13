odgi_path="/home/ctools/odgi/bin/odgi"

for NODE in $(seq 1 1 6355)
do
$odgi_path position -i pub.graph.odgi -g $NODE -r "HV0c" >> pangenome_mapping
echo "Done with"
echo $NODE
done
