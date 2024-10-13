#!/bin/bash
vg_path="/home/projects/MAAG/Magpie/Magpie/vgan_corrected/dep/vg/bin/vg"

#$vg_path convert -o -g pub.graph.gfa > pub.graph.og
#$vg_path mod -X 6 pub.graph.og -t 60 > new_pub.graph.og
#mv new_pub.graph.og pub.graph.og
#$vg_path view pub.graph.og > pub.graph.gfa
#python circularize_paths.py pub.graph.og graph_circ.og
#mv graph_circ.og pub.graph.og
#$vg_path paths -Lv pub.graph.og > graph_paths
#./sort.sh graph_paths
#$vg_path snarls pub.graph.og > pub.graph.snarls
#$vg_path index -j pub.graph.dist pub.graph.og
#$vg_path index -x pub.graph.xg pub.graph.og
#$vg_path convert -v pub.graph.og -t 60 > pub.graph.vg
#$vg_path gbwt -g pub.graph.giraffe.gbz --gbz-format -G pub.graph.gfa
#$vg_path gbwt -o pub.graph.gbwt -g pub.graph.gg -Z pub.graph.giraffe.gbz
$vg_path minimizer -o pub.graph.min -g pub.graph.gbwt -d pub.graph.dist pub.graph.og -k 20 -w 10
$vg_path rymer -o pub.graph.ry -g pub.graph.gbwt -d pub.graph.dist pub.graph.og -k 20 -w 10
