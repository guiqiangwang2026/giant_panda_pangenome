# build a graph-pangenome
minigraph --inv no -xggs -t {threads} assembly/sample1.fa assembly/sample2.fa > graph/{asb}_graph.gfa

# extract node information
awk '$1~/S/ {{ split($5,chr,":"); split($6,pos,":"); split($7,arr,":"); print $2,length($3),chr[3],pos[3],arr[3] }}' {asb}_graph.gfa > {graph/{asb}_graph_len.tsv

# extract edge(link) information
awk '$1 == "L"' {asb}_graph.gfa > graph/{asb}_graph_link.tsv
