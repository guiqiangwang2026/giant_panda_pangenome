### 1. Graph-based pan-genome

```shell
# build pangenome
bash 01.minigraph_build_Pangeno.sh

# pan and core node analysis
Rscriipt 02.Pan_core.R
```

### 2. SVs calling

```shell
# extract sv
gfatools bubble panda_graph.gfa > bubble/panda_bubble.tsv
```

### 3. SVs Statistics

```shell
# Biallelic type bubble
awk '$5==2 {{ print $1,$2,$4,$5,$12 }}' panda_bubble.tsv > panda_biallelic_bubble.tsv

# Biallelic SV
scripts/03.get_bialsv.py -a graph/{asb}_graph_len.tsv -a analysis/bubble/{asb}_biallelic_bubble.tsv  > panda_biallelic_sv.tsv

# multiallelic type bubble
awk '$5>2 && $5 < 8 {{ print $1,$2,$4,$5,$12 }}' panda_bubble.tsv > panda_multiallelic_bubble.tsv

# multiallelic SV
scripts/04.get_multisv.py -a graph/{asb}_graph_len.tsv -a analysis/bubble/{asb}_multiallelic_bubble.tsv -a analysis/colour_node/{asb}_nodecol.tsv > panda_multiallelic_sv.tsv

# Convert panda_bubble.tsv to bed
awk '{{ print $1,$2,$2+1,$1"_"$2 }}' OFS="\t" panda_bubble.tsv > panda_bubble.bed
```

### 4. SVs annotation

#### 4.1 Prepare files

```shell
# After perfomed panda-graphs
mkdir /analysis/bubble/01.sv_analysis/03.annotation
cd /analysis/bubble/01.sv_analysis/03.annotation
cp ../../../../graph/*_graph.gfa ./
cp ../../../../assembly/panda_annota.gff ./
cp ../../*_path_trace.tsv panda_path_trace.tsv.ori
cp TE.gff ./rep.gff
```

#### 4.2 Pre-process files

```shell
cat panda_path_trace.tsv.ori  | perl -alne '$F[2]="$F[2]_$F[0]"; print join "\t", @F' > panda_path_trace.tsv

cat panda_graph.gfa | grep '^S' | perl -alne 'splice @F,2,1; print join "\t", @F' > thin.gfa

perl 05.prase_gff.pl panda_annota.gff panda_annota.gff.gene

perl 06.prase_gff.stat.pl panda_annota.gff.gene panda_annota.gff.gene.stat
```

#### 4.3 Compute SV location (exon / intron / intergenetic)

```shell
# SV location information
perl 07.trace_add_pos.pl thin.gfa panda_path_trace.tsv 1.trace_add_pos.pl.out

# Calculate the length of reflen and compare it with altlen, then output the larger of the two values
cat 1.trace_add_pos.pl.out | perl -alne 'BEGIN{print join "\t", qw/svid type chr start end altlen reflen maxlen/} $reflen = $F[4]-$F[3]; $maxlen=$reflen>$F[5] ? $reflen : $F[5]; print join "\t", @F, $reflen, $maxlen' > 1.trace_add_pos.pl.out.stat

# Retain variants with a length > 50
cat 1.trace_add_pos.pl.out.stat | perl -alne 'print if $.==1 or $F[-1]>50'  > 1.trace_add_pos.pl.out.stat.svonly
cat 1.trace_add_pos.pl.out | perl -alne 'print join "\t", @F[2,3], $F[4]+1, $F[0], $F[1], @F[3,4,5]' > 2.bed

# The intersection of SV and genes
bedtools intersect -a 2.bed -b panda_annota.gff.gene -wb  > 3.bed.gff4
Rscript 08.R
```

#### 4.4 Compute SV location (5k flank region)

```shell
# Obtain 5KB flank regions of the gene
bash 09.prase_flank.sh
sed -i “s/ /\t/” 3.5kbflk.bed

# The intersection of SV and 5k flank region
bedtools intersect -a 2.bed -b 3.5kbflk.bed -wb > 5.5kbflk.bed
```

