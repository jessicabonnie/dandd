#!/bin/bash
ml anaconda
conda activate ~/miniconda3/envs/dashing

species=$1
approach=$2
kval=$3
upperout=$4

datadir=/home/jbonnie1/scr16_blangme2/jessica/data/${species}
outdir=${upperout}/${species}
apout=${outdir}/${approach}
dashing=/home/jbonnie1/lib/dashing/dashing
outprefix=${outdir}_${approach}_k${kval}

mkdir -p ${apout}
cd ${datadir}
filelist=$(cd ${datadir} && ls *gz)
# rm ${outdir}_${approach}.out
for fasta in ${filelist[@]}; do
echo Now acting on $fasta
if [ $approach == 'kmc' ]; then
/usr/bin/time -o ${outprefix}.out -v kmc -ci1 -cx15 -k ${kval} -fm ${datadir}/${fasta} ${apout}/${fasta}.kmc ${outdir}/kmc
elif [ $approach == 'dashing' ]; then
/usr/bin/time -o ${outprefix}.out -v ${dashing} sketch -k ${kval} -S 20 ${fasta} -P ${apout}
fi
done


fastas=($(ls ${datadir}/*gz))
nfasta=${#fastas[@]}

get_seeded_random()
{
  seed="$1"
  openssl enc -aes-256-ctr -pass pass:"$seed" -nosalt \
    </dev/zero 2>/dev/null
}

alias 'ubtime=/usr/bin/time '

union_and_count()
{
  prefix=$1
  karg=$2
  ${dashing} union -p 10 -z -o ${prefix}_k${karg}.hll -F ${prefix}.txt
  ${dashing} card -p10 --presketched ${prefix}_k${karg}.hll > ${prefix}_k${karg}.card
  
}
export -f union_and_count

rand_ind=( $(shuf -i1-${nfasta} --random-source=<(get_seeded_random 42)) )
declare -p rand_ind

resorted=( $(for i in ${rand_ind[@]}; do echo ${fastas[$i]}; done) )
sketched=( $(for i in ${resorted[@]}; do echo ${apout}/$(basename $i).w.${kval}.spacing.20.hll; done ))
echo ${resorted[@]}
#for howmany in {1..10}; do
for howmany in $(seq $nfasta); do
echo $howmany
#rm ${outprefix}_c${howmany}.out

if [ $approach == 'kmc' ]; then
echo ${resorted[@]:0:${howmany}} | sed 's/ /\n/g' > ${apout}/c${howmany}.txt
/usr/bin/time -o ${outprefix}_c${howmany}.out -v kmc -ci1 -cx15 -k $kval -fm @${apout}/c${howmany}.txt ${apout}/c${howmany}_k${kval}.kmc ${outdir}/kmc
elif [ $approach == 'dashing' ]; then
echo ${sketched[@]:0:${howmany}} | sed 's/ /\n/g' > ${apout}/c${howmany}.txt
#/usr/bin/time -v ${dashing} union -p 10 -z -o ${apout}/c${howmany}.hll -F ${apout}/c${howmany}.txt >> ${outprefix}_c${howmany}.out 2>&1
#/usr/bin/time -v union_and_count ${apout}/c${howmany} $kval >> ${outprefix}_c${howmany}.out 2>&1
/usr/bin/time -o ${outprefix}_c${howmany}.out -v ${dashing} union -p 10 -z -o ${apout}/c${howmany}_k${kval}.hll -F ${apout}/c${howmany}.txt
/usr/bin/time -a -o ${outprefix}_c${howmany}.out ${dashing} card -p 10 --presketched ${apout}/c${howmany}_k${kval}.hll > ${apout}/c${howmany}_k${kval}.card

fi

# if [ $approach == 'kmc' ]; then
# /usr/bin/time -v kmc -ci1 -cx15 -fm @${apout}/c${howmany}.txt ${apout}/k${kval}_c${howmany}.kmc ${outdir}/kmc >> ${outdir}_${approach}_c${howmany}.out 2>&1
# elif [ $approach == 'dashing' ]; then
# /usr/bin/time -v ${dashing} union -p 10 -z -o ${apout}/c${howmany}.hll -F ${apout}/c${howmany}.txt >> ${outdir}_${approach}_c${howmany}.out 2>&1
# fi

done
