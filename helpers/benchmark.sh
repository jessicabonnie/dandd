#!/bin/bash
ml anaconda
conda activate ~/miniconda3/envs/dashing

species=$1
approach=$2
kval=$3
upperout=$4

datadir=/data/blangme2/jessica/${species}
outdir=${upperout}/${species}
apout=${outdir}/${approach}
dashing=/home/jbonnie1/mini/dashing/dashing
outprefix=${outdir}_${approach}_k${kval}
summary=${outdir}_${approach}_k${kval}.csv
#nthreads=8

mkdir -p ${apout}
cd ${datadir}
fastalist=$(cd ${datadir} && ls *gz)
# rm ${outdir}_${approach}.out
#create header of output table
echo "Method,Stage,Fasta,Sketch,Order,Count,Card,MaxResSetSize_kb, WallClock_hms,SystemTime_sec,UserTime_sec" > ${summary}

index=1
## For each fasta create a kmc database or sketch
echo "NOW CREATING INDIVIDUAL DB/SKETCHES"
stage=1

for fasta in ${fastalist[@]}; do
  echo Now acting on $fasta
  echo $fasta
  cardloc=${apout}/${fasta}.card

  if [ $approach == 'kmc' ]; then
    outsketch=${apout}/${fasta}.kmc
    /usr/bin/time -o ${outprefix}.out -v sh -c "kmc -v -k${kval} -ci1 -fm ${datadir}/${fasta} ${apout}/${fasta}.kmc ${outdir}/kmc > ${cardloc}"
    card=$(grep "No. of unique counted k-mers" ${cardloc} | awk '{print $NF}')

  elif [ $approach == 'dashing' ]; then
    outsketch=${apout}/${fasta}.w.${kval}.spacing.20.hll
    /usr/bin/time -o ${outprefix}.out -v sh -c "${dashing} sketch  -k ${kval} -S 20 --prefix ${apout} ${fasta}"
    ${dashing} card  --presketched ${outsketch} > ${cardloc}
    card=$(awk 'NR==2{print $NF}' ${cardloc})
  fi
#capture benchmarking values
mrss=$(grep "Maximum resident" ${outprefix}.out | awk '{print $NF}')
## TODO switch to wallclock
wctime=$(grep "wall clock" ${outprefix}.out | awk '{print $NF}')
systime=$(grep "System time" ${outprefix}.out | awk '{print $NF}')
utime=$(grep "User time" ${outprefix}.out | awk '{print $NF}')

    echo "${approach},${stage},${fasta},${outsketch},$(($index-1)),1,${card},${mrss},${wctime},${systime},${utime}" >> ${summary}
    index+=1
done

# create an array from the fastas and also count what is there
fastas=($(ls ${datadir}/*gz))
nfasta=${#fastas[@]}

get_seeded_random()
#function to set a random seed when shuffle is called
{
  seed="$1"
  openssl enc -aes-256-ctr -pass pass:"$seed" -nosalt \
    </dev/zero 2>/dev/null
}
# 
# sum(){
#   #function to extract and sum one column from a space delimited file
#   filename=$1
#   col=$2
#   cut -d ' ' -f${col} ${filename}| paste -sd+ | bc
# }



## create an array of randomized indices using a predetermined seed
rand_ind=( $(shuf -i0-$((${nfasta}-1)) --random-source=<(get_seeded_random 42)) )
declare -p rand_ind
echo rand indices ${rand_ind[@]}

## resort the fastas, the sketches, and the kmc databases using the randomized indices
resorted=( $(for i in ${rand_ind[@]}; do echo ${fastas[${i}-1]}; done) )

#sketched=( $(for i in ${resorted[@]}; do echo ${apout}/$(basename $i).w.${kval}.spacing.20.hll; done ))
#kmcdbs=( $(for i in ${resorted[@]}; do echo ${apout}/$(basename $i).kmc; done ))

#Depending on which approach, the randomized sketch list will be slightly different due to suffixes
if [ $approach == 'kmc' ]; then 
echo "I'm inside appraoch check"
  sketched=( $(for i in ${resorted[@]}; do echo ${apout}/$(basename ${i}).kmc; done ))
else sketched=( $(for i in ${resorted[@]}; do echo ${apout}/$(basename $i).w.${kval}.spacing.20.hll; done ))
fi
#declare -p sketched

echo ${resorted[@]}
echo ${sketched[@]}
curunion=${sketched[0]}

#for howmany in {1..10}; do

# Make step 1 different from the others ... can be incorporated into the main forloop later
#if [ $approach == 'kmc' ]; then
# capture which database is being used at count 1
#echo ${kmcdbs[0]}  > ${apout}/c1.txt
#curunion=${sketched[0]}
#newunion=${apout}/c1_k${kval}.kmc
#newunion=${apout}/c1_k${kval}.kmc
#cp ${curunion}.kmc_pre ${newunion}.kmc_pre
#cp ${curunion}.kmc_pre ${newunion}.kmc_pre
#tmpdir=${apout}/tmp1
#curunion=${apout}/c1_k${kval}.kmc

#/usr/bin/time -o ${outprefix}_c1.out -v sh -c "kmc -v -ci1 -k${kval} -fkmc ${curunion} ${newunion} ${tmpdir} > ${apout}/c1_k${kval}.card"
# count the unique kmers in the previously built database of c1 
#/usr/bin/time -o ${outprefix}_c1.out -v kmc union -ci1 -k${kval} -fm ${curfasta} ${curunion} ${outdir}/kmc > ${apout}/c1_k${kval}.card
#elif [ $approach == 'dashing' ]; then
#echo ${sketched[0]}  > ${apout}/c1.txt
#curunion=${sketched[0]}
#/usr/bin/time -v -o ${outprefix}_c1.out sh -c "${dashing} card -p 10 --presketched ${curunion} > ${apout}/c1_k${kval}.card"
#fi

## Progressive union rather than a series of global unions which is unfair to kmc
echo "START OF PROGRESSIVE UNION"
stage=2
fsum=NA
#ordering=rand_ind[@]

#echo "Stage,Fasta,Sketch,Order,Count,Card" > ${summary}
for howmany in $(seq 1 $nfasta); do
  echo $howmany
  
  #Things to do at the start of the loop regardless of approach
  these=${sketched[@]:0:${howmany}}
  echo ${these[@]} | sed 's/ /\n/g' > ${apout}/c${howmany}.txt
  #echo ${sketched[@]:0:${howmany}} | sed 's/ /\n/g' > ${apout}/c${howmany}.txt
  cardloc=${apout}/c${howmany}_k${kval}.card
  timeout=${outprefix}_c${howmany}.out

  
  if [ $approach == 'kmc' ]; then
    newunion=${apout}/c${howmany}_k${kval}.kmc
    echo "Random Sketches in New Union:"
    echo ${sketched[@]:0:${howmany}}
    #newhist=${newunion}.hist
    
    if [[ "$howmany" -eq '1' ]]; 
    then
    echo "INSIDE howmany is 1"
      # capture which database is being used at count 1
      echo First Random Sketch: 
      curunion=${sketched[0]}
      
      # using the previous created sketch, get the cardinality for the first db in the random ordering
      cmd="kmc_tools info ${curunion} > ${cardloc}"
      echo ${cmd}
      /usr/bin/time -o ${timeout} -v sh -c "${cmd}"
      #cut -f2 ${newhist} | paste -sd+ | bc > ${cardloc}
      
    else
      echo "INSIDE HOW MANY IS GREATER"
      echo Current Union ${curunion}
      echo New Union ${newunion}
      newguy=${sketched[${howmany}-1]}
      echo Next addition ${newguy}
      
      # benchmark timing for progressive union of sketches
      # TODO : update to use info command 

      cmd="kmc_tools  simple ${curunion} ${newguy} union ${newunion}; kmc_tools info ${newunion}> ${cardloc}" 
      echo ${cmd}
      /usr/bin/time -o ${timeout} -v sh -c "$cmd"
       # kmc_tools  transform ${newunion} histogram ${newhist}
       #    cut -f2 ${newunion}.hist | paste -sd+ | bc > ${cardloc}
      
      #A new union is made when sketch number > 1 so pass that along to the start of the next loop
      curunion=${newunion}
      #/usr/bin/time -o ${outprefix}_c${howmany}.out -v kmc -ci1 -k${kval} -fm @${apout}/c${howmany}.txt ${apout}/c${howmany}_k${kval}.kmc ${outdir}/kmc > ${apout}/c${howmany}_k${kval}.card
      #rm -r ${tmpdir}
    fi
    card=$(awk 'NR==2{print $NF}' ${cardloc})
    echo $card
    
  elif [ $approach == 'dashing' ]; then
    newunion=${apout}/c${howmany}_k${kval}.hll
    if [[ "$howmany" -eq 1 ]]; then
      /usr/bin/time -v -o ${timeout} sh -c "${dashing} card --presketched ${curunion} > ${cardloc}"
    else
    /usr/bin/time -v -o ${timeout} sh -c "${dashing} union -z -o ${newunion} ${curunion} ${sketched[${howmany}-1]}; ${dashing} card --presketched ${newunion} > ${cardloc}"
    curunion=${newunion}
    fi
    card=$(awk 'NR==2{print $NF}' ${cardloc})
  fi
  #capture benchmarking values
  mrss=$(grep "Maximum resident" ${timeout} | awk '{print $NF}')
  wctime=$(grep "wall clock" ${timeout} | awk '{print $NF}')
  systime=$(grep "System time" ${timeout}| awk '{print $NF}')
  utime=$(grep "User time" ${timeout}| awk '{print $NF}')
  if [[ "$howmany" -eq '1' ]]; then wctime=NA; mrss=NA;fi

  echo "${approach},${stage},${fsum},${outsketch},${rand_ind[@]},${howmany},${card},${mrss},${wctime},${systime},${utime}" >> ${summary}


done

## Benchmark a full union creation from sketches
echo "NOW STARTING FULL UNION FROM SKETCHES"
fullunion=${apout}/fullunionq_k${kval}
stage=3
cardloc=${fullunion}.card
timeout=${outprefix}_fullunionq_k${kval}.out
#echo "${approach},${stage},${fastalist},${fullunion},NA,${nfasta},${card}" > ${summary}

if [ $approach == 'kmc' ]; then

  #Create an instruction file for kmc_tools complex
  echo "INPUT:" > ${apout}/complex_union.txt
  echo ${sketched[@]} | sed 's/ /\n/g' | awk '{print "input"NR" = ",$0, "-ci1"}'>> ${apout}/complex_union.txt
  echo "OUTPUT:" >> ${apout}/complex_union.txt
  string="${fullunion} = input1"
  for i in $(seq 2 $nfasta); do string="${string} + input$i"; done
  echo ${string} >> ${apout}/complex_union.txt
  
  #Benchmarch the process

  cmd="kmc_tools complex ${apout}/complex_union.txt; kmc_tools info ${fullunion} > ${cardloc}" 
  echo ${cmd}
  /usr/bin/time -o ${timeout} -v sh -c "${cmd}"
  card=$(awk 'NR==2{print $NF}' ${cardloc})
  #kmc_tools  transform ${fullunion} histogram ${fullunion}.hist; cut -f2 ${fullunion}.hist | paste -sd+ | bc > ${cardloc}
  
elif [ $approach == 'dashing' ]; then
cmd="${dashing} union  -z -o ${fullunion} ${sketched[@]}"
echo ${cmd}
/usr/bin/time -v -o ${timeout} sh -c "${cmd}"

${dashing} card  -S 20 --presketched ${fullunion} > ${cardloc}
  
card=$(awk 'NR==2{print $NF}' ${cardloc})
fi
mrss=$(grep "Maximum resident" ${timeout} | awk '{print $NF}')
wctime=$(grep "wall clock" ${timeout} | awk '{print $NF}')
systime=$(grep "System time" ${timeout} | awk '{print $NF}')
utime=$(grep "User time" ${timeout} | awk '{print $NF}')
echo "${approach},${stage},NA,${fullunion},NA,${nfasta},${card},${mrss},${wctime},${systime},${utime}" >> ${summary}

exit

# Benchmark full union creation from fastas

echo "NOW STARTING FULL UNION FROM FASTAS"
stage=4
echo ${fastas[@]} | sed 's/ /\n/g' > ${apout}/all_fastas.txt
fullunionf=${apout}/fullunionf_k${kval}
cardloc=${fullunionf}.card
timeout=${outprefix}_fullunionf_k${kval}.out

if [ $approach == 'kmc' ]; then
  #echo "kmc -k${kval} -ci1 -fm \@${apout}/all_fastas.txt ${fullunionf} ${outdir}/kmcfullf 2> ${fullunionf}.card"
  /usr/bin/time -o ${timeout} -v sh -c "kmc  -k${kval} -ci1 -fm @${apout}/all_fastas.txt ${fullunionf} ${outdir}/kmc > ${cardloc}"

elif [ $approach == 'dashing' ]; then
cmd="${dashing} hll -k ${kval}  -S 20 ${apout}/all_fastas.txt  > ${cardloc}"
echo ${cmd}
  /usr/bin/time -v -o ${timeout} sh -c "$cmd"
  card=$(awk 'NR==2{print $NF}' ${cardloc})
  echo cardloc $cardloc
fi
mrss=$(grep "Maximum resident" ${timeout} | awk '{print $NF}')
wctime=$(grep "wall clock" ${timeout} | awk '{print $NF}')
systime=$(grep "System time" ${timeout} | awk '{print $NF}')
utime=$(grep "User time" ${timeout} | awk '{print $NF}')
echo "${approach},${stage},all,${fullunionf},NA,${nfasta},${card},${mrss},${wctime},${systime},${utime}" >> ${summary}

#/usr/bin/time -o ${outprefix}_c${howmany}.out -v kmc -ci1 -k${kval} -fm @${apout}/c${howmany}.txt ${apout}/c${howmany}_k${kval}.kmc ${outdir}/kmc > ${apout}/c${howmany}_k${kval}.card
