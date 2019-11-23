#! /bin/bash
function execute_cases() {
	Casename=$1
	# read -p 'ddd' var
	# np series
	for np in ${npseq[@]}
	do
		for logRA in ${logRAseq[@]}
		do
			for npt in ${nptseq[@]}
			do
				for ind in $(seq 0 $[${#Periodseq[@]}-1])
				do
					P=${Periodseq[ind]}
					DT=${DTseq[ind]}
					# modify setups in INPUT.txt
					sed "s/RA=1E5/RA=1E${logRA}/g" $inputfile_backup > $inputfile
					sed -i "s/DT=0.1/DT=${DT}/g"                       $inputfile
					sed -i "s/NT=1000/NT=${NT}/g"                      $inputfile
					sed -i "s/Period=10/Period=${P}/g"                $inputfile
					sed -i "s/Amplitude=0/Amplitude=${Amplitude}/g"  $inputfile
					sed -i "s/Npt=16/Npt=${npt}/g"  $inputfile

					title="${Casename}_RA=1E${logRA},npt=${npt},np=${np},P=${P}"
					TMPDIR=$outputdir/$title
					if [ ! -d $TMPDIR ]; then
						mkdir $TMPDIR && mkdir $TMPDIR/OUTPUT
						cp -r $workdir/INPUT $TMPDIR
					fi

					# run parallel program
					cmd='mpiexec -n $np -wdir $TMPDIR $exe'	
					eval $cmd >$TMPDIR/logfile 2>&1
					echo "CASE${NCases} has been resolved"
					((NCases+=1))					
				done
			done
		done
	done
	#compname="${Casename}.tar.gz"
	#tar -czf ${outputdir}/${compname} ${outputdir}/${Casename}* --remove-files
}




workdir=$(pwd)
exe=$workdir/build/bin/PIDTS
inputfile_backup=$workdir/INPUT/INPUT_BACKUP.txt
inputfile=$workdir/INPUT/INPUT.txt
outputdir=$workdir/OUTPUT

NCases=0

# test cases
logRAseq=(3 4 5 6)
npseq=(1 4 8 12 16 20)
nptseq=(1 4 16 32 64 128)
Amplitude=0.8
Periodseq=10
DTseq=0.1
NT=1000

# run
execute_cases TestCase

echo "SUM: ${NCases} cases"



















