#/usr/bin/bash

# script to do a quick liftOver between 2 genomes using BLAT
# uses GridEngine dispatching (limits number on concurrently queued jobs)
# query (NEW) lifted to
# target (OLD) ifted from

usage(){
	echo ""
	echo "#######################################"
	echo "### EASY CREATION OF LIFTOVER FILES ###"
	echo "#######################################"
	echo ""
	echo "Usage: $0 new_assembly old_assembly [blat,liftup,chain,net/smallnet]"
	echo ""
	exit 1
}
[[ $# -eq 0 ]] && usage

# CONFIGURATION
QUERY=${1}
TARGET=${2}
BLAT_DONE=${3}

# DIRECTORIES
BIN=./bin/
WORKDIR=./tmp/
CHAINDIR=./chain/
NETDIR=./net/
LOGDIR=./logs/
OVERDIR=./over/
# OTHER PARAMETERS
CHUNKSIZE=200000
MERS=11
# GRIDENGINE (and throttling)
QUEUE=vert
#QUEUE=medium_jobs.q
QSUB="qsub -cwd -b y -V -q ${QUEUE} -o /dev/null "
STOPSUBMIT=60
MANYSUBMIT=10
# BLAT PARAMETERS
MIN_IDENT=98
MIN_SCORE=100
# progress files (hidden)
STEP_SPLIT=./done.split
STEP_TWOBIT=./done.twobit
STEP_OOC=./done.ooc

# create directories
if [[ ! -d ${WORKDIR} ]];
	then
	mkdir ${WORKDIR}
fi
if [[ ! -d ${OVERDIR} ]];
	then
	mkdir ${OVERDIR}
fi
if [[ ! -d ${NETDIR} ]];
	then
	mkdir ${NETDIR}
fi
if [[ ! -d ${CHAINDIR} ]];
	then
	mkdir ${CHAINDIR}
fi
if [[ ! -d ${LOGDIR} ]];
	then
	mkdir ${LOGDIR}
fi

#### POST_QUEUE
if [ $3 ];
then
	if [ $BLAT_DONE == 'blat' ];
	then
		echo "PREPARING BLAT JOBS"
		CHRFILES=$( find $WORKDIR -maxdepth 1 -type f -name 'chr*_*.fa' | sort )
		COUNTER=0
		for QUERY_CHR in $CHRFILES
		do
			OUTFILE=${WORKDIR}`echo ${QUERY_CHR} | sed 's/.*\/\(.*\)$/\1\.psl/'`
			JOB="${BIN}blat ${TARGET}.2bit ${QUERY_CHR} ${OUTFILE} -ooc=${WORKDIR}${MERS}.ooc -tileSize=${MERS} -minScore=${MIN_SCORE} -minIdentity=${MIN_IDENT} -fastMap"
			JOBNAME=`echo ${QUERY_CHR} | perl -pe 's/\D+(\d+)\D+(\d+)/BL\1\_\2/'`
			if [[ ! -e ${OUTFILE} ]];
				then
				$QSUB -N $JOBNAME -o ${LOGDIR}${JOBNAME}.stdout -e ${LOGDIR}${JOBNAME}.stderr ${JOB}
				((++COUNTER))
			else
				echo "${OUTFILE} exists -> skipped"
			fi
			if [[ $COUNTER -ge $MANYSUBMIT ]];
				then
				echo '-> QUEUE SUBMISSION THROTTLING'
				sleep 1
				SIZESUBMIT=`qstat | wc -l`
				echo "${SIZESUBMIT} running (${STOPSUBMIT} limit)"
				while [[ ${SIZESUBMIT} -ge ${STOPSUBMIT} ]];
					do
					sleep 60
					SIZESUBMIT=`qstat | wc -l`
					echo "$SIZESUBMIT jobs waiting in queue"
				done
				COUNTER=0
			fi
		done
		echo "execute ${SUBMIT} wait for completion then rerun with:"
		echo "     $0 $1 $2 liftup"
		exit 0
	fi

	if [ $BLAT_DONE == 'liftup' ];
	then
		echo "Lifting up chunks"
		CHRFILES=$( find $WORKDIR -maxdepth 1 -type f -name 'chr*.lft' | sort | perl -pe 's/\.lft//g' )
		for CHR in $CHRFILES
		do
			echo "  $CHR"
			if [[ ! -e ${CHR}.psl ]]; then
				${BIN}liftUp -pslQ ${CHR}.psl ${CHR}.lft warn ${CHR}_*.psl && rm -f ${CHR}_*.psl
			fi
		done
	    echo "Rerun with:"
		echo "     $0 $1 $2 chain"
		exit 0
	fi

	if [ $BLAT_DONE == 'chain' ];
	then
		PSLFILES=$( find $WORKDIR -maxdepth 1 -type f -follow -name 'chr*.psl' | grep -v 'fa' | sort -r )
		for PSL in $PSLFILES
	    do
	    	echo ${PSL}
	    	OUTFILE=${CHAINDIR}`echo ${PSL} | sed 's/.*\/\(.*\)\.psl$/\1\.chain/'`
			JOB="${BIN}axtChain -linearGap=medium -psl ${PSL} ${TARGET}.2bit ${QUERY}.2bit ${OUTFILE}"
   			JOBNAME=`echo ${PSL} | perl -pe 's/\D+(\d+)\D+/CH\1/'`
			if [[ ! -e ${OUTFILE} ]];
				then
				#echo "$QSUB -N $JOBNAME -o ${LOGDIR}${JOBNAME}.stdout -e ${LOGDIR}${JOBNAME}.stderr ${JOB}"
				$QSUB -N $JOBNAME -o ${LOGDIR}${JOBNAME}.stdout -e ${LOGDIR}${JOBNAME}.stderr ${JOB}
				((++COUNTER))
			else
				echo "${OUTFILE} exists -> skipped"
			fi
			if [[ $COUNTER -ge $MANYSUBMIT ]];
				then
				sleep 1
				SIZESUBMIT=`qstat | wc -l`
				echo "${SIZESUBMIT} running (${STOPSUBMIT} limit)"
				while [[ ${SIZESUBMIT} -ge ${STOPSUBMIT} ]];
					do
					sleep 60
					SIZESUBMIT=`qstat | wc -l`
					echo "$SIZESUBMIT jobs waiting in queue"
				done
				COUNTER=0
			fi
	    done
	    echo "execute ${SUBMIT} wait for completion then rerun with:"
		echo "     $0 $1 $2 net"
		exit 0
	fi

	if [ $BLAT_DONE == 'net' ];
	then
		echo "Sorting chains"
		${BIN}chainMergeSort ${CHAINDIR}*.chain | ${BIN}chainSplit chain stdin
		tar cvfz ${CHAINDIR}unsorted_chains.tgz ${CHAINDIR}/chr*.chain && rm ${CHAINDIR}/chr*.chain

		echo "Netting"
		CHAINFILES=$( find $CHAINDIR -maxdepth 1 -type f -follow -name '*.chain' | sort )
		for i in $CHAINFILES
		do
			echo "netting ${i} ..."
			NETFILE=${NETDIR}`echo ${i} | sed 's/.*\/\(.*\).chain$/\1\.net/'`
			OVERFILE=${OVERDIR}`echo ${i} | sed 's/.*\/\(.*\).chain$/\1\.chain/'`
			${BIN}chainNet ${i} ${TARGET}.sizes ${QUERY}.sizes ${NETFILE} /dev/null
			${BIN}netChainSubset ${NETFILE} $i ${OVERFILE}
		done
		echo "Finishing"
		cat ${OVERDIR}*.chain > liftOver.chain
		echo "DONE"
		exit 0
	fi

	if [ $BLAT_DONE == 'smallnet' ];
	then
		echo "Sorting chains"
		${BIN}chainMergeSort ${CHAINDIR}*.chain | ${BIN}chainSplit chainMerge stdin -lump 50
		tar cvfz ${CHAINDIR}unsorted_chains.tgz ${CHAINDIR}/chr*.chain && rm ${CHAINDIR}/chr*.chain
		
		cat chainMerge/*.chain > all.chain
		${BIN}chainSort all.chain all.sorted.chain

		echo "Netting and finishing"
		${BIN}chainNet all.sorted.chain ${TARGET}.sizes ${QUERY}.sizes all.net /dev/null
		${BIN}netChainSubset all.net all.chain all.liftOver

		echo "DONE"
		exit 0
	fi
fi

### PREPARE
echo "SPLITTING QUERY"
if [ ! -e ${STEP_SPLIT} ];
then
	QUERYSEQS=`grep -c "^>" ${QUERY}`
	${BIN}faSplit sequence ${QUERY} ${QUERYSEQS} ${WORKDIR}chr
	CHRFILES=$( find $WORKDIR -maxdepth 1 -type f -follow -name 'chr*.fa' | sort | perl -pe 's/\.fa//g' )
	for CHR in $CHRFILES
	do
		echo "  ${CHR}"
		LIFT=${CHR}.lft
		${BIN}faSplit -lift=${LIFT} size ${CHR}.fa ${CHUNKSIZE} ${CHR}_
	done
	touch ${STEP_SPLIT}
else
	echo "...done already"
fi
echo "CREATING 2bit INDIZES AND GETTING SIZES"
if [ ! -e ${STEP_TWOBIT} ];
then
	${BIN}faToTwoBit ${TARGET} ${TARGET}.2bit
	${BIN}faToTwoBit ${QUERY} ${QUERY}.2bit
	${BIN}twoBitInfo ${TARGET}.2bit ${TARGET}.sizes
	${BIN}twoBitInfo ${QUERY}.2bit ${QUERY}.sizes
	touch ${STEP_TWOBIT}
else
	echo "...done already"
fi
echo "CREATING ooc FILE (-mers)"
if [ ! -e ${STEP_OOC} ];
then
	${BIN}blat ${TARGET}.2bit /dev/null /dev/null -tileSize=${MERS} -makeOoc=${WORKDIR}${MERS}.ooc -repMatch=1024
	touch ${STEP_OOC}
else
	echo "...done already"
fi
echo "Everthing ready. Run those =>  blat,liftup,chain,net"


