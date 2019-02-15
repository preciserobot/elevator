# creates liftOver file from 2 genomes (FASTA) and a pairwise genome alignment (PSL)

QUERY=$1
TARGET=$2
PSL=$3
BIN=bin

# create 2bit index and size
if [[ ! -e ${QUERY}.sizes ]]; then
	${BIN}/faToTwoBit ${QUERY} ${QUERY}.2bit
	${BIN}/twoBitInfo ${QUERY}.2bit ${QUERY}.sizes
fi
if [[ ! -e ${TARGET}.sizes ]]; then
	${BIN}/faToTwoBit ${TARGET} ${TARGET}.2bit
	${BIN}/twoBitInfo ${TARGET}.2bit ${TARGET}.sizes
fi

# chaining
echo "CHAINING $PSL ..."
CHAINDIR=chain_`basename $PSL .psl`
mkdir ${CHAINDIR}
${BIN}/axtChain -linearGap=medium -psl $PSL ${TARGET}.2bit ${QUERY}.2bit ${CHAINDIR}/`basename $PSL .psl`.chain

# merge chains
echo "MERGING ${PSL} ..."
CHAINMERGEDIR=chainMerge_`basename $PSL .psl`
mkdir ${CHAINMERGEDIR}
${BIN}/chainMergeSort ${CHAINDIR}/*.chain | ${BIN}/chainSplit ${CHAINMERGEDIR} stdin -lump=50
ALLCHAINS=`basename $PSL .psl`.allchains
cat ${CHAINMERGEDIR}/*.chain > $ALLCHAINS
${BIN}/chainSort ${ALLCHAINS} ${ALLCHAINS}.sorted

# netting
echo "NETTING ${PSL} ..."
NETDIR=net_`basename $PSL .psl`
mkdir ${NETDIR}
${BIN}/chainNet ${ALLCHAINS}.sorted ${TARGET}.sizes ${QUERY}.sizes ${NETDIR}/all.net /dev/null
${BIN}/netChainSubset ${NETDIR}/all.net ${ALLCHAINS}.sorted `basename $PSL .psl`.liftOver

# cleanup (silent)
if [[ -e `basename $PSL .psl`.liftOver ]]; then
	rm -rf ${NETDIR} ${CHAINMERGEDIR} ${CHAINDIR}
fi
