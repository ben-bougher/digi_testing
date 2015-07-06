# plotCMP.sh
# Convert binary array of CMP gather to .su file and generate wiggle plot
# Input: CMP gather, number of samples per trace (ns), output file name
# Output: PS image of CMP wiggle plot
infile=$1
numSamples=$2
outfile=$3
#--------------------------------------------------------------------------
# Check if user supplied correct number of arguments
if [ $# -ne 3 ]; then
	echo " "
	echo "USAGE $0 [infile ns outfile]"
	echo " "
	
	exit 1
fi

suaddhead < $infile ftn=0 ns=$numSamples > tmp.su

supswigp d1=1.0 < tmp.su > $outfile.ps


