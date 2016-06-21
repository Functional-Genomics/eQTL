#!/bin/bash


function usage()
{
        cat<<EOF
Usage:

filter_consensus_vcf.sh -v min_vaf  -i <invcf.gz> -o <outvcf.gz> -c <chromosome> -g <minGQ> -d <minDP>  

Options:

        -i    input vcf [mandatory]
	-o    output vcf [mandatory]
        -c    chromosome [INT; mandatory]
        -g    define minGQ [default=50]
        -d    define minDP [default=20]
        -f    retaining also variants without PASS in FORMAT field [default=true]
	-v    minimum vaf
	-h    show this help message

EOF

}

if [ $# -eq 0 ]
then
        usage
        exit 1
fi


usePass=true
#useProb=false
minGQ=50
minDP=20
minvaf=0.1

while getopts ":v:i:o:c:g:d:fh" opt; do
        case $opt in
                v)
                        minvaf=$OPTARG
                        ;;
                i)
                        input=$OPTARG
                        ;;
                o)
                        output=$OPTARG
                        ;;
                c)
                        chr=$OPTARG
                        ;;
                g)
                        minGQ=$OPTARG
                        ;;
                d)
                        minDP=$OPTARG
                        ;;
                f)
                        usePass=false
                        ;;
                h)
                        usage
                        exit 1
                        ;;
                *)
                        echo "Invalid option:"
                        exit 1
                        ;;
         esac

done

if [ "$input-" == "-" ]
then
	echo "ERROR: missing input file"
	usage
	exit 1
fi

if [ ! -e "$input" ] 
then
	echo "ERROR: $input file not found"
	exit 1
fi

if [ "$output-" == "-" ] 
then
	echo "ERROR: missing output file"
	exit 1
fi

if [ "$chr-" == "-" ]
then
        echo "ERROR: missing chromosome"
        exit 1
fi

filter_somatic_vcf.py $input $minvaf $chr | bgzip -c > $output && \
tabix -p vcf $output

exit 0


