#!/bin/bash


function usage()
{
        cat<<EOF
Usage:

merged_somatic_filter_vcf.sh -i <invcf.gz> -c <chromosome> -g <minGQ> -d <minDP> -f

Options:

        -i    input vcf [mandatory]
        -c    chromosome [INT; mandatory]
        -g    define minGQ [default=50]
        -d    define minDP [default=20]
        -f    retaining also variants with LOW AND HIGH support [optional; default=filter only HIGH SUPPORT]
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

while getopts ":i:o:c:g:d:fh" opt; do
        case $opt in

                i)
                        input=$OPTARG
                        ;;
#                o)
#                        output=$OPTARG
#                        ;;
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

#if [ "$output-" == "-" ] 
#then
#	echo "ERROR: missing output file"
#	exit 1
#fi

if [ "$chr-" == "-" ]
then
        echo "ERROR: missing chromosome"
        exit 1
fi

set -e
if [ $usePass == true ]
then
	echo 'keeping only HIGH SUPPORT variants according to the variant calling pipeline filters' > /dev/stderr
	filter_somatic_vcf.py $input 0.0 $chr y
	echo 'Done' > /dev/stderr

else
	echo 'retaining variants with LOW and HIGH SUPPORT FILTER field' > /dev/stderr
	filter_somatic_vcf.py $input 0.0 $chr n
        echo 'Done' > /dev/stderr

fi


exit 0


