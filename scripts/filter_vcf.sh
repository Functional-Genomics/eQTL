#!/bin/bash


function usage()
{
        cat<<EOF
Usage:

filter_vcf.sh -i <invcf.gz> -o <outvcf.gz> -c <chromosome> -g <minGQ> -d <minDP>  

Options:

        -i    input vcf [mandatory]
	-o    output vcf [mandatory]
        -c    chromosome [INT; mandatory]
        -g    define minGQ [default=50]
        -d    define minDP [default=20]
        -f    retaining also variants without PASS in FORMAT field [default=true]
	#-v   specify variant calling type [Sanger | Annai | ... ]
        #-p   select only those variants with >=0.95 probability of the most probable genotype [default=false]
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


if [ $usePass == true ]
then
	echo 'keeping only PASS variants according to the variant calling pipeline filters'

	vcftools --gzvcf $input --minGQ $minGQ --minDP $minDP --recode --recode-INFO-all --remove-filtered-all --chr $chr --stdout | bgzip -c > $output &&

	tabix -p vcf $output

	echo 'Done'

else
	echo 'retaining variants without PASS in FORMAT field'

        vcftools --gzvcf $input --minGQ $minGQ --minDP $minDP --recode --recode-INFO-all --chr $chr --stdout | bgzip -c > $output &&

        tabix -p vcf $output

        echo 'Done'

fi



#line I was using for Sanger pipeline vcf files

#name=$(echo $input | sed 's/\.svcp.*gz//g'); vcftools --gzvcf $input --minGQ $minGQ --minDP $minDP --recode --recode-INFO-all --remove-filtered-all --chr $chr --stdout | cut -f 1,2,3,4,5,6,7,8,9,11 | grep -v '##SAMPLE=<ID=NORMAL' | sed 's/TUMOUR/'$name'/g' | vcffilter -f " TP > 0.95 | TP = 0.95 | SP > 0.95 | SP = 0.95 " | bgzip -c > ${name}_chr${chr}_somatic_snv_mn.vcf.gz


exit 0


