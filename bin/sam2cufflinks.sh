
#
# To add the custom field of the strand of transcription to sam files for cufflink
# 
# USAGE:
# samtools view -S alignment.sam | ./sam2cufflinks.sh <sense_mate>
# 
# Choose <sense_mate> between:
# - "MATE1_SENSE" 
# - "MATE2_SENSE"
#
# It writes on stdout
# 


# comments on the flag:
# 0x10 --> reverse complemented 
# 0x40 --> the read is MATE1
# 0x50 --> the read is MATE2

mate=$1


if [[ $mate == "MATE1_SENSE" ]]; then
	awk 'BEGIN{FS="\t";OFS="\t"}{ 
		if (/^@/) {print $0; next}
		if (and($2,0x10) && and($2,0x40)) 
		{cufftag="XS:A:-"}
		if (and($2,0x10) && and($2,0x80))
		{cufftag="XS:A:+"}
		if (!and($2,0x10) && and($2,0x40))
		{cufftag="XS:A:+"}
		if (!and($2,0x10) && and($2,0x80))
		{cufftag="XS:A:-"};
		{print $0, cufftag}
	}'
fi

if [[ $mate == "MATE2_SENSE" ]]; then
	awk 'BEGIN{FS="\t";OFS="\t"}{ 
		if (/^@/) {print $0; next}
		if (and($2,0x10) && and($2,0x40)) 
		{cufftag="XS:A:+"}
		if (and($2,0x10) && and($2,0x80))
		{cufftag="XS:A:-"}
		if (!and($2,0x10) && and($2,0x40))
		{cufftag="XS:A:-"}
		if (!and($2,0x10) && and($2,0x80))
		{cufftag="XS:A:+"}
		{print $0, cufftag}
	}'
fi

if [[ $mate == 'NONE' ]]; then
	cat -
fi
