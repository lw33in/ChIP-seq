#!/usr/bin/perl -w

# An input file containing a list of sample fastqs
# This is Encode3 pipeline


use strict;
use File::Basename;
use POSIX qw(strftime);
use Cwd qw(abs_path);

my $usage = "\nUsage: ENCODE3Pipeline.pl <sampleFastqList.txt>\n\nOptions:\n\n \n\nExamples:\nENCODE3Pipeline.pl /export/PipelineRuns/Encode3-ChIP-Seq/test/fastqList.txt\n\n  The input FastqList should have the following format:\n\nUSC_ID  SEP_ID  DISEASE  MARK  PEAK_TYPE(SPP use \"NA\"  MACS2 use \"narrow\", \"broad\", or \"gapped\") PEAK_SW(\"MACS2\" or \"SPP\") CONTROL_TAG_ALIGN_FILE  CASE_FASTQ_LOCATION (paired end fastqs should be separated by a comma) \nExample input line:\nUSC7788    SEP044    H293   H3K4me3      narrow  MACS2   /export/ChIPseq/ENCODE3.MACS2.from.fastq/Input.filt.nodup.srt.SE.tagAlign.gz      /export/ChIP-Seq/CombinedFastqs/H3K4me3.fastq.gz\n\n"; 

my $abspath = abs_path();
my $tstamp = strftime("%Y-%m-%d-%H-%M",localtime);
my $input = shift || die $usage;

# BWA INDEX PREFIX
my $BWA_INDEX_NAME = "/home/ENCODE3.pipeline.test/pipeline.important.files/encode_hg19_mf/male.hg19.fa";

# SOFTWARE PACKAGES
# Current default 
my $SAMTOOLS = "/home/ENCODE3.pipeline.test/pipeline.important.files/samtools-0.1.19";

my $PICARD = "/home/ENCODE3.pipeline.test/pipeline.important.files/picard-tools-1.46";
## Current default is picard-tools-1.90
#my $PICARD = "/export/uec-gs1/knowles/uec-01/knowles/software/picard/default";

# use bedtools 2.19.1
# Current default is BEDTools-Version-2.2.19
my $BEDTOOLS = "/usr/bin";

# Current default is UCSC_Browser_Tools
my $UCSC_BROWSER_TOOLS = "/home/ENCODE3.pipeline.test/pipeline.important.files/UCSC_Browser_Tools/";

# Current default is bwa 0.7.12
my $BWA = "/home/ENCODE3.pipeline.test/pipeline.important.files/bwa-0.7.12";

my $MACS2 = "/home/rcf-proj/swn/swnVE/bin";

# Current igvtools 
my $IGVTOOLS = "/home/ENCODE3.pipeline.test/bioinformaticstools/IGVTools";

#CHRSIZEFILE
my $CHRSIZEFILE="/home/ENCODE3.pipeline.test/bioinformaticstools/male.hg19.genome";

my $OFPREFIX;
my $FINAL_BAM_PREFIX;
my $CHIP_TA_PREFIX;
my $PEAK_OUTPUT_DIR;
my $SUBSAMPLED_TA_FILE;
my $USC_ID;
my $SEP_ID;
my $disease;
my $mark;
my $peakType;
my $peakSoftware;
my $controlFile;
my $fastq;
my $fastq_1;
my $fastq_2;


my @inputFileList;

open (INPUT_LIST, "<$input") or die "Couldn't open $input for reading\n";
# Full path of input file. 
while (<INPUT_LIST>) {
   # skip blank lines and comment lines beginning with hash
   next if m/^(#|\s*$)/;
   chomp;
   my $numFields = scalar @{[split '\s+', $_]};
   die "\nERROR!:  Missing column on line:\n$_\n\nInput line should contain the following 8 columns separated by whitespace:\nUSC_ID  SEP_ID  DISEASE  MARK  PEAK_TYPE(SPP use NA  MACS2 use narrow, broad, or gapped) PEAK_SW(macs2 or spp) CONTROL_TAG_ALIGN_FILE             CASE_FASTQ_LOCATION (paired end fastqs should be separated by a comma)\n\n" if ($numFields != 8);
   push (@inputFileList, $_);
}
close INPUT_LIST; 

chomp (@inputFileList);
die "\nERROR: No input files found\n" if !@inputFileList;

foreach (@inputFileList) {
    
    my @fields = split '\s+', $_;
    $USC_ID = $fields[0];
    $SEP_ID = $fields[1];
    $disease = $fields[2];
    $mark = $fields[3];
    $peakType = $fields[4];
    $peakSoftware = $fields[5];
    $controlFile = $fields[6];
    $fastq = $fields[7];
    $fastq_1 = $fastq;

    # Check for commas designating paired end
    # If there is a comma in $fastq assume it's paired end 
    if ($fastq =~ /\,/) {
        @fields = split ',', $fastq;
        $fastq_1 = $fields[0];
        $fastq_2 = $fields[1];
print "\nPAIRED END:\nR1 = $fastq_1\nR2 = $fastq_2\nUSC_ID = $USC_ID\nSEP_ID = $SEP_ID\ndisease = $disease\nmark = $mark\npeakType = $peakType\npeakSoftware = $peakSoftware\n";
    } else {
print "\nSINGLE END:\n$fastq\nUSC_ID = $USC_ID\nSEP_ID = $SEP_ID\ndisease = $disease\nmark = $mark\npeakType = $peakType\npeakSoftware = $peakSoftware\n";
    }

    # VARS FOR PBS SCRIPTS
    $OFPREFIX="${SEP_ID}_${USC_ID}_${mark}_${disease}";
    $FINAL_BAM_PREFIX="${OFPREFIX}.filt.nodup.srt";

     # Create sample dir
    `mkdir $OFPREFIX`;

     my $bwaJob;
     my $libraryComplexityJob;
     my $tagAlignJob;
     my $peaksJob;
     my $tracksJob;
     my $statsJob;

     # FOR SINGLE END
     if (!$fastq_2) {
         printf "Creating slurm scripts for Single End processing for $fastq\n\n";

         # STEP 1 BWA:
         createSLURM_BWA_SE("$fastq");
         $bwaJob = "$OFPREFIX-Encode3Pipeline-bwa-singleEnd.slurm";

         # STEP 2A LIBRARY COMPLEXITY:
         createSLURM_LibraryComplexity_SE();
         $libraryComplexityJob = "$OFPREFIX-Encode3Pipeline-computeLibraryComplexity-singleEnd.slurm";

         # STEP 2B TAG ALIGN:
         createSLURM_TagAlign_SE();
         $tagAlignJob = "$OFPREFIX-Encode3Pipeline-tagAlign-singleEnd.slurm";

         # STEP 3 MACS2 OR SPP
         if ($peakSoftware =~ /MACS2/i) {
             $PEAK_OUTPUT_DIR="MACS2_PEAKS";
             $CHIP_TA_PREFIX="${FINAL_BAM_PREFIX}.MACS2";
             createSLURM_MACS2("SE");
             $peaksJob = "$OFPREFIX-Encode3Pipeline-MACS2.slurm";
         } elsif ($peakSoftware =~ /SPP/i) { 
             $CHIP_TA_PREFIX="${FINAL_BAM_PREFIX}.SPP";
             $PEAK_OUTPUT_DIR="SPP_PEAKS";
             createSLURM_SPP("SE");
             $peaksJob = "$OFPREFIX-Encode3Pipeline-SPP.slurm";
         } else { 
             warn "WARNING: NO PEAK SOFTWARE FOUND! PLEASE CHECK CONFIG FILE!\n";
         }

         # STEP 4 TRACKS
         createSLURM_tracks("SE");
         $tracksJob = "$OFPREFIX-Encode3Pipeline-tracks.slurm";

         # STEP 5 Get pipeline stats
         createSLURM_GetStats();
         $statsJob = "$OFPREFIX-Encode3Pipeline-GetStats.slurm";

     } # END SINGLE END

     # FOR PAIRED END
     else {  #if ($fastq_2)
        printf "Creating slurm scripts for Paired End processing\n\n";

         # STEP 1 BWA:
         createSLURM_BWA_PE($fastq_1, $fastq_2);
         $bwaJob = "$OFPREFIX-Encode3Pipeline-bwa-pairedEnd.slurm";

         # STEP 2A LIBRARY COMPLEXITY:
         createSLURM_LibraryComplexity_PE();
         $libraryComplexityJob = "$OFPREFIX-Encode3Pipeline-computeLibraryComplexity-pairedEnd.slurm";
     
         # STEP 2B TAG ALIGN:
         createSLURM_TagAlign_PE();
         $tagAlignJob = "$OFPREFIX-Encode3Pipeline-tagAlign-pairedEnd.slurm";

         # STEP 3 MACS2 OR SPP
         if ($peakSoftware =~ /MACS2/i) {
                $PEAK_OUTPUT_DIR="MACS2_PEAKS";
                $CHIP_TA_PREFIX="${FINAL_BAM_PREFIX}.MACS2";
                createSLURM_MACS2("PE2SE");
                $peaksJob = "$OFPREFIX-Encode3Pipeline-MACS2.slurm";
         } elsif ($peakSoftware =~ /SPP/i) {
                $CHIP_TA_PREFIX="${FINAL_BAM_PREFIX}.SPP";
                $PEAK_OUTPUT_DIR="SPP_PEAKS";
                createSLURM_SPP("PE2SE");
                $peaksJob = "$OFPREFIX-Encode3Pipeline-SPP.slurm";
         } else { 
                warn "WARNING: NO PEAK SOFTWARE FOUND! PLEASE CHECK CONFIG FILE!\n";
         }
         
         # STEP 4 TRACKS
         createSLURM_tracks("PE2SE");
         $tracksJob = "$OFPREFIX-Encode3Pipeline-tracks.slurm";

         # STEP 5 Get pipeline stats
         createSLURM_GetStats();
         $statsJob = "$OFPREFIX-Encode3Pipeline-GetStats.slurm";

     }  # END PAIRED END

    my $bwaJobID;
    my $bwaJobIDNumber;
    if ($bwaJob) {
         $bwaJobID = `sbatch $bwaJob`;
	 $bwaJobIDNumber = @{[$bwaJobID =~ m/\w+/g]}[3];
         chomp ($bwaJobID);
         if($bwaJobID) {
             print "Submitting $bwaJob: $bwaJobIDNumber\n";
         } else {
             print "Problem submitting $bwaJob\n";
         }
    }

    if ($libraryComplexityJob) {
        my $libraryComplexityJobID = `sbatch --dependency=afterok:$bwaJobIDNumber $libraryComplexityJob`;
        my $libraryComplexityJobIDNumber = @{[$libraryComplexityJobID =~ m/\w+/g]}[3];
        chomp ($libraryComplexityJobID);
        if($libraryComplexityJobID) {
            print "Submitting $libraryComplexityJob: $libraryComplexityJobIDNumber\n";
        } else {
            print "Problem submitting $libraryComplexityJob\n";
        }
    }

    my $tagAlignJobID;
    my $tagAlignJobIDNumber;
    if ($tagAlignJob) {
         $tagAlignJobID = `sbatch --dependency=afterok:$bwaJobIDNumber $tagAlignJob`;
         $tagAlignJobIDNumber = @{[$tagAlignJobID =~ m/\w+/g]}[3];
         chomp ($tagAlignJobID);
         if($tagAlignJobID) {
             print "Submitting $tagAlignJob: $tagAlignJobIDNumber\n";
         } else {
             print "Problem submitting $tagAlignJob\n";
         }
    }

    my $peaksJobID;
    my $peaksJobIDNumber;
    if ($peaksJob) {
         $peaksJobID = `sbatch --dependency=afterok:$tagAlignJobIDNumber $peaksJob`;
         $peaksJobIDNumber = @{[$peaksJobID =~ m/\w+/g]}[3];
	 chomp ($peaksJobID);
         if($peaksJobID) {
             print "Submitting $peaksJob: $peaksJobIDNumber\n";
         } else {
             print "Problem submitting $peaksJob\n";
         }
    }

    if ($tracksJob) {
         my $tracksJobID = `sbatch --dependency=afterok:$peaksJobIDNumber $tracksJob`;
	 chomp ($tracksJobID);
         if($tracksJobID) {
             print "Submitting $tracksJob: $tracksJobID\n";
         } else {
             print "Problem submitting $tracksJob\n";
         }
    }

    if ($statsJob) {
         my $statsJobID = `sbatch --dependency=afterok:$peaksJobIDNumber $statsJob`;
         chomp ($statsJobID);
         if($statsJobID) {
             print "Submitting $statsJob: $statsJobID\n";
         } else {
             print "Problem submitting $statsJob\n";
         }
    }

} # foreach (@inputFileList)
exit(0);

###
# Parts taken from /export/ChIPseq/ENCODE3/For_3reps/SE/FIN/TF_SE_using_MACS2_rep1.sh
###
sub createSLURM_BWA_SE{
    my @list = @_;
    my $fastq = $list[0];

    # Create PBS script to run BWA alignnment and run post-alignment filtering
    open(SLURMJOB, " > $OFPREFIX-Encode3Pipeline-bwa-singleEnd.slurm") || die "$!\n";
    print SLURMJOB <<PBS;
#!/bin/bash
#SBATCH --job-name=$OFPREFIX-Encode3Pipeline-bwa-singleEnd
#SBATCH --ntasks=16 --mem-per-cpu=16GB
#SBATCH --time=23:59:59
#SBATCH --output=$OFPREFIX-Encode3Pipeline-bwa-singleEnd.out
#SBATCH --error=$OFPREFIX-Encode3Pipeline-bwa-singleEnd.err
#
cd \$SLURM_SUBMIT_DIR
cd $OFPREFIX
source /home/rcf-proj/swn/swnVE/bin/activate
#source /usr/usc/python/2.7.8/setup.sh
#R_LIBS_USER=/home/bioinformatic_tools/R-3.2.0/R-3.2.0/library

#===========================
## 1a. Read alignment (BWA)
###
# Run BWA alignment
###

# works fine both with gz file and unzipped file
FASTQ_FILE_1=$fastq

NTHREADS=16

# map reads to create raw SAM file 
SAI_FILE_1="${OFPREFIX}.sai"
RAW_BAM_PREFIX="${OFPREFIX}.raw.srt"
RAW_BAM_FILE="\${RAW_BAM_PREFIX}.bam" 
RAW_BAM_FILE_MAPSTATS="\${RAW_BAM_PREFIX}.flagstat.qc" # QC File

$BWA/bwa aln -q 5 -l 32 -k 2 -t \${NTHREADS} ${BWA_INDEX_NAME} \${FASTQ_FILE_1} > \${SAI_FILE_1}
$BWA/bwa samse ${BWA_INDEX_NAME} \${SAI_FILE_1} \${FASTQ_FILE_1} | $SAMTOOLS/samtools view -Su - | samtools sort - \${RAW_BAM_PREFIX}
#rm \${SAI_FILE_1}
$SAMTOOLS/samtools flagstat \${RAW_BAM_FILE} > \${RAW_BAM_FILE_MAPSTATS}

#===============================
## 1b. Post-alignment filtering 
###
# Run Post-alignment filtering
###

FILT_BAM_PREFIX="${OFPREFIX}.filt.srt"
FILT_BAM_FILE="\${FILT_BAM_PREFIX}.bam"
MAPQ_THRESH=30

$SAMTOOLS/samtools view -F 1804 -q \${MAPQ_THRESH} -b \${RAW_BAM_FILE} > \${FILT_BAM_FILE}

# Mark duplicates 
TMP_FILT_BAM_FILE="\${FILT_BAM_PREFIX}.dupmark.bam"
MARKDUP="$PICARD/MarkDuplicates.jar"

DUP_FILE_QC="\${FILT_BAM_PREFIX}.dup.qc" # QC file
java -Xmx30G -jar \${MARKDUP} INPUT=\${FILT_BAM_FILE} OUTPUT=\${TMP_FILT_BAM_FILE} METRICS_FILE=\${DUP_FILE_QC} VALIDATION_STRINGENCY=LENIENT ASSUME_SORTED=true REMOVE_DUPLICATES=false
mv \${TMP_FILT_BAM_FILE} \${FILT_BAM_FILE}

# remove duplictes and index final position sorted BAM
FINAL_BAM_FILE="${FINAL_BAM_PREFIX}.bam" # To be stored 
FINAL_BAM_INDEX_FILE="${FINAL_BAM_PREFIX}.bai" # To be stored 
FINAL_BAM_FILE_MAPSTATS="${FINAL_BAM_PREFIX}.flagstat.qc" # QC file
$SAMTOOLS/samtools view -F 1804 -b \${FILT_BAM_FILE} > \${FINAL_BAM_FILE}

# Index Final BAM file
$SAMTOOLS/samtools index \${FINAL_BAM_FILE} \${FINAL_BAM_INDEX_FILE}
$SAMTOOLS/samtools flagstat \${FINAL_BAM_FILE} > \${FINAL_BAM_FILE_MAPSTATS}
PBS

close(SLURMJOB) || die "$!\n";
} #end sub createSLURM_BWA_SE

###
sub createSLURM_LibraryComplexity_SE{

    # Create PBS script to compute library complexity
    open(SLURMJOB, " > $OFPREFIX-Encode3Pipeline-computeLibraryComplexity-singleEnd.slurm") || die "$!\n";
    print SLURMJOB <<PBS;
#!/bin/bash
#SBATCH --job-name=$OFPREFIX-Encode3Pipeline-computeLibraryComplexity-singleEnd
#SBATCH --ntasks=16 --mem-per-cpu=16GB
#SBATCH --time=23:59:59
#SBATCH --output=$OFPREFIX-Encode3Pipeline-computeLibraryComplexity-singleEnd.out
#SBATCH --error=$OFPREFIX-Encode3Pipeline-computeLibraryComplexity-singleEnd.err
#
cd \$SLURM_SUBMIT_DIR
cd $OFPREFIX
source /home/rcf-proj/swn/swnVE/bin/activate
#source /usr/usc/python/2.7.8/setup.sh
#R_LIBS_USER=/home/bioinformatic_tools/R-3.2.0/R-3.2.0/library

# Compute library complexity
# sort by position and strand

FILT_BAM_PREFIX=${OFPREFIX}.filt.srt
FILT_BAM_FILE=\${FILT_BAM_PREFIX}.bam

# Obtain unique count statistics
# use bedtools 2.19.1
PBC_FILE_QC="${FINAL_BAM_PREFIX}.pbc.qc"
$BEDTOOLS/bedtools bamtobed -i \${FILT_BAM_FILE} | awk 'BEGIN{OFS="\\t"}{print \$1,\$2,\$3,\$6}' | grep -v 'chrM' | sort | uniq -c | awk 'BEGIN{mt=0;m0=0;m1=0;m2=0} (\$1==1){m1=m1+1} (\$1==2){m2=m2+1} {m0=m0+1} {mt=mt+\$1} END{printf "%d\\t%d\\t%d\\t%d\\t%f\\t%f\\t%f\\n",mt,m0,m1,m2,m0/mt,m1/m0,m1/m2}' > \${PBC_FILE_QC}

#rm \${FILT_BAM_FILE}
# PBC File output
# TotalReadPairs [tab] DistinctReadPairs [tab] OneReadPair [tab] TwoReadPairs [tab] NRF=Distinct/Total [tab] PBC1=OnePair/Distinct [tab] PBC2=OnePair/TwoPair
PBS

close(SLURMJOB) || die "$!\n";

} #end createSLURM_LibraryComplexity_SE


sub createSLURM_TagAlign_SE{

my $NREADS=15000000;
my $f = $NREADS/1000000;
$SUBSAMPLED_TA_FILE="${OFPREFIX}.filt.nodup.sample.$f.SE.tagAlign.gz";

     # Create PBS script to convert BAM to tagAlign (BED 3+3 format)
     open(SLURMJOB, " > $OFPREFIX-Encode3Pipeline-tagAlign-singleEnd.slurm") || die "$!\n";
     print SLURMJOB <<PBS;
#!/bin/bash
#SBATCH --job-name=$OFPREFIX-Encode3Pipeline-tagAlign-singleEnd
#SBATCH --ntasks=16 --mem-per-cpu=16GB
#SBATCH --time=23:59:59
#SBATCH --output=$OFPREFIX-Encode3Pipeline-tagAlign-singleEnd.out
#SBATCH --error=$OFPREFIX-Encode3Pipeline-tagAlign-singleEnd.err
#
cd \$SLURM_SUBMIT_DIR
cd $OFPREFIX
source /home/rcf-proj/swn/swnVE/bin/activate
#source /usr/usc/python/2.7.8/setup.sh
#R_LIBS_USER=/home/bioinformatic_tools/R-3.2.0/R-3.2.0/library

#===========================
### 2a. convert BAM to tagAlign (BED 3+3 format)
##

# bedtools 2.19.1, gawk shuf
# module add bedtools/2.19.1

# Create tagAlign file 
FINAL_BAM_FILE="${FINAL_BAM_PREFIX}.bam"
FINAL_TA_FILE="${FINAL_BAM_PREFIX}.SE.tagAlign.gz"

$BEDTOOLS/bedtools bamtobed -i \${FINAL_BAM_FILE} | awk 'BEGIN{OFS="\\t"}{\$4="N";\$5="1000";print \$0}' | gzip -c > \${FINAL_TA_FILE}

# Subsample tagAlign file
NREADS=15000000
#SUBSAMPLED_TA_FILE="${OFPREFIX}.filt.nodup.sample.\$((\$NREADS/1000000)).SE.tagAlign.gz"
zcat \${FINAL_TA_FILE} | grep -v "chrM" | shuf -n \${NREADS} | gzip -c > $SUBSAMPLED_TA_FILE
PBS

close(SLURMJOB) || die "$!\n";

} #end createSLURM_TagAlign


sub createSLURM_MACS2{
    my @list = @_;
    my $ta = $list[0];

     # Create PBS script to convert SE BAM to tagAlign (BED 3+3 format)
     open(SLURMJOB, " > $OFPREFIX-Encode3Pipeline-MACS2.slurm") || die "$!\n";
     print SLURMJOB <<PBS;
#!/bin/bash
#SBATCH --job-name=$OFPREFIX-Encode3Pipeline-MACS2
#SBATCH --ntasks=16 --mem-per-cpu=16GB
#SBATCH --time=23:59:59
#SBATCH --output=$OFPREFIX-Encode3Pipeline-MACS2.out
#SBATCH --error=$OFPREFIX-Encode3Pipeline-MACS2.err
#
cd \$SLURM_SUBMIT_DIR
cd $OFPREFIX
source /home/rcf-proj/swn/swnVE/bin/activate
#source /usr/usc/python/2.7.8/setup.sh
#R_LIBS_USER=/home/bioinformatic_tools/R-3.2.0/R-3.2.0/library

#=========================
## 2b. Calculate Cross-correlation QC scores

#FINAL_TA_FILE="${FINAL_BAM_PREFIX}.$ta.tagAlign.gz"
FINAL_BAM_FILE="${FINAL_BAM_PREFIX}.bam" 

NTHREADS=16

CC_SCORES_FILE="$SUBSAMPLED_TA_FILE.cc.qc"
CC_PLOT_FILE="$SUBSAMPLED_TA_FILE.cc.plot.pdf"
# CC_SCORE FILE format
# Filename <tab> numReads <tab> estFragLen <tab> corr_estFragLen <tab> PhantomPeak <tab> corr_phantomPeak <tab> argmin_corr <tab> min_corr <tab> phantomPeakCoef <tab> relPhantomPeakCoef <tab> QualityTag
source /usr/usc/R/3.2.5/setup.sh
/usr/usc/R/3.2.5/bin/Rscript /home/bioinformaticstools/phantompeakqualtools/run_spp_nodups.R -c=$SUBSAMPLED_TA_FILE -p=\${NTHREADS} -filtchr=chrM -savp=\${CC_PLOT_FILE} -out=\${CC_SCORES_FILE}
#sed r 's/,[^\t]+//g' \${CC_SCORES_FILE} > temp
#mv temp \${CC_SCORES_FILE}
###############################################
#123456789
# without savr and savd - really fast #
read FRAGLEN <<< \$(sed -r 's\/,[^\\t]+\/\/g' \$CC_SCORES_FILE | awk 'BEGIN { OFS="\\t"} {print \$3}')
#123456789
PEAK_OUTPUT_DIR=$PEAK_OUTPUT_DIR
mkdir \$PEAK_OUTPUT_DIR
CHIP_TA_PREFIX="${FINAL_BAM_PREFIX}.MACS2"
##################################################
#========================
# 6. For Histone Marks
# Peak calling and signal tracks using MACSv2 for histone marks
GENOMESIZE='hs' # for human 
#GENOMESIZE='mm' #for mouse

PEAKTYPE=$peakType

if [ "\$PEAKTYPE" == "narrow" ];
then
    #========= # Generate narrow peaks and preliminary signal tracks ============
    $MACS2/macs2 callpeak -t \${FINAL_BAM_FILE} -c $controlFile -f BAMPE -n \${PEAK_OUTPUT_DIR}/\${CHIP_TA_PREFIX} -g \${GENOMESIZE} -p 1e-2 --nomodel --shift 0 --extsize \${FRAGLEN} --keep-dup all -B --SPMR
    # Sort by Col8 in descending order and replace long peak names in Column 4 with Peak_<peakRank>
    sort -k 8gr,8gr \${PEAK_OUTPUT_DIR}/\${CHIP_TA_PREFIX}_peaks.narrowPeak | awk 'BEGIN{OFS="\\t"}{\$4="Peak_"NR ; print \$0}' | gzip -c > \${PEAK_OUTPUT_DIR}/\${CHIP_TA_PREFIX}.narrowPeak.gz
    # remove additional files
    #rm -f \${PEAK_OUTPUT_DIR}/\${CHIP_TA_PREFIX}_peaks.xls \${PEAK_OUTPUT_DIR}/\${CHIP_TA_PREFIX}_peaks.narrowPeak \${PEAK_OUTPUT_DIR}/\${CHIP_TA_PREFIX}_summits.bed
fi 

if [[ "\$PEAKTYPE" == "broad" || "\$PEAKTYPE" == "gapped" ]]; 
then
    #========= # Generate Broad and Gapped Peaks ================================
    $MACS2/macs2 callpeak -t \${FINAL_BAM_FILE} -c $controlFile -n \${PEAK_OUTPUT_DIR}/\${CHIP_TA_PREFIX} -g \${GENOMESIZE} -p 1e-2 --broad --nomodel --shift 0 --extsize \${FRAGLEN} --keep-dup all -B --SPMR

    if [ "\$PEAKTYPE" == "broad" ];
then
         # Sort by Col8 (for broadPeak) in descending order and replace long peak names in Column 4 with Peak_<peakRank>
         sort -k 8gr,8gr \${PEAK_OUTPUT_DIR}/\${CHIP_TA_PREFIX}_peaks.broadPeak | awk 'BEGIN{OFS="\\t"}{\$4="Peak_"NR ; print \$0}' | gzip -c > \${PEAK_OUTPUT_DIR}/\${CHIP_TA_PREFIX}.broadPeak.gz
fi

    if [ "\$PEAKTYPE" == "gapped" ];
then
         # Sort by Col 14(for gappedPeak) in descending order and replace long peak names in Column 4 with Peak_<peakRank>
         sort -k 14gr,14gr \${PEAK_OUTPUT_DIR}/\${CHIP_TA_PREFIX}_peaks.gappedPeak | awk 'BEGIN{OFS="\\t"}{\$4="Peak_"NR ; print \$0}' | gzip -c > \${PEAK_OUTPUT_DIR}/\${CHIP_TA_PREFIX}.gappedPeak.gz
fi

    # remove additional files
    #rm -f \${PEAK_OUTPUT_DIR}/\${CHIP_TA_PREFIX}_peaks.xls \${PEAK_OUTPUT_DIR}/\${CHIP_TA_PREFIX}_peaks.broadPeak \${PEAK_OUTPUT_DIR}/\${CHIP_TA_PREFIX}_peaks.gappedPeak \${PEAK_OUTPUT_DIR}/\${CHIP_TA_PREFIX}_summits.bed
#===========================================
fi
PBS

close(SLURMJOB) || die "$!\n";

} #end createSLURM_MACS2

sub createSLURM_SPP{

    my @list = @_;
    my $ta = $list[0];

    open(SLURMJOB, " > $OFPREFIX-Encode3Pipeline-SPP.slurm") || die "$!\n";
    print SLURMJOB <<PBS;
#!/bin/bash
#SBATCH --job-name=$OFPREFIX-Encode3Pipeline-SPP
#SBATCH --ntasks=16 --mem-per-cpu=16GB
#SBATCH --time=23:59:59
#SBATCH --output=$OFPREFIX-Encode3Pipeline-SPP.out
#SBATCH --error=$OFPREFIX-Encode3Pipeline-SPP.err
#
cd \$SLURM_SUBMIT_DIR
cd $OFPREFIX
source /home/rcf-proj/swn/swnVE/bin/activate
#source /usr/usc/python/2.7.8/setup.csh
#R_LIBS_USER=/home/bioinformatic_tools/R-3.2.0/R-3.2.0/library

### TRUE REPLICATE ###
#FINAL_TA_FILE="${FINAL_BAM_PREFIX}.$ta.tagAlign.gz"
FINAL_BAM_FILE="${FINAL_BAM_PREFIX}.bam"
CHIP_TA_PREFIX="${OFPREFIX}.filt.nodup.srt.SPP"
PEAK_OUTPUT_DIR=$PEAK_OUTPUT_DIR
mkdir \$PEAK_OUTPUT_DIR
NTHREADS=30
CC_SCORES_FILE="$SUBSAMPLED_TA_FILE.cc.qc"
CC_PLOT_FILE="${SUBSAMPLED_TA_FILE}.cc.plot.pdf"
/usr/usc/R/3.2.5/bin/Rscript /home/rcf-proj/ENCODE3.pipeline.test/bioinformaticstools/phantompeakqualtools/run_spp_nodups.R -c=$SUBSAMPLED_TA_FILE -p=\${NTHREADS} -filtchr=chrM -savp=\${CC_PLOT_FILE} -out=\${CC_SCORES_FILE}
read FRAGLEN <<< \$(sed -r 's\/,[^\\t]+\/\/g' \$CC_SCORES_FILE | awk 'BEGIN { OFS="\\t"} {print \$3}')

## 3a. Peak calling - SPP
/usr/usc/R/3.2.5/bin/Rscript /home/ENCODE3.pipeline.test/bioinformaticstools/phantompeakqualtools/run_spp_nodups.R -c=\${FINAL_TA_FILE} -i=$controlFile -npeak=300000 -odir=\${PEAK_OUTPUT_DIR} -speak=\${FRAGLEN} -savr -savp -rf -out=\${PEAK_OUTPUT_DIR}/\${CHIP_TA_PREFIX}.ccscores -p=\${NTHREADS}

## For tracks call peaks using MACS2 narrowPeak
GENOMESIZE='hs' # for human
#GENOMESIZE='mm' #for mouse
CHIP_TA_PREFIX="${OFPREFIX}.filt.nodup.srt.MACS2"
$MACS2/macs2 callpeak -t \${FINAL_BAM_FILE} -c $controlFile -f BAMPE -n \${PEAK_OUTPUT_DIR}/\${CHIP_TA_PREFIX} -g \${GENOMESIZE} -p 1e-2 --nomodel --shift 0 --extsize \${FRAGLEN} --keep-dup all -B --SPMR
sort -k 8gr,8gr \${PEAK_OUTPUT_DIR}/\${CHIP_TA_PREFIX}_peaks.narrowPeak | awk 'BEGIN{OFS="\\t"}{\$4="Peak_"NR ; print \$0}' | gzip -c > \${PEAK_OUTPUT_DIR}/\${CHIP_TA_PREFIX}.narrowPeak.gz
PBS

close(SLURMJOB) || die "$!\n";

} #end createSLURM_SPP

sub createSLURM_tracks{

    my @list = @_;
    my $ta = $list[0];

    my $CHIP_TA_PREFIX;
    my $PEAK_OUTPUT_DIR;
    if($peakSoftware =~ /MACS2/) {
       $PEAK_OUTPUT_DIR="MACS2_PEAKS";
       $CHIP_TA_PREFIX="${FINAL_BAM_PREFIX}.MACS2";
    } elsif ($peakSoftware =~ /SPP/) {
       $PEAK_OUTPUT_DIR="SPP_PEAKS";
       $CHIP_TA_PREFIX="${FINAL_BAM_PREFIX}.MACS2";
    }

    #Create PBS script to write tracks     
    open(SLURMJOB, " > $OFPREFIX-Encode3Pipeline-tracks.slurm") || die "$!\n";
    print SLURMJOB <<PBS;
#!/bin/bash
#SBATCH --job-name=$OFPREFIX-Encode3Pipeline-tracks
#SBATCH --ntasks=16 --mem-per-cpu=16GB
#SBATCH --time=23:59:59
#SBATCH --output=$OFPREFIX-Encode3Pipeline-tracks.out
#SBATCH --error=$OFPREFIX-Encode3Pipeline-tracks.err
#
cd \$SLURM_SUBMIT_DIR
cd $OFPREFIX
source /home/rcf-proj/swn/swnVE/bin/activate
#source /usr/usc/python/2.7.8/setup.sh
#R_LIBS_USER=/home/bioinformatic_tools/R-3.2.0/R-3.2.0/library

FINAL_TA_FILE="${FINAL_BAM_PREFIX}.$ta.tagAlign.gz"
FINAL_BAM_FILE="${FINAL_BAM_PREFIX}.bam"
PEAK_OUTPUT_DIR=$PEAK_OUTPUT_DIR
CHIP_TA_PREFIX=$CHIP_TA_PREFIX

# For Fold enrichment signal tracks ============================================ 
# This file is a tab delimited file with 2 columns Col1 (chromosome name), Col2 (chromosome size in bp).
$MACS2/macs2 bdgcmp -t \${PEAK_OUTPUT_DIR}/\${CHIP_TA_PREFIX}_treat_pileup.bdg -c \${PEAK_OUTPUT_DIR}/\${CHIP_TA_PREFIX}_control_lambda.bdg --outdir \${PEAK_OUTPUT_DIR} -o \${CHIP_TA_PREFIX}_FE.bdg -m FE
# Remove coordinates outside chromosome sizes (stupid MACS2 bug)
$BEDTOOLS/slopBed -i \${PEAK_OUTPUT_DIR}/\${CHIP_TA_PREFIX}_FE.bdg -g ${CHRSIZEFILE} -b 0 | awk '{if (\$3 != -1) print \$0}' | $UCSC_BROWSER_TOOLS/bedClip stdin ${CHRSIZEFILE} \${PEAK_OUTPUT_DIR}/\${CHIP_TA_PREFIX}.fc.signal.bedgraph
#rm -f \${PEAK_OUTPUT_DIR}/\${CHIP_TA_PREFIX}_FE.bdg
#123456789
# sortBed before bedGraphToBigWig - got error MACS2_PEAKS/SEP044.USC971.H3K27Ac.CON.filt.nodup.srt.MACS2.fc.signal.bedgraph is not case-sensitive sorted at line 45083.  Please use "sort -k1,1 -k2,2n" with LC_COLLATE=C,  or bedSort and try again
LC_COLLATE=C sort -k1,1 -k2,2n \${PEAK_OUTPUT_DIR}/\${CHIP_TA_PREFIX}.fc.signal.bedgraph > \${PEAK_OUTPUT_DIR}/\${CHIP_TA_PREFIX}.fc.signal.sorted.bedgraph
# 123456789
# Convert bedgraph to bigwig using UCSC_Browser_Tools
$UCSC_BROWSER_TOOLS/bedGraphToBigWig \${PEAK_OUTPUT_DIR}/\${CHIP_TA_PREFIX}.fc.signal.sorted.bedgraph ${CHRSIZEFILE} \${PEAK_OUTPUT_DIR}/\${CHIP_TA_PREFIX}.fc.signal.sorted.bw
#rm -f \${PEAK_OUTPUT_DIR}/\${CHIP_TA_PREFIX}.fc.signal.bedgraph
# Convert bedgraph to tdf file using igv tools
/usr/bin/java -Xmx30g -jar $IGVTOOLS/igvtools.jar toTDF -z 7 \${PEAK_OUTPUT_DIR}/\${CHIP_TA_PREFIX}.fc.signal.sorted.bedgraph \${PEAK_OUTPUT_DIR}/\${CHIP_TA_PREFIX}.fc.signal.sorted.tdf hg19 

#=============================== # For -log10(p-value) signal tracks ==============================
# Compute sval = min(no. of reads in ChIP, no. of reads in control) / 1,000,000 
#chipReads=\$(zcat \${FINAL_TA_FILE} | wc -l | awk '{printf "%f", \$1/1000000}');
chipReads=\$(\${FINAL_BAM_FILE} | wc -l | awk '{printf "%f", \$1/1000000}');
# 123456789
controlReads=\$(zcat $controlFile | wc -l | awk '{printf "%f", \$1/1000000}');
sval=\$(echo "\${chipReads} \${controlReads}" | awk '\$1>\$2{printf "%f",\$2} \$1<=\$2{printf "%f",\$1}');
$MACS2/macs2 bdgcmp -t \${PEAK_OUTPUT_DIR}/\${CHIP_TA_PREFIX}_treat_pileup.bdg -c \${PEAK_OUTPUT_DIR}/\${CHIP_TA_PREFIX}_control_lambda.bdg --outdir \${PEAK_OUTPUT_DIR} -o \${CHIP_TA_PREFIX}_ppois.bdg -m ppois -S \${sval}
# Remove coordinates outside chromosome sizes (stupid MACS2 bug)
$BEDTOOLS/slopBed -i \${PEAK_OUTPUT_DIR}/\${CHIP_TA_PREFIX}_ppois.bdg -g ${CHRSIZEFILE} -b 0 | awk '{if (\$3 != -1) print \$0}' | $UCSC_BROWSER_TOOLS/bedClip stdin ${CHRSIZEFILE} \${PEAK_OUTPUT_DIR}/\${CHIP_TA_PREFIX}.pval.signal.bedgraph
rm -rf \${PEAK_OUTPUT_DIR}/\${CHIP_TA_PREFIX}_ppois.bdg
# 123456789
LC_COLLATE=C sort -k1,1 -k2,2n \${PEAK_OUTPUT_DIR}/\${CHIP_TA_PREFIX}.pval.signal.bedgraph > \${PEAK_OUTPUT_DIR}/\${CHIP_TA_PREFIX}.pval.signal.sorted.bedgraph

# 123456789
# Convert bedgraph to bigwig
$UCSC_BROWSER_TOOLS/bedGraphToBigWig \${PEAK_OUTPUT_DIR}/\${CHIP_TA_PREFIX}.pval.signal.sorted.bedgraph ${CHRSIZEFILE} \${PEAK_OUTPUT_DIR}/\${CHIP_TA_PREFIX}.pval.signal.sorted.bw
#rm -f \${PEAK_OUTPUT_DIR}/\${CHIP_TA_PREFIX}.pval.signal.bedgraph 
#rm -f \${PEAK_OUTPUT_DIR}/\${CHIP_TA_PREFIX}_treat_pileup.bdg \${peakFile}_control_lambda.bdg
/usr/bin/java -Xmx30g -jar $IGVTOOLS/igvtools.jar toTDF -z 7 \${PEAK_OUTPUT_DIR}/\${CHIP_TA_PREFIX}.pval.signal.sorted.bedgraph \${PEAK_OUTPUT_DIR}/\${CHIP_TA_PREFIX}.pval.signal.sorted.tdf hg19
PBS

close(SLURMJOB) || die "$!\n";

} #end createSLURM_tracks

sub createSLURM_BWA_PE{
    my @list = @_;
    my $fastq_1 = $list[0];
    my $fastq_2 = $list[1];

# Create PBS script to run BWA alignnment and run post-alignment filtering
   open(SLURMJOB, " > $OFPREFIX-Encode3Pipeline-bwa-pairedEnd.slurm") || die "$!\n";
   print SLURMJOB <<PBS;
#!/bin/bash
#SBATCH --job-name=$OFPREFIX-Encode3Pipeline-bwa-pairedEnd
#SBATCH --ntasks=16 --mem-per-cpu=16GB
#SBATCH --time=23:59:59
#SBATCH --output=$OFPREFIX-Encode3Pipeline-bwa-pairedEnd.out
#SBATCH --error=$OFPREFIX-Encode3Pipeline-bwa-pairedEnd.err
#
cd \$SLURM_SUBMIT_DIR
cd $OFPREFIX
source /home/rcf-proj/swn/swnVE/bin/activate
#source /usr/usc/python/2.7.8/setup.csh
#R_LIBS_USER=/home/bioinformatic_tools/R-3.2.0/R-3.2.0/library

###
# Run BWA alignment
###

FASTQ_FILE_1=$fastq_1
FASTQ_FILE_2=$fastq_2
NTHREADS=16

# map reads to create raw SAM file 
SAI_FILE_1="${OFPREFIX}_1.sai" 
SAI_FILE_2="${OFPREFIX}_2.sai" 
RAW_SAM_FILE="${OFPREFIX}.raw.sam.gz" 
TMP_RAW_BAM="${OFPREFIX}.tmp.raw.srt.bam" 
RAW_BAM_PREFIX="${OFPREFIX}.raw.srt" 
RAW_BAM_FILE="\${RAW_BAM_PREFIX}.bam" # To be stored 
mkdir TMP 
TMP=TMP 
RANDOM=RANDOM 
BADCIGAR_FILE="\${TMP}/badReads\${RANDOM}.tmp" 
RAW_BAM_FILE_MAPSTATS="\${RAW_BAM_PREFIX}.flagstat.qc"
$BWA/bwa aln -q 5 -l 32 -k 2 -t \${NTHREADS} ${BWA_INDEX_NAME} \${FASTQ_FILE_1} > \${SAI_FILE_1}
$BWA/bwa aln -q 5 -l 32 -k 2 -t \${NTHREADS} ${BWA_INDEX_NAME} \${FASTQ_FILE_2} > \${SAI_FILE_2} 
$BWA/bwa sampe ${BWA_INDEX_NAME} \${SAI_FILE_1} \${SAI_FILE_2} \${FASTQ_FILE_1} \${FASTQ_FILE_2} | gzip -c > \${RAW_SAM_FILE} 
rm \${SAI_FILE_1} 
rm \${SAI_FILE_2} 

# Find bad CIGAR read names 
zcat \${RAW_SAM_FILE} | awk 'BEGIN {FS="\\t" ; OFS="\\t"} ! \/^@\/ && \$6!="*" { cigar=\$6;gsub("[09]+D","",cigar); n = split(cigar,vals,"[AZ]"); s = 0; for (i=1;i<=n;i++) s=s+vals[i];seqlen=length(\$10) ; if (s!=seqlen) print \$1"\t" ; }' | sort | uniq > \${BADCIGAR_FILE} 
# Remove bad CIGAR read pairs 
if [[ \$(cat \${BADCIGAR_FILE} | wc -l) -gt 0 ]] 
then 
        zcat \${RAW_SAM_FILE} | grep -v -F -f \${BADCIGAR_FILE} | $SAMTOOLS/samtools view -Su - > \${TMP_RAW_BAM} 
        $SAMTOOLS/samtools sort \${TMP_RAW_BAM} \${RAW_BAM_PREFIX} 
else 
        $SAMTOOLS/samtools view -Su \${RAW_SAM_FILE} | $SAMTOOLS/samtools sort - \${RAW_BAM_PREFIX} 
fi 
rm \${BADCIGAR_FILE} \${TMP_RAW_BAM} \${RAW_SAM_FILE} 
$SAMTOOLS/samtools flagstat \${RAW_BAM_FILE} > \${RAW_BAM_FILE_MAPSTATS} 
#=============================== 

###
## Run Post-alignment filtering
####

FILT_BAM_PREFIX="${OFPREFIX}.filt.srt"
FILT_BAM_FILE="\${FILT_BAM_PREFIX}.bam"
TMP_FILT_UNSORTED_BAM="tmp.\${FILT_BAM_PREFIX}.unsorted.bam"
TMP_FILT_BAM_PREFIX="tmp.\${FILT_BAM_PREFIX}.nmsrt"
TMP_FILT_BAM_FILE="\${TMP_FILT_BAM_PREFIX}.bam"
MAPQ_THRESH=30
$SAMTOOLS/samtools view -F 1804 -f 2 -q \${MAPQ_THRESH} -u \${RAW_BAM_FILE} > \${TMP_FILT_UNSORTED_BAM}
$SAMTOOLS/samtools sort -n \${TMP_FILT_UNSORTED_BAM} \${TMP_FILT_BAM_PREFIX} # Will produce name sorted BAM

# Remove orphan reads (pair was removed)
# and read pairs mapping to different chromosomes
# Obtain position sorted BAM
$SAMTOOLS/samtools fixmate -r \${TMP_FILT_BAM_FILE} ${OFPREFIX}.fixmate.tmp
$SAMTOOLS/samtools view -F 1804 -f 2 -u ${OFPREFIX}.fixmate.tmp > ${OFPREFIX}.fixmate.tmp2.bam
$SAMTOOLS/samtools sort ${OFPREFIX}.fixmate.tmp2.bam \${FILT_BAM_PREFIX}
rm \${TMP_FILT_UNSORTED_BAM}
rm ${OFPREFIX}.fixmate.tmp
rm ${OFPREFIX}.fixmate.tmp2.bam
rm \${TMP_FILT_BAM_FILE}

# Mark duplicates
# use picard-tools/1.92
TMP_FILT_BAM_FILE="\${FILT_BAM_PREFIX}.dupmark.bam"
MARKDUP="$PICARD/MarkDuplicates.jar"
DUP_FILE_QC="\${FILT_BAM_PREFIX}.dup.qc" # QC file
java -Xmx30G -jar \${MARKDUP} INPUT=\${FILT_BAM_FILE} OUTPUT=\${TMP_FILT_BAM_FILE} METRICS_FILE=\${DUP_FILE_QC} VALIDATION_STRINGENCY=LENIENT ASSUME_SORTED=true REMOVE_DUPLICATES=false
mv \${TMP_FILT_BAM_FILE} \${FILT_BAM_FILE}

# remove duplictes and index final position sorted BAM
FINAL_BAM_FILE="${FINAL_BAM_PREFIX}.bam" # To be stored
FINAL_BAM_INDEX_FILE="${FINAL_BAM_PREFIX}.bai" # To be stored
FINAL_BAM_FILE_MAPSTATS="${FINAL_BAM_PREFIX}.flagstat.qc" # QC file
FINAL_NMSRT_BAM_PREFIX="${OFPREFIX}.filt.nmsrt.nodup"
FINAL_NMSRT_BAM_FILE="\${FINAL_NMSRT_BAM_PREFIX}.bam" # To be stored
$SAMTOOLS/samtools view -F 1804 -f 2 -b \${FILT_BAM_FILE} > \${FINAL_BAM_FILE}
$SAMTOOLS/samtools sort -n \${FINAL_BAM_FILE} \${FINAL_NMSRT_BAM_PREFIX}

# Index Final BAM file
$SAMTOOLS/samtools index \${FINAL_BAM_FILE} \${FINAL_BAM_INDEX_FILE}
$SAMTOOLS/samtools flagstat \${FINAL_BAM_FILE} > \${FINAL_BAM_FILE_MAPSTATS}

rm -rf \${TMP}
PBS

close(SLURMJOB) || die "$!\n";
} #end createSLURM_BWA_PE 

sub createSLURM_LibraryComplexity_PE{

    # Create PBS script to compute library complexity
    open(SLURMJOB, " > $OFPREFIX-Encode3Pipeline-computeLibraryComplexity-pairedEnd.slurm") || die "$!\n";
    print SLURMJOB <<PBS;
#!/bin/bash
#SBATCH --job-name=$OFPREFIX-Encode3Pipeline-computeLibraryComplexity-pairedEnd
#SBATCH --ntasks=16 --mem-per-cpu=16GB
#SBATCH --time=23:59:59
#SBATCH --output=$OFPREFIX-Encode3Pipeline-computeLibraryComplexity-pairedEnd.out
#SBATCH --error=$OFPREFIX-Encode3Pipeline-computeLibraryComplexity-pairedEnd.err
#
cd \$SLURM_SUBMIT_DIR
cd $OFPREFIX
source /home/rcf-proj/swn/swnVE/bin/activate
#source /usr/usc/python/2.7.8/setup.sh
#R_LIBS_USER=/home/bioinformatic_tools/R-3.2.0/R-3.2.0/library

# Compute library complexity
# sort by position and strand

FILT_BAM_PREFIX=${OFPREFIX}.filt.srt
FILT_BAM_FILE=\${FILT_BAM_PREFIX}.bam
PBC_FILE_QC="${FINAL_BAM_PREFIX}.pbc.qc"

$SAMTOOLS/samtools sort -n \${FILT_BAM_FILE} ${OFPREFIX}.srt.tmp
$BEDTOOLS/bedtools bamtobed -bedpe -i ${OFPREFIX}.srt.tmp.bam | awk 'BEGIN{OFS="\\t"}{print \$1,\$2,\$4,\$6,\$9,\$10}' | grep -v 'chrM' | sort | uniq -c | awk 'BEGIN{mt=0;m0=0;m1=0;m2=0} (\$1==1){m1=m1+1} (\$1==2){m2=m2+1} {m0=m0+1} {mt=mt+\$1} END{printf "%d\\t%d\\t%d\\t%d\\t%f\\t%f\\t%f\\n",mt,m0,m1,m2,m0/mt,m1/m0,m1/m2}' > \${PBC_FILE_QC}

rm ${OFPREFIX}.srt.tmp.bam
rm \${FILT_BAM_FILE}
# PBC File output
# # TotalReadPairs [tab] DistinctReadPairs [tab] OneReadPair [tab] TwoReadPairs [tab] NRF=Distinct/Total [tab] PBC1=OnePair/Distinct [tab] P
# BC2=OnePair/TwoPair
PBS

close(SLURMJOB) || die "$!\n";
} #end createSLURM_LibraryComplexity_PE 

sub createSLURM_TagAlign_PE{

    my $NREADS=15000000;
    my $f = $NREADS/1000000;
    $SUBSAMPLED_TA_FILE="${OFPREFIX}.filt.nodup.sample.$f.MATE1.tagAlign.gz";

    open(SLURMJOB, " > $OFPREFIX-Encode3Pipeline-tagAlign-pairedEnd.slurm") || die "$!\n";
    print SLURMJOB <<PBS;
#!/bin/bash
#SBATCH --job-name=$OFPREFIX-Encode3Pipeline-tagAlign-pairedEnd
#SBATCH --ntasks=16 --mem-per-cpu=16GB
#SBATCH --time=23:59:59
#SBATCH --output=$OFPREFIX-Encode3Pipeline-tagAlign-pairedEnd.out
#SBATCH --error=$OFPREFIX-Encode3Pipeline-tagAlign-pairedEnd.err
#
cd \$SLURM_SUBMIT_DIR
cd $OFPREFIX
source /home/rcf-proj/swn/swnVE/bin/activate
#source /usr/usc/python/2.7.8/setup.sh
#R_LIBS_USER=/home/bioinformatic_tools/R-3.2.0/R-3.2.0/library

#===========================
#### 2a. convert BAM to tagAlign (BED 3+3 format)##
#
# bedtools 2.19.1, gawk shuf
# module add bedtools/2.19.1

# Create tagAlign file containing both read pairs
FINAL_TA_FILE="${FINAL_BAM_PREFIX}.PE2SE.tagAlign.gz"
FINAL_BAM_FILE="${FINAL_BAM_PREFIX}.bam"

$BEDTOOLS/bedtools bamtobed -i \${FINAL_BAM_FILE} | awk 'BEGIN{OFS="\\t"}{\$4="N";\$5="1000";print \$0}' | gzip -c > \${FINAL_TA_FILE}

#=============
# Create BEDPE file 
#=============
FINAL_NMSRT_BAM_PREFIX="${OFPREFIX}.filt.nmsrt.nodup"
FINAL_NMSRT_BAM_FILE="\${FINAL_NMSRT_BAM_PREFIX}.bam" 
FINAL_BEDPE_FILE="\${FINAL_NMSRT_BAM_PREFIX}.bedpe.gz"
$BEDTOOLS/bedtools bamtobed -bedpe -mate1 -i \${FINAL_NMSRT_BAM_FILE} | gzip -c > \${FINAL_BEDPE_FILE}

# =================================
# Subsample tagAlign file
NREADS=$NREADS
#SUBSAMPLED_TA_FILE="${OFPREFIX}.filt.nodup.sample.\$((\$NREADS/1000000)).MATE1.tagAlign.gz"
zcat \${FINAL_BEDPE_FILE} | grep -v "chrM" | shuf -n \${NREADS} | awk 'BEGIN{OFS="\\t"}{print \$1,\$2,\$3,"N","1000",\$9}' | gzip -c > $SUBSAMPLED_TA_FILE
PBS

close(SLURMJOB) || die "$!\n";

} #end createSLURM_TagAlign_PE

sub createSLURM_GetStats{

    open(SLURMJOB, " > $OFPREFIX-Encode3Pipeline-GetStats.slurm") || die "$!\n";
    print SLURMJOB <<PBS;
#!/bin/bash
#SBATCH --job-name=$OFPREFIX-Encode3Pipeline-GetStats
#SBATCH --ntasks=16 --mem-per-cpu=16GB
#SBATCH --time=23:59:59
#SBATCH --output=$OFPREFIX-Encode3Pipeline-GetStats.out
#SBATCH --error=$OFPREFIX-Encode3Pipeline-GetStats.err
#
cd \$SLURM_SUBMIT_DIR
cd $OFPREFIX
source /home/rcf-proj/swn/swnVE/bin/activate
#source /usr/usc/python/2.7.8/setup.sh
#R_LIBS_USER=/home/bioinformatic_tools/R-3.2.0/R-3.2.0/library

TOTAL_READ_COUNT="\$(zcat $fastq_1 | echo \$((`wc -l`/4)))"
UNIQ_COUNT="\$(cat $OFPREFIX.filt.nodup.srt.flagstat.qc | awk 'NR==1{map=\$1}END{print map}')"
PBC_STATS="\$(cat $OFPREFIX.filt.nodup.srt.pbc.qc | awk '{print \$1" "\$2" "\$5" "\$6" "\$7}')"
CCQC_STATS="\$(cat *tagAlign.gz.cc.qc  | awk '{print \$3" "\$9" "\$10}' | awk '{split(\$1,a,","); print a[1]" "\$2" "\$3}')" 
echo "TOTAL_READS (from fastq)\tUNIQ_READS\tTOTAL_READS(PBC_STATS)\tDISTINCT_READS(PBC_STATS)\tNRF=DISTINCT/TOTAL(PBC_STATS)\tPBC1=OnePair/Distinct\tPBC2=OnePair/TwoPair\tEST_FRAG_LENGTH\tPHANTOM_PEAK_COEF(CCQC_STATS)\tREL_PHANTOM_PEAK_COEF(CCQC_STATS)" >> $OFPREFIX-PipelineStats
echo "\$TOTAL_READ_COUNT	\$UNIQ_COUNT	\$PBC_STATS	\$CCQC_STATS" >> $OFPREFIX-PipelineStats
PBS

close(SLURMJOB) || die "$!\n";

} #end createSLURM_GetStats

