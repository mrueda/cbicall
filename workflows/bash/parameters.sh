#   $VERSION taken from CBICall

DATADIR=/media/mrueda/2TBS
DBDIR=$DATADIR/Databases
NGSUTILS=$DATADIR/NGSutils

# ENV
export TMPDIR=$DATADIR/tmp
export LC_ALL=C
export GATK_DISABLE_AUTO_S3_UPLOAD=true # does not work

MEM=8G
ARCH=$(uname -m)  
#JAVA=/usr/bin/java # Java 9
#module load java/1.7.0_21

if [ "$ARCH" == "aarch64" ]; then
    JAVA=/usr/lib/jvm/java-8-openjdk-arm64/bin/java
    BWA=$NGSUTILS/bwa-0.7.18_arm64/bwa           # Needs ~6g RAM
    SAM=$NGSUTILS/samtools-0.1.19_arm64/samtools # x4 faster than v1.3
    BED=$NGSUTILS/bedtools2_arm64/bin/bedtools

else
    JAVA=/usr/lib/jvm/java-8-openjdk-amd64/bin/java
    BWA=$NGSUTILS/bwa-0.7.18/bwa           # Needs ~6g RAM
    SAM=$NGSUTILS/samtools-0.1.19/samtools # x4 faster than v1.3
    BED=$NGSUTILS/bedtools2/bin/bedtools
fi

PIC="$JAVA  -Xmx$MEM -Djava.io.tmpdir=$TMPDIR -jar $NGSUTILS/picard/build/libs/picard.jar"
GATK="$JAVA -Xmx$MEM -Djava.io.tmpdir=$TMPDIR -jar $NGSUTILS/gatk/3.5/GenomeAnalysisTK.jar"
MTOOLBOXDIR=$NGSUTILS/MToolBox-master/MToolBox

# GATK bundle, human genome hg19
BUNDLE=$DBDIR/GATK_bundle/b37
REF=$BUNDLE/references_b37_Homo_sapiens_assembly19.fasta
REFGZ=$BUNDLE/references_b37_Homo_sapiens_assembly19.fasta.gz
dbSNP=$DBDIR/dbSNP/human_9606_b144_GRCh37p13/All_20160408.vcf.gz
MILLS_INDELS=$BUNDLE/b37_Mills_and_1000G_gold_standard.indels.b37.vcf.gz
KG_INDELS=$BUNDLE/b37_1000G_phase1.indels.b37.vcf.gz
HAPMAP=$BUNDLE/b37_hapmap_3.3.b37.vcf.gz
OMNI=$BUNDLE/b37_1000G_omni2.5.b37.vcf.gz

# training sets for variant recalibration
SNP_RES="-resource:hapmap,known=false,training=true,truth=true,prior=15.0 $HAPMAP \
         -resource:omni,known=false,training=true,truth=false,prior=12.0 $OMNI \
         -resource:dbsnp,known=true,training=false,truth=false,prior=6.0 $dbSNP "
INDEL_RES="-resource:mills,known=true,training=true,truth=true,prior=12.0 $MILLS_INDELS "

# Agilent SureSelect Whole Exome
EXOM=$DBDIR/Agilent_SureSelect/hg19/bed

# parameters for UnifiedGenotyper
# down sampling, default=250
DCOV=1000
# call_conf,emit_conf, default=30
UG_CALL=50
UG_EMIT=10
