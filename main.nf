#!/usr/bin/env nextflow


//help function for the tool
def show_help (){
  log.info digenoma_Header()
  log.info tool_header()
    log.info"""

    Usage:

    The typical command for running the pipeline is as follows:

    nextflow run digenoma-lab/k-count -r v1.0 --reads '*_R{1,2}.fastq.gz'  -profile singularity

    Mandatory arguments:
      --reads [file]                Path to input data

    Input alternatives:
      --reads_csv                   file with tabular data for each sample to process [sampleID fwd_path rev_path]      --svs                         file with tabular data for each sample with structural variants in bedpe format [sampleID sv_path]
      --kmer                        k-mer size [def: 21]
      --ploidy                      ploidy of the sample [def: 2]
      --remove_kmc_db               Remove kmc k-mer databases [def: false]
      Run examples:

      nextflow run  digenoma-lab/k-count -r v1.0 -singularity --reads 'test_dataset/reads/*.R{1,2}.fastq.gz' --ref_fa test_dataset/genome.fa --ref_gtf test_dataset/genome.gtf
      or
      nextflow run  digenoma-lab/k-count -r v1.0 -singularity --bams '/path/to/bams/' --ref_fa test_dataset/genome.fa --ref_gtf test_dataset/genome.gtf

      """.stripIndent()
}

// we star coding the pipeline

// Show help message
if (params.help) exit 0, show_help()

//print header and tool name
log.info digenoma_Header()
log.info tool_header()


//log.info params.reads_csv
mode="reads"
//expect a file with header "sampleID fwd_path rev_path"
//see file ./test_dataset/sample_fwrev.txt
if(params.reads_csv) {
      Channel.fromPath(file(params.reads_csv)).splitCsv(header: true, sep: '\t', strip: true)
                      .map{row -> [ row.sampleID, [file(row.fwd), file(row.rev)]]}
                      .ifEmpty{exit 1, "params.reads_csv was empty - no input files supplied" }
                      .set{reads_kmc}

}else{
    //expect a regular expresion like '*_{1,2}.fastq.gz'
    Channel.fromFilePairs(params.reads, size: 2 )
        .ifEmpty{exit 1, "Cannot find any reads matching: ${params.reads}\nNB: Path needs to be enclosed in quotes!" }
        .set{reads_kmc}
}




//funtion to count k-mers for a set of reads
process ckmers{
  tag "${sample}"
  cpus params.cpu
  memory params.mem+'G'
  //we can remove this to don't keep the bam files
  publishDir "${params.outdir}/kmc", mode: 'copy'
  //publishDir "$params.outdir/$sampleId/counts", pattern: "*_counts.txt"

  input:
      set val(sample), file(reads) from reads_kmc
  output:
      //star bam files
      set val(sample), file("${sample}.hist") into gs_in
      //star mapping stats and gene counts *.{tsv,txt}
      set val(sample), file("${sample}.{kmc_pre,kmc_suf}") optional true into kmc_output
  //when: !(params.bams)

  script:
  def del_kmc=""
  if(params.remove_kmc_db){
    del_kmc="rm -f ${sample}.kmc_pre ${sample}.kmc_suf"
  }
  if(params.debug){
    """
    touch ${sample}.hist
    touch ${sample}.kmc_pre
    touch ${sample}.kmc_suf
    ${del_kmc}
    """
  //BAMs are given as input
  }else{
  """
  ls ${reads} > reads.txt
  mkdir ${sample}-tmp
  # We count the 21 mers in 12Gb memory
  kmc -k${params.kmer} -t${params.cpu} -m${params.mem} -ci1 -cs10000 @reads.txt ${sample} ${sample}-tmp
  # We create the output histogram
  kmc_tools transform ${sample} histogram ${sample}.hist -cx10000
  ${del_kmc}
  """
  //normal reads are given as input
  }
}

//}
/*
 * run arriba fusion with genomic SVs
 * In case of the Variant Call Format, the file must comply with the VCF specification for structural variants.
 * In particular, Arriba requires that the SVTYPE field be present in the INFO column and specify one of the four values BND, DEL, DUP, INV.
 * In addition, for all SVTYPEs other than BND, the END field must be present and specify the second breakpoint of the structural variant.
 * Structural variants with single breakends are silently ignored.
*/
process genomescope {
    tag "${sample}-GS"
    cpus 1
    memory params.mem+'G'

    publishDir "${params.outdir}/genomescope/", mode: 'copy'

    input:
        set sample, file(hist) from gs_in
    output:
        path("${sample}-GS")

    script:
    if(params.debug){
      """
      mkdir ${sample}-GS
      touch ${sample}-GS/${sample}-GS.png
      """
    }else{
    """
    genomescope2 -i ${hist} -o ${sample}-GS -k ${params.kmer}
    """
  }
}

//header for the IARC tools
// the logo was generated using the following page
// http://patorjk.com/software/taag  (ANSI logo generator)
def digenoma_Header (){
     return  """

██████  ██  ██████  ███████ ███    ██  ██████  ███    ███  █████      ██       █████  ██████
██   ██ ██ ██       ██      ████   ██ ██    ██ ████  ████ ██   ██     ██      ██   ██ ██   ██
██   ██ ██ ██   ███ █████   ██ ██  ██ ██    ██ ██ ████ ██ ███████     ██      ███████ ██████
██   ██ ██ ██    ██ ██      ██  ██ ██ ██    ██ ██  ██  ██ ██   ██     ██      ██   ██ ██   ██
██████  ██  ██████  ███████ ██   ████  ██████  ██      ██ ██   ██     ███████ ██   ██ ██████

Nextflow pipelines for genomics.

"""
}

//this use ANSI colors to make a short tool description
//useful url: http://www.lihaoyi.com/post/BuildyourownCommandLinewithANSIescapecodes.html
def tool_header (){
        return """
        K-count: (${workflow.manifest.version})
        """
}
