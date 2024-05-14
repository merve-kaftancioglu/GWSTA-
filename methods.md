the bash scripts: 
""" !/bin/bash
SBATCH --job-name=bwa_single
SBATCH --nodes=1
SBATCH --account=users
SBATCH --qos=cuda
SBATCH --partition=cuda
SBATCH --cpus-per-task=16
SBATCH -o /cta/users/merve.kaftancioglu/alternative_scenario/bwa_as/out-err/%x-%j-slurm.out
SBATCH -e /cta/users/merve.kaftancioglu/alternative_scenario/bwa_as/out-err/%x-%j-slurm.err

 Load BWA module
module load bwa-0.7.17-gcc-9.2.0-xwkclqy

 Define directories
input_dir="/cta/users/merve.kaftancioglu/alternative_scenario/trimmomatic_as/output"
output_dir="/cta/users/merve.kaftancioglu/alternative_scenario/bwa_as/output"
reference_genome="/cta/users/merve.kaftancioglu/alternative_scenario/bwa_as/reference/human_g1k_v37.fasta"

 Check if BWA index exists
if [ ! -f "$reference_genome.bwt" ]; then
  echo "Error: BWA index not found. Please create the index using: bwa mem index $reference_genome"
  exit 1
fi

sample_ids=("h3k27ac_ms_b_ENCSR295QZX" "h3k27ac_ms_b_ENCSR538URI" "h3k27ac_ms_b_ENCSR617GMR" "h3k4me1_ms_b_ENCSR084FXT" "h3k4me1_ms_b_ENCSR480ITK" "h3k4me1_ms_b_ENCSR759DHN" "h3k4me1_ms_b_ENCSR900IAM" "h3k4me1_ms_b_ENCSR931JYC" "h3k4me3_ms_b_ENCSR260CRI" "h3k4me3_ms_b_ENCSR461QMZ" "h3k4me3_ms_b_ENCSR530AVK" "h3k4me3_ms_b_ENCSR713ZYF" "h3k4me3_ms_b_ENCSR923EGO" "h3k9me3_ms_b_ENCSR445DNH" "h3k9me3_ms_b_ENCSR486SZN" "h3k9me3_ms_b_ENCSR983MKX" "h3k27ac_ms_b_ENCSR364OIK" "h3k27ac_ms_b_ENCSR541PLO" "h3k36me3_ms_cd4_ENCSR471WZE" "h3k36me3_ms_cd4_ENCSR473TGK" "h3k4me1_ms_nk_ENCSR482UGV" "h3k4me1_ms_nk_ENCSR718LCI" "h3k4me1_ms_nk_ENCSR806HUT" "h3k4me3_ms_nk_ENCSR234BAD" "h3k4me3_ms_nk_ENCSR703TWY" "h3k4me3_ms_nk_ENCSR755ZGS" "h3k9me3_ms_nk_ENCSR061ATV" "h3k9me3_ms_nk_ENCSR611YIA" "h3k9me3_ms_nk_ENCSR667HUI" "h3k27ac_ms_nk_ENCSR246TTM")

 Loop through all sample subdirectories
for sample_dir in "$input_dir"/*; do
   Extract sample ID from directory name (assuming it matches)
  sample_id=$(basename "$sample_dir")

   Check if it's a directory (skip non-directories)
  if [[ -d "$sample_dir" ]]; then
    echo "Processing sample: $sample_id"

     Create output directory for the sample (if it doesn't exist)
    mkdir -p "$output_dir/$sample_id"

     Find FASTQ files within the sample directory
     Adjust the pattern below if your file names differ
    fastq_files=( "$sample_dir"/*.fastq.gz )

     Check if any FASTQ files found
    if [[ ${fastq_files[@]} -eq 0 ]]; then
      echo "Warning: No FASTQ files found in $sample_dir"
      continue  Skip to next iteration if no FASTQ files
    fi

     BWA mem command with paired-end or single-end handling
    if [[ ${fastq_files[@]} -eq 2 ]]; then
       Paired-end reads (assuming *_1.fastq.gz and *_2.fastq.gz naming)
      bwa mem $reference_genome ${fastq_files[0]} ${fastq_files[1]} > "$output_dir/$sample_id.sam"
    else
       Single-end read (assuming only one FASTQ file)
       Extract only the sample name from the directory path
      sample_name=$(basename "$sample_dir" /)
      bwa mem $reference_genome ${fastq_files[0]} > "$output_dir/$sample_name.sam"
    fi
  fi
done

""" !/bin/bash

 Load STAR module (if applicable)
module load star-2.7.0e-gcc-9.2.0-vynasg3   Replace with your module name (if used)

 Define function to check STAR binary
check_star_binary() {
   Try to locate the STAR executable
  star_path=$(which star 2>/dev/null)

   Check if STAR was found
  if [[ -z "$star_path" ]]; then
    echo "ERROR: STAR binary not found!"
    echo "Please ensure STAR is installed and accessible in your environment."
    exit 1   Exit script with error code
  else
    echo "STAR binary found at: $star_path"
  fi
}

 Check for STAR binary
check_star_binary


""" !/bin/bash
SBATCH --job-name=bwa_single
SBATCH --nodes=1
SBATCH --account=users
SBATCH --qos=cuda
SBATCH --partition=cuda
SBATCH --cpus-per-task=16
SBATCH -o /cta/users/merve.kaftancioglu/alternative_scenario/bwa_as/out-err/%x-%j-slurm.out
SBATCH -e /cta/users/merve.kaftancioglu/alternative_scenario/bwa_as/out-err/%x-%j-slurm.err

 Load BWA module
module load bwa-0.7.17-gcc-9.2.0-xwkclqy

 Define directories
input_dir="/cta/users/merve.kaftancioglu/alternative_scenario/trimmomatic_as/output"
output_dir="/cta/users/merve.kaftancioglu/alternative_scenario/bwa_as/output"
reference_genome="/cta/users/merve.kaftancioglu/alternative_scenario/bwa_as/reference/human_g1k_v37.fasta"

 Check if BWA index exists
if [ ! -f "$reference_genome.bwt" ]; then
  echo "Error: BWA index not found. Please create the index using: bwa mem index $reference_genome"
  exit 1
fi

sample_ids=("h3k27ac_ms_nk_ENCSR469CFL" "h3k27ac_ms_nk_ENCSR746AIX" "h3k27me3_ms_b_ENCSR009ZRH" "h3k27me3_ms_b_ENCSR182NLA" "h3k27me3_ms_b_ENCSR272YVX" "h3k27me3_ms_b_ENCSR649FUX" "h3k27me3_ms_b_ENCSR842WWX" "h3k36me3_ms_b_ENCSR089VQL" "h3k36me3_ms_b_ENCSR119PSR" "h3k36me3_ms_b_ENCSR238WFK" "h3k36me3_ms_b_ENCSR987OPY" "h3k4me1_ms_cd4_ENCSR036YVP" "h3k4me1_ms_cd4_ENCSR043JIA" "h3k4me1_ms_cd4_ENCSR093VUP" "h3k4me1_ms_cd4_ENCSR124IGC" "h3k4me1_ms_cd4_ENCSR315MYP" "h3k4me1_ms_cd4_ENCSR458QEY" "h3k4me1_ms_cd4_ENCSR485FBT" "h3k4me1_ms_cd4_ENCSR641ZFV" "h3k4me1_ms_cd4_ENCSR687CQX" "h3k4me1_ms_cd8_ENCSR231XAP" "h3k4me1_ms_cd8_ENCSR572XTB" "h3k4me1_ms_cd8_ENCSR788TEF" "h3k4me1_ms_cd8_ENCSR815YZL" "h3k4me1_ms_cd8_ENCSR861FEC")

 Loop through all sample subdirectories
for sample_dir in "$input_dir"/*; do
   Extract sample ID from directory name (assuming it matches)
  sample_id=$(basename "$sample_dir")

   Check if it's a directory (skip non-directories)
  if [[ -d "$sample_dir" ]]; then
    echo "Processing sample: $sample_id"

     Create output directory for the sample (if it doesn't exist)
    mkdir -p "$output_dir/$sample_id"

     Find FASTQ files within the sample directory
     Adjust the pattern below if your file names differ
    fastq_files=( "$sample_dir"/*.fastq.gz )

     Check if any FASTQ files found
    if [[ ${fastq_files[@]} -eq 0 ]]; then
      echo "Warning: No FASTQ files found in $sample_dir"
      continue  Skip to next iteration if no FASTQ files
    fi

     BWA mem command with paired-end or single-end handling
    if [[ ${fastq_files[@]} -eq 2 ]]; then
       Paired-end reads (assuming *_1.fastq.gz and *_2.fastq.gz naming)
      bwa mem $reference_genome ${fastq_files[0]} ${fastq_files[1]} > "$output_dir/$sample_id.sam"
    else
       Single-end read (assuming only one FASTQ file)
       Extract only the sample name from the directory path
      sample_name=$(basename "$sample_dir" /)
      bwa mem $reference_genome ${fastq_files[0]} > "$output_dir/$sample_name.sam"
    fi
  fi
done

""" !/bin/bash    
files_to_merge=("/home/merve/Downloads/h3k27ac_ms_b_ENCSR295QZX_links.txt"      "/home/merve/Downloads/h3k27ac_ms_b_ENCSR538URI_links.txt"      "/home/merve/Downloads/h3k27ac_ms_b_ENCSR617GMR_links.txt"      "/home/merve/Downloads/h3k4me1_ms_b_ENCSR084FXT_links.txt"      "/home/merve/Downloads/h3k4me1_ms_b_ENCSR480ITK_links.txt"      "/home/merve/Downloads/h3k4me1_ms_b_ENCSR759DHN_links.txt"      "/home/merve/Downloads/h3k4me1_ms_b_ENCSR900IAM_links.txt"      "/home/merve/Downloads/h3k4me1_ms_b_ENCSR931JYC_links.txt"      "/home/merve/Downloads/h3k4me3_ms_b_ENCSR260CRI_links.txt"      "/home/merve/Downloads/h3k4me3_ms_b_ENCSR461QMZ_links.txt"      "/home/merve/Downloads/h3k4me3_ms_b_ENCSR530AVK_links.txt"      "/home/merve/Downloads/h3k4me3_ms_b_ENCSR713ZYF_links.txt"      "/home/merve/Downloads/h3k4me3_ms_b_ENCSR923EGO_links.txt"      "/home/merve/Downloads/h3k9me3_ms_b_ENCSR445DNH_links.txt"      "/home/merve/Downloads/h3k9me3_ms_b_ENCSR486SZN_links.txt"      "/home/merve/Downloads/h3k9me3_ms_b_ENCSR983MKX_links.txt"      "/home/merve/Downloads/h3k27ac_ms_b_ENCSR364OIK_links.txt"      "/home/merve/Downloads/h3k27ac_ms_b_ENCSR541PLO_links.txt"      "/home/merve/Downloads/h3k36me3_ms_cd4_ENCSR471WZE_links.txt"      "/home/merve/Downloads/h3k36me3_ms_cd4_ENCSR473TGK_links.txt"      "/home/merve/Downloads/h3k4me1_ms_nk_ENCSR482UGV_links.txt"      "/home/merve/Downloads/h3k4me1_ms_nk_ENCSR718LCI_links.txt"      "/home/merve/Downloads/h3k4me1_ms_nk_ENCSR806HUT_links.txt"      "/home/merve/Downloads/h3k4me3_ms_nk_ENCSR234BAD_links.txt"      "/home/merve/Downloads/h3k4me3_ms_nk_ENCSR703TWY_links.txt"      "/home/merve/Downloads/h3k4me3_ms_nk_ENCSR755ZGS_links.txt"      "/home/merve/Downloads/h3k9me3_ms_nk_ENCSR061ATV_links.txt"      "/home/merve/Downloads/h3k9me3_ms_nk_ENCSR611YIA_links.txt"      "/home/merve/Downloads/h3k9me3_ms_nk_ENCSR667HUI_links.txt"      "/home/merve/Downloads/h3k27ac_ms_nk_ENCSR246TTM_links.txt"      "/home/merve/Downloads/h3k27ac_ms_nk_ENCSR469CFL_links.txt"      "/home/merve/Downloads/h3k27ac_ms_nk_ENCSR746AIX_links.txt"      "/home/merve/Downloads/h3k27me3_ms_b_ENCSR009ZRH_links.txt"      "/home/merve/Downloads/h3k27me3_ms_b_ENCSR182NLA_links.txt"      "/home/merve/Downloads/h3k27me3_ms_b_ENCSR272YVX_links.txt"      "/home/merve/Downloads/h3k27me3_ms_b_ENCSR649FUX_links.txt"      "/home/merve/Downloads/h3k27me3_ms_b_ENCSR842WWX_links.txt"      "/home/merve/Downloads/h3k36me3_ms_b_ENCSR089VQL_links.txt"      "/home/merve/Downloads/h3k36me3_ms_b_ENCSR119PSR_links.txt"      "/home/merve/Downloads/h3k36me3_ms_b_ENCSR238WFK_links.txt"      "/home/merve/Downloads/h3k36me3_ms_b_ENCSR987OPY_links.txt"      "/home/merve/Downloads/h3k4me1_ms_cd4_ENCSR036YVP_links.txt"      "/home/merve/Downloads/h3k4me1_ms_cd4_ENCSR043JIA_links.txt"      "/home/merve/Downloads/h3k4me1_ms_cd4_ENCSR093VUP_links.txt"      "/home/merve/Downloads/h3k4me1_ms_cd4_ENCSR124IGC_links.txt"      "/home/merve/Downloads/h3k4me1_ms_cd4_ENCSR315MYP_links.txt"      "/home/merve/Downloads/h3k4me1_ms_cd4_ENCSR458QEY_links.txt"      "/home/merve/Downloads/h3k4me1_ms_cd4_ENCSR485FBT_links.txt"      "/home/merve/Downloads/h3k4me1_ms_cd4_ENCSR641ZFV_links.txt"      "/home/merve/Downloads/h3k4me1_ms_cd4_ENCSR687CQX_links.txt"      "/home/merve/Downloads/h3k4me1_ms_cd8_ENCSR231XAP_links.txt"      "/home/merve/Downloads/h3k4me1_ms_cd8_ENCSR572XTB_links.txt"      "/home/merve/Downloads/h3k4me1_ms_cd8_ENCSR788TEF_links.txt"      "/home/merve/Downloads/h3k4me1_ms_cd8_ENCSR815YZL_links.txt"      "/home/merve/Downloads/h3k4me1_ms_cd8_ENCSR861FEC_links.txt"      "/home/merve/Downloads/h3k4me1_ms_cd8_ENCSR940PHE_links.txt"      "/home/merve/Downloads/h3k4me3_ms_cd4_ENCSR180NCM_links.txt"      "/home/merve/Downloads/h3k4me3_ms_cd4_ENCSR341QLC_links.txt"      "/home/merve/Downloads/h3k4me3_ms_cd4_ENCSR482TGI_links.txt"      "/home/merve/Downloads/h3k4me3_ms_cd4_ENCSR486XJK_links.txt"      "/home/merve/Downloads/h3k4me3_ms_cd4_ENCSR496LKR_links.txt"      "/home/merve/Downloads/h3k4me3_ms_cd4_ENCSR603LTN_links.txt"      "/home/merve/Downloads/h3k4me3_ms_cd4_ENCSR802MXQ_links.txt"      "/home/merve/Downloads/h3k4me3_ms_cd4_ENCSR878YHM_links.txt"      "/home/merve/Downloads/h3k4me3_ms_cd4_ENCSR954ZLD_links.txt"      "/home/merve/Downloads/h3k4me3_ms_cd8_ENCSR231ZZH_links.txt"      "/home/merve/Downloads/h3k4me3_ms_cd8_ENCSR278QHR_links.txt"      "/home/merve/Downloads/h3k4me3_ms_cd8_ENCSR516CKJ_links.txt"      "/home/merve/Downloads/h3k4me3_ms_cd8_ENCSR535YYH_links.txt"      "/home/merve/Downloads/h3k4me3_ms_cd8_ENCSR741XAE_links.txt"      "/home/merve/Downloads/h3k4me3_ms_cd8_ENCSR848XJL_links.txt"      "/home/merve/Downloads/h3k9me3_ms_cd4_ENCSR057IZD_links.txt"      "/home/merve/Downloads/h3k9me3_ms_cd4_ENCSR550DPT_links.txt"      "/home/merve/Downloads/h3k9me3_ms_cd4_ENCSR677OEF_links.txt"      "/home/merve/Downloads/h3k9me3_ms_cd4_ENCSR729ZXH_links.txt"      "/home/merve/Downloads/h3k9me3_ms_cd4_ENCSR851RJV_links.txt"      "/home/merve/Downloads/h3k9me3_ms_cd4_ENCSR919DFZ_links.txt"      "/home/merve/Downloads/h3k9me3_ms_cd4_ENCSR953STZ_links.txt"      "/home/merve/Downloads/h3k9me3_ms_cd4_ENCSR959VZU_links.txt"      "/home/merve/Downloads/h3k9me3_ms_cd8_ENCSR101USF_links.txt"      "/home/merve/Downloads/h3k9me3_ms_cd8_ENCSR354GNT_links.txt"      "/home/merve/Downloads/h3k9me3_ms_cd8_ENCSR377QYB_links.txt"      "/home/merve/Downloads/h3k9me3_ms_cd8_ENCSR733LCG_links.txt"      "/home/merve/Downloads/h3k9me3_ms_cd8_ENCSR980KTW_links.txt"      "/home/merve/Downloads/h3k27ac_ms_cd4_ENCSR200SSJ_links.txt"      "/home/merve/Downloads/h3k27ac_ms_cd4_ENCSR322MTA_links.txt"      "/home/merve/Downloads/h3k27ac_ms_cd4_ENCSR331WMS_links.txt"      "/home/merve/Downloads/h3k27ac_ms_cd4_ENCSR350UKV_links.txt"      "/home/merve/Downloads/h3k27ac_ms_cd4_ENCSR474PYR_links.txt"      "/home/merve/Downloads/h3k27ac_ms_cd4_ENCSR520QDR_links.txt"      "/home/merve/Downloads/h3k27ac_ms_cd4_ENCSR540XNK_links.txt"      "/home/merve/Downloads/h3k27ac_ms_cd4_ENCSR705VSO_links.txt"      "/home/merve/Downloads/h3k27ac_ms_cd4_ENCSR832UMM_links.txt"      "/home/merve/Downloads/h3k27ac_ms_cd8_ENCSR078ATS_links.txt"      "/home/merve/Downloads/h3k27ac_ms_cd8_ENCSR348YRH_links.txt"      "/home/merve/Downloads/h3k27ac_ms_cd8_ENCSR458TOW_links.txt"      "/home/merve/Downloads/h3k27ac_ms_cd8_ENCSR476IPR_links.txt"      "/home/merve/Downloads/h3k27ac_ms_cd8_ENCSR787HDF_links.txt"      "/home/merve/Downloads/h3k27ac_ms_cd8_ENCSR923JIB_links.txt"      "/home/merve/Downloads/h3k27me3_ms_nk_ENCSR469QVG_links.txt"      "/home/merve/Downloads/h3k27me3_ms_nk_ENCSR565WDW_links.txt"      "/home/merve/Downloads/h3k36me3_ms_nk_ENCSR158VSE_links.txt"      "/home/merve/Downloads/h3k36me3_ms_nk_ENCSR245KON_links.txt"      "/home/merve/Downloads/h3k36me3_ms_nk_ENCSR530YDY_links.txt"      "/home/merve/Downloads/h3k36me_ms_cd4_ENCSR532PXR_links.txt"      "/home/merve/Downloads/h3k27me3_ms_cd4_ENCSR277XYX_links.txt"      "/home/merve/Downloads/h3k27me3_ms_cd4_ENCSR526TNC_links.txt"      "/home/merve/Downloads/h3k27me3_ms_cd4_ENCSR592EKF_links.txt"      "/home/merve/Downloads/h3k27me3_ms_cd4_ENCSR613UFD_links.txt"      "/home/merve/Downloads/h3k27me3_ms_cd4_ENCSR740SDR_links.txt"      "/home/merve/Downloads/h3k27me3_ms_cd4_ENCSR779JLY_links.txt"      "/home/merve/Downloads/h3k27me3_ms_cd4_ENCSR993CTA_links.txt"      "/home/merve/Downloads/h3k27me3_ms_cd8_ENCSR116FVG_links.txt"      "/home/merve/Downloads/h3k27me3_ms_cd8_ENCSR122JCM_links.txt"      "/home/merve/Downloads/h3k27me3_ms_cd8_ENCSR216ZVA_links.txt"      "/home/merve/Downloads/h3k27me3_ms_cd8_ENCSR284IKS_links.txt"      "/home/merve/Downloads/h3k27me3_ms_cd8_ENCSR385BOZ_links.txt"      "/home/merve/Downloads/h3k27me3_ms_cd8_ENCSR521SFR_links.txt"      "/home/merve/Downloads/h3k36me3_ms_cd4_ENCSR276NGH_links.txt"      "/home/merve/Downloads/h3k36me3_ms_cd4_ENCSR330CQU_links.txt"      "/home/merve/Downloads/h3k36me3_ms_cd4_ENCSR482VIB_links.txt"      "/home/merve/Downloads/h3k36me3_ms_cd4_ENCSR785PKM_links.txt"      "/home/merve/Downloads/h3k36me3_ms_cd4_ENCSR865FTW_links.txt"      "/home/merve/Downloads/h3k36me3_ms_cd4_ENCSR898VJE_links.txt"      "/home/merve/Downloads/h3k36me3_ms_cd8_ENCSR239OWL_links.txt"      "/home/merve/Downloads/h3k36me3_ms_cd8_ENCSR435MSB_links.txt"      "/home/merve/Downloads/h3k36me3_ms_cd8_ENCSR631VIW_links.txt"      "/home/merve/Downloads/h3k36me3_ms_cd8_ENCSR652WFO_links.txt"      "/home/merve/Downloads/h3k36me3_ms_cd8_ENCSR757FGN_links.txt"      "/home/merve/Downloads/h3k4me1_normal_b_ENCSR156BXM_links.txt"      "/home/merve/Downloads/h3k4me3_normal_b_ENCSR791CAF_links.txt"      "/home/merve/Downloads/h3k9me3_normal_b_ENCSR445LTM_links.txt"      "/home/merve/Downloads/h3k27ac_normal_b_ENCSR685KZA_links.txt"      "/home/merve/Downloads/h3k4me1_normal_nk_ENCSR277YKG_links.txt"      "/home/merve/Downloads/h3k4me3_normal_nk_ENCSR394JFQ_links.txt"      "/home/merve/Downloads/h3k9me3_normal_nk_ENCSR025UNZ_links.txt"      "/home/merve/Downloads/h3k27ac_normal_nk_ENCSR977FMZ_links.txt"      "/home/merve/Downloads/h3k27me3_normal_b_ENCSR589LHR_links.txt"      "/home/merve/Downloads/h3k36me3_normal_b_ENCSR831AXK_links.txt"      "/home/merve/Downloads/h3k4me1_normal_cd4_ENCSR102SOR_links.txt"      "/home/merve/Downloads/h3k4me1_normal_cd8_ENCSR217SHH_links.txt"      "/home/merve/Downloads/h3k4me3_normal_cd4_ENCSR537KJA_links.txt"      "/home/merve/Downloads/h3k4me3_normal_cd8_ENCSR123ZAT_links.txt"      "/home/merve/Downloads/h3k9me3_normal_cd4_ENCSR433EWI_links.txt"      "/home/merve/Downloads/h3k9me3_normal_cd8_ENCSR294HTM_links.txt"      "/home/merve/Downloads/h3k27ac_normal_cd4_ENCSR819NCZ_links.txt"      "/home/merve/Downloads/h3k27ac_normal_cd8_ENCSR976RWL_links.txt"      "/home/merve/Downloads/h3k27me3_normal_nk_ENCSR639NIG_links.txt"      "/home/merve/Downloads/h3k36me3_normal_nk_ENCSR056GJY_links.txt"      "/home/merve/Downloads/h3k27me3_normal_cd4_ENCSR068YVZ_links.txt"      "/home/merve/Downloads/h3k27me3_normal_cd8_ENCSR720QRX_links.txt"      "/home/merve/Downloads/h3k36me3_normal_cd4_ENCSR864OKB_links.txt"      "/home/merve/Downloads/h3k36me3_normal_cd8_ENCSR303SQG_links.txt"      "/home/merve/Downloads/h3k27me3_normal_cd4_ENCSR002UMT_links.txt"      "/home/merve/Downloads/h3k4me3_normal_cd4_ENCSR935ELX_links.txt")
output_file="merged_text.txt"
cat "${files_to_merge[@]}" > "/home/merve/Downloads/$output_file"            
echo "Merged files into: $output_file"

""" !/bin/bash
SBATCH --job-name=bwa_single
SBATCH --nodes=1
SBATCH --account=users
SBATCH --qos=cuda
SBATCH --partition=cuda
SBATCH --cpus-per-task=16
SBATCH -o /cta/users/merve.kaftancioglu/alternative_scenario/bwa_as/out-err/%x-%j-slurm.out
SBATCH -e /cta/users/merve.kaftancioglu/alternative_scenario/bwa_as/out-err/%x-%j-slurm.err

 Load BWA module
module load bwa-0.7.17-gcc-9.2.0-xwkclqy

 Define directories
input_dir="/cta/users/merve.kaftancioglu/alternative_scenario/trimmomatic_as/output"
output_dir="/cta/users/merve.kaftancioglu/alternative_scenario/bwa_as/output"
reference_genome="/cta/users/merve.kaftancioglu/alternative_scenario/bwa_as/reference/human_g1k_v37.fasta"

 Check if BWA index exists
if [ ! -f "$reference_genome.bwt" ]; then
  echo "Error: BWA index not found. Please create the index using: bwa mem index $reference_genome"
  exit 1
fi

sample_ids=("h3k4me1_ms_cd8_ENCSR940PHE" "h3k4me3_ms_cd4_ENCSR180NCM" "h3k4me3_ms_cd4_ENCSR341QLC" "h3k4me3_ms_cd4_ENCSR482TGI" "h3k4me3_ms_cd4_ENCSR486XJK" "h3k4me3_ms_cd4_ENCSR496LKR" "h3k4me3_ms_cd4_ENCSR603LTN" "h3k4me3_ms_cd4_ENCSR802MXQ" "h3k4me3_ms_cd4_ENCSR878YHM" "h3k4me3_ms_cd4_ENCSR954ZLD" "h3k4me3_ms_cd8_ENCSR231ZZH" "h3k4me3_ms_cd8_ENCSR278QHR" "h3k4me3_ms_cd8_ENCSR516CKJ" "h3k4me3_ms_cd8_ENCSR535YYH" "h3k4me3_ms_cd8_ENCSR741XAE" "h3k4me3_ms_cd8_ENCSR848XJL" "h3k9me3_ms_cd4_ENCSR057IZD" "h3k9me3_ms_cd4_ENCSR550DPT" "h3k9me3_ms_cd4_ENCSR677OEF" "h3k9me3_ms_cd4_ENCSR729ZXH" "h3k9me3_ms_cd4_ENCSR851RJV" "h3k9me3_ms_cd4_ENCSR919DFZ" "h3k9me3_ms_cd4_ENCSR953STZ" "h3k9me3_ms_cd4_ENCSR959VZU" "h3k9me3_ms_cd8_ENCSR101USF" "h3k9me3_ms_cd8_ENCSR354GNT" "h3k9me3_ms_cd8_ENCSR377QYB" "h3k9me3_ms_cd8_ENCSR733LCG" "h3k9me3_ms_cd8_ENCSR980KTW" "h3k27ac_ms_cd4_ENCSR200SSJ")

 Loop through all sample subdirectories
for sample_dir in "$input_dir"/*; do
   Extract sample ID from directory name (assuming it matches)
  sample_id=$(basename "$sample_dir")

   Check if it's a directory (skip non-directories)
  if [[ -d "$sample_dir" ]]; then
    echo "Processing sample: $sample_id"

     Create output directory for the sample (if it doesn't exist)
    mkdir -p "$output_dir/$sample_id"

     Find FASTQ files within the sample directory
     Adjust the pattern below if your file names differ
    fastq_files=( "$sample_dir"/*.fastq.gz )

     Check if any FASTQ files found
    if [[ ${fastq_files[@]} -eq 0 ]]; then
      echo "Warning: No FASTQ files found in $sample_dir"
      continue  Skip to next iteration if no FASTQ files
    fi

     BWA mem command with paired-end or single-end handling
    if [[ ${fastq_files[@]} -eq 2 ]]; then
       Paired-end reads (assuming *_1.fastq.gz and *_2.fastq.gz naming)
      bwa mem $reference_genome ${fastq_files[0]} ${fastq_files[1]} > "$output_dir/$sample_id.sam"
    else
       Single-end read (assuming only one FASTQ file)
       Extract only the sample name from the directory path
      sample_name=$(basename "$sample_dir" /)
      bwa mem $reference_genome ${fastq_files[0]} > "$output_dir/$sample_name.sam"
    fi
  fi
done

""" !/bin/bash
SBATCH --job-name=bwa_single
SBATCH --nodes=1
SBATCH --account=users
SBATCH --qos=cuda
SBATCH --partition=cuda
SBATCH --cpus-per-task=16
SBATCH -o /cta/users/merve.kaftancioglu/alternative_scenario/bwa_as/out-err/%x-%j-slurm.out
SBATCH -e /cta/users/merve.kaftancioglu/alternative_scenario/bwa_as/out-err/%x-%j-slurm.err

 Load BWA module
module load bwa-0.7.17-gcc-9.2.0-xwkclqy

 Define directories
input_dir="/cta/users/merve.kaftancioglu/alternative_scenario/trimmomatic_as/output"
output_dir="/cta/users/merve.kaftancioglu/alternative_scenario/bwa_as/output"
reference_genome="/cta/users/merve.kaftancioglu/alternative_scenario/bwa_as/reference/human_g1k_v37.fasta"

 Check if BWA index exists
if [ ! -f "$reference_genome.bwt" ]; then
  echo "Error: BWA index not found. Please create the index using: bwa mem index $reference_genome"
  exit 1
fi

sample_ids=("h3k27ac_ms_cd4_ENCSR322MTA" "h3k27ac_ms_cd4_ENCSR331WMS" "h3k27ac_ms_cd4_ENCSR350UKV" "h3k27ac_ms_cd4_ENCSR474PYR" "h3k27ac_ms_cd4_ENCSR520QDR" "h3k27ac_ms_cd4_ENCSR540XNK" "h3k27ac_ms_cd4_ENCSR705VSO" "h3k27ac_ms_cd4_ENCSR832UMM" "h3k27ac_ms_cd8_ENCSR078ATS" "h3k27ac_ms_cd8_ENCSR348YRH" "h3k27ac_ms_cd8_ENCSR458TOW" "h3k27ac_ms_cd8_ENCSR476IPR" "h3k27ac_ms_cd8_ENCSR787HDF" "h3k27ac_ms_cd8_ENCSR923JIB" "h3k27me3_ms_nk_ENCSR469QVG" "h3k27me3_ms_nk_ENCSR565WDW" "h3k36me3_ms_nk_ENCSR158VSE" "h3k36me3_ms_nk_ENCSR245KON" "h3k36me3_ms_nk_ENCSR530YDY" "h3k36me_ms_cd4_ENCSR532PXR" "h3k27me3_ms_cd4_ENCSR277XYX" "h3k27me3_ms_cd4_ENCSR526TNC" "h3k27me3_ms_cd4_ENCSR592EKF" "h3k27me3_ms_cd4_ENCSR613UFD" "h3k27me3_ms_cd4_ENCSR740SDR" "h3k27me3_ms_cd4_ENCSR779JLY" "h3k27me3_ms_cd4_ENCSR993CTA" "h3k27me3_ms_cd8_ENCSR116FVG" "h3k27me3_ms_cd8_ENCSR122JCM" "h3k27me3_ms_cd8_ENCSR216ZVA")

 Loop through all sample subdirectories
for sample_dir in "$input_dir"/*; do
   Extract sample ID from directory name (assuming it matches)
  sample_id=$(basename "$sample_dir")

   Check if it's a directory (skip non-directories)
  if [[ -d "$sample_dir" ]]; then
    echo "Processing sample: $sample_id"

     Create output directory for the sample (if it doesn't exist)
    mkdir -p "$output_dir/$sample_id"

     Find FASTQ files within the sample directory
     Adjust the pattern below if your file names differ
    fastq_files=( "$sample_dir"/*.fastq.gz )

     Check if any FASTQ files found
    if [[ ${fastq_files[@]} -eq 0 ]]; then
      echo "Warning: No FASTQ files found in $sample_dir"
      continue  Skip to next iteration if no FASTQ files
    fi

     BWA mem command with paired-end or single-end handling
    if [[ ${fastq_files[@]} -eq 2 ]]; then
       Paired-end reads (assuming *_1.fastq.gz and *_2.fastq.gz naming)
      bwa mem $reference_genome ${fastq_files[0]} ${fastq_files[1]} > "$output_dir/$sample_id.sam"
    else
       Single-end read (assuming only one FASTQ file)
       Extract only the sample name from the directory path
      sample_name=$(basename "$sample_dir" /)
      bwa mem $reference_genome ${fastq_files[0]} > "$output_dir/$sample_name.sam"
    fi
  fi
done

""" !/bin/bash
SBATCH --job-name=bwa_single
SBATCH --nodes=1
SBATCH --account=users
SBATCH --qos=cuda
SBATCH --partition=cuda
SBATCH --cpus-per-task=16
SBATCH -o /cta/users/merve.kaftancioglu/alternative_scenario/bwa_as/out-err/%x-%j-slurm.out
SBATCH -e /cta/users/merve.kaftancioglu/alternative_scenario/bwa_as/out-err/%x-%j-slurm.err

 Load BWA module
module load bwa-0.7.17-gcc-9.2.0-xwkclqy

 Define directories
input_dir="/cta/users/merve.kaftancioglu/alternative_scenario/trimmomatic_as/output"
output_dir="/cta/users/merve.kaftancioglu/alternative_scenario/bwa_as/output"
reference_genome="/cta/users/merve.kaftancioglu/alternative_scenario/bwa_as/reference/human_g1k_v37.fasta"

 Check if BWA index exists
if [ ! -f "$reference_genome.bwt" ]; then
  echo "Error: BWA index not found. Please create the index using: bwa mem index $reference_genome"
  exit 1
fi

sample_ids=("h3k27me3_ms_cd8_ENCSR284IKS" "h3k27me3_ms_cd8_ENCSR385BOZ" "h3k27me3_ms_cd8_ENCSR521SFR" "h3k36me3_ms_cd4_ENCSR276NGH" "h3k36me3_ms_cd4_ENCSR330CQU" "h3k36me3_ms_cd4_ENCSR482VIB" "h3k36me3_ms_cd4_ENCSR785PKM" "h3k36me3_ms_cd4_ENCSR865FTW" "h3k36me3_ms_cd4_ENCSR898VJE" "h3k36me3_ms_cd8_ENCSR239OWL" "h3k36me3_ms_cd8_ENCSR435MSB" "h3k36me3_ms_cd8_ENCSR631VIW" "h3k36me3_ms_cd8_ENCSR652WFO" "h3k36me3_ms_cd8_ENCSR757FGN" "h3k4me1_normal_b_ENCSR156BXM" "h3k4me3_normal_b_ENCSR791CAF" "h3k9me3_normal_b_ENCSR445LTM" "h3k27ac_normal_b_ENCSR685KZA" "h3k4me1_normal_nk_ENCSR277YKG" "h3k4me3_normal_nk_ENCSR394JFQ" "h3k9me3_normal_nk_ENCSR025UNZ" "h3k27ac_normal_nk_ENCSR977FMZ" "h3k27me3_normal_b_ENCSR589LHR" "h3k36me3_normal_b_ENCSR831AXK")

 Loop through all sample subdirectories
for sample_dir in "$input_dir"/*; do
   Extract sample ID from directory name (assuming it matches)
  sample_id=$(basename "$sample_dir")

   Check if it's a directory (skip non-directories)
  if [[ -d "$sample_dir" ]]; then
    echo "Processing sample: $sample_id"

     Create output directory for the sample (if it doesn't exist)
    mkdir -p "$output_dir/$sample_id"

     Find FASTQ files within the sample directory
     Adjust the pattern below if your file names differ
    fastq_files=( "$sample_dir"/*.fastq.gz )

     Check if any FASTQ files found
    if [[ ${fastq_files[@]} -eq 0 ]]; then
      echo "Warning: No FASTQ files found in $sample_dir"
      continue  Skip to next iteration if no FASTQ files
    fi

     BWA mem command with paired-end or single-end handling
    if [[ ${fastq_files[@]} -eq 2 ]]; then
       Paired-end reads (assuming *_1.fastq.gz and *_2.fastq.gz naming)
      bwa mem $reference_genome ${fastq_files[0]} ${fastq_files[1]} > "$output_dir/$sample_id.sam"
    else
       Single-end read (assuming only one FASTQ file)
       Extract only the sample name from the directory path
      sample_name=$(basename "$sample_dir" /)
      bwa mem $reference_genome ${fastq_files[0]} > "$output_dir/$sample_name.sam"
    fi
  fi
done

""" !/bin/bash
SBATCH --job-name=bwa_single
SBATCH --nodes=1
SBATCH --account=users
SBATCH --qos=cuda
SBATCH --partition=cuda
SBATCH --cpus-per-task=16
SBATCH -o /cta/users/merve.kaftancioglu/alternative_scenario/bwa_as/out-err/%x-%j-slurm.out
SBATCH -e /cta/users/merve.kaftancioglu/alternative_scenario/bwa_as/out-err/%x-%j-slurm.err

 Load BWA module
module load bwa-0.7.17-gcc-9.2.0-xwkclqy

 Define directories
input_dir="/cta/users/merve.kaftancioglu/alternative_scenario/trimmomatic_as/output"
output_dir="/cta/users/merve.kaftancioglu/alternative_scenario/bwa_as/output"
reference_genome="/cta/users/merve.kaftancioglu/alternative_scenario/bwa_as/reference/human_g1k_v37.fasta"

 Check if BWA index exists
if [ ! -f "$reference_genome.bwt" ]; then
  echo "Error: BWA index not found. Please create the index using: bwa mem index $reference_genome"
  exit 1
fi

sample_ids=("h3k27me3_normal_cd4_ENCSR002UMT" "h3k4me3_normal_cd4_ENCSR935ELX" "h3k4me1_normal_cd4_ENCSR102SOR" "h3k4me1_normal_cd8_ENCSR217SHH" "h3k4me3_normal_cd4_ENCSR537KJA" "h3k4me3_normal_cd8_ENCSR123ZAT" "h3k9me3_normal_cd4_ENCSR433EWI" "h3k9me3_normal_cd8_ENCSR294HTM" "h3k27ac_normal_cd4_ENCSR819NCZ" "h3k27ac_normal_cd8_ENCSR976RWL" "h3k27me3_normal_nk_ENCSR639NIG" "h3k36me3_normal_nk_ENCSR056GJY" "h3k27me3_normal_cd4_ENCSR068YVZ" "h3k27me3_normal_cd8_ENCSR720QRX" "h3k36me3_normal_cd4_ENCSR864OKB" "h3k36me3_normal_cd8_ENCSR303SQG")

 Loop through all sample subdirectories
for sample_dir in "$input_dir"/*; do
   Extract sample ID from directory name (assuming it matches)
  sample_id=$(basename "$sample_dir")

   Check if it's a directory (skip non-directories)
  if [[ -d "$sample_dir" ]]; then
    echo "Processing sample: $sample_id"

     Create output directory for the sample (if it doesn't exist)
    mkdir -p "$output_dir/$sample_id"

     Find FASTQ files within the sample directory
     Adjust the pattern below if your file names differ
    fastq_files=( "$sample_dir"/*.fastq.gz )

     Check if any FASTQ files found
    if [[ ${fastq_files[@]} -eq 0 ]]; then
      echo "Warning: No FASTQ files found in $sample_dir"
      continue  Skip to next iteration if no FASTQ files
    fi

     BWA mem command with paired-end or single-end handling
    if [[ ${fastq_files[@]} -eq 2 ]]; then
       Paired-end reads (assuming *_1.fastq.gz and *_2.fastq.gz naming)
      bwa mem $reference_genome ${fastq_files[0]} ${fastq_files[1]} > "$output_dir/$sample_id.sam"
    else
       Single-end read (assuming only one FASTQ file)
       Extract only the sample name from the directory path
      sample_name=$(basename "$sample_dir" /)
      bwa mem $reference_genome ${fastq_files[0]} > "$output_dir/$sample_name.sam"
    fi
  fi
done

""" !/bin/bash
SBATCH --job-name=bwa_single
SBATCH --nodes=1
SBATCH --account=users
SBATCH --qos=cuda
SBATCH --partition=cuda
SBATCH --cpus-per-task=16
SBATCH -o /cta/users/merve.kaftancioglu/alternative_scenario/bwa_as/out-err/%x-%j-slurm.out
SBATCH -e /cta/users/merve.kaftancioglu/alternative_scenario/bwa_as/out-err/%x-%j-slurm.err

 Load BWA module
module load bwa-0.7.17-gcc-9.2.0-xwkclqy

 Define directories (consider adjusting paths)
input_dir="/cta/users/merve.kaftancioglu/alternative_scenario/trimmomatic_as/output"
output_dir="/cta/users/merve.kaftancioglu/alternative_scenario/bwa_as/output"
reference_genome="/cta/users/merve.kaftancioglu/alternative_scenario/bwa_as/reference/human_g1k_v37.fasta"

 Check if BWA index exists
if [ ! -f "$reference_genome.bwt" ]; then
  echo "Error: BWA index not found. Please create the index using: bwa mem index $reference_genome"
  exit 1
fi

 Define sample IDs (modify as needed)
sample_ids=("h3k27ac_ms_b_ENCSR295QZX" "h3k27ac_ms_b_ENCSR538URI" "h3k27ac_ms_b_ENCSR617GMR" "h3k4me1_ms_b_ENCSR084FXT" "h3k4me1_ms_b_ENCSR480ITK" "h3k4me1_ms_b_ENCSR759DHN" "h3k4me1_ms_b_ENCSR900IAM" "h3k4me1_ms_b_ENCSR931JYC" "h3k4me3_ms_b_ENCSR260CRI" "h3k4me3_ms_b_ENCSR461QMZ" "h3k4me3_ms_b_ENCSR530AVK" "h3k4me3_ms_b_ENCSR713ZYF" "h3k4me3_ms_b_ENCSR923EGO" "h3k9me3_ms_b_ENCSR445DNH" "h3k9me3_ms_b_ENCSR486SZN" "h3k9me3_ms_b_ENCSR983MKX" "h3k27ac_ms_b_ENCSR364OIK" "h3k27ac_ms_b_ENCSR541PLO" "h3k36me3_ms_cd4_ENCSR471WZE" "h3k36me3_ms_cd4_ENCSR473TGK" "h3k4me1_ms_nk_ENCSR482UGV" "h3k4me1_ms_nk_ENCSR718LCI" "h3k4me1_ms_nk_ENCSR806HUT" "h3k4me3_ms_nk_ENCSR234BAD" "h3k4me3_ms_nk_ENCSR703TWY" "h3k4me3_ms_nk_ENCSR755ZGS" "h3k9me3_ms_nk_ENCSR061ATV" "h3k9me3_ms_nk_ENCSR611YIA" "h3k9me3_ms_nk_ENCSR667HUI" "h3k27ac_ms_nk_ENCSR246TTM")

 Loop through all sample subdirectories
for sample_dir in "$input_dir"/*; do
   Check if directory name matches any sample ID
  if [[ " ${sample_ids[@]} " =~ " $(basename "$sample_dir") " ]]; then
    echo "Processing sample: $(basename "$sample_dir")"

     Create output directory for the sample (if it doesn't exist)
    mkdir -p "$output_dir/$(basename "$sample_dir")"

     Find FASTQ files within the sample directory
     Adjust the pattern below if your file names differ
    fastq_files=( "$sample_dir"/*.fastq.gz )

     Check if any FASTQ files found
    if [[ ${fastq_files[@]} -eq 0 ]]; then
      echo "Warning: No FASTQ files found in $sample_dir"
      continue  Skip to next iteration if no FASTQ files
    fi

     BWA mem command with paired-end or single-end handling
    if [[ ${fastq_files[@]} -eq 2 ]]; then
       Paired-end reads (assuming *_1.fastq.gz and *_2.fastq.gz naming)
      bwa mem $reference_genome ${fastq_files[0]} ${fastq_files[1]} > "$output_dir/$(basename "$sample_dir").sam"
    else
       Single-end read (assuming only one FASTQ file)
       Extract only the sample name from the directory path
      sample_name=$(basename "$sample_dir" /)
      bwa mem $reference_genome ${fastq_files[0]} > "$output_dir/$sample_name.sam"
    fi
  fi
done

""" !/bin/bash
SBATCH --job-name=bwa_single
SBATCH --nodes=1
SBATCH --account=users
SBATCH --qos=cuda
SBATCH --partition=cuda
SBATCH --cpus-per-task=16
SBATCH -o /cta/users/merve.kaftancioglu/alternative_scenario/bwa_as/out-err/%x-%j-slurm.out
SBATCH -e /cta/users/merve.kaftancioglu/alternative_scenario/bwa_as/out-err/%x-%j-slurm.err

 Load BWA module
module load bwa-0.7.17-gcc-9.2.0-xwkclqy

 Define directories (consider adjusting paths)
input_dir="/cta/users/merve.kaftancioglu/alternative_scenario/trimmomatic_as/output"
output_dir="/cta/users/merve.kaftancioglu/alternative_scenario/bwa_as/output"
reference_genome="/cta/users/merve.kaftancioglu/alternative_scenario/bwa_as/reference/human_g1k_v37.fasta"

 Check if BWA index exists
if [ ! -f "$reference_genome.bwt" ]; then
  echo "Error: BWA index not found. Please create the index using: bwa mem index $reference_genome"
  exit 1
fi

 Define sample IDs (modify as needed)
sample_ids=("h3k27ac_ms_nk_ENCSR469CFL" "h3k27ac_ms_nk_ENCSR746AIX" "h3k27me3_ms_b_ENCSR009ZRH" "h3k27me3_ms_b_ENCSR182NLA" "h3k27me3_ms_b_ENCSR272YVX" "h3k27me3_ms_b_ENCSR649FUX" "h3k27me3_ms_b_ENCSR842WWX" "h3k36me3_ms_b_ENCSR089VQL" "h3k36me3_ms_b_ENCSR119PSR" "h3k36me3_ms_b_ENCSR238WFK" "h3k36me3_ms_b_ENCSR987OPY" "h3k4me1_ms_cd4_ENCSR036YVP" "h3k4me1_ms_cd4_ENCSR043JIA" "h3k4me1_ms_cd4_ENCSR093VUP" "h3k4me1_ms_cd4_ENCSR124IGC" "h3k4me1_ms_cd4_ENCSR315MYP" "h3k4me1_ms_cd4_ENCSR458QEY" "h3k4me1_ms_cd4_ENCSR485FBT" "h3k4me1_ms_cd4_ENCSR641ZFV" "h3k4me1_ms_cd4_ENCSR687CQX" "h3k4me1_ms_cd8_ENCSR231XAP" "h3k4me1_ms_cd8_ENCSR572XTB" "h3k4me1_ms_cd8_ENCSR788TEF" "h3k4me1_ms_cd8_ENCSR815YZL" "h3k4me1_ms_cd8_ENCSR861FEC")

 Loop through all sample subdirectories
for sample_dir in "$input_dir"/*; do
   Check if directory name matches any sample ID
  if [[ " ${sample_ids[@]} " =~ " $(basename "$sample_dir") " ]]; then
    echo "Processing sample: $(basename "$sample_dir")"

     Create output directory for the sample (if it doesn't exist)
    mkdir -p "$output_dir/$(basename "$sample_dir")"

     Find FASTQ files within the sample directory
     Adjust the pattern below if your file names differ
    fastq_files=( "$sample_dir"/*.fastq.gz )

     Check if any FASTQ files found
    if [[ ${fastq_files[@]} -eq 0 ]]; then
      echo "Warning: No FASTQ files found in $sample_dir"
      continue  Skip to next iteration if no FASTQ files
    fi

     BWA mem command with paired-end or single-end handling
    if [[ ${fastq_files[@]} -eq 2 ]]; then
       Paired-end reads (assuming *_1.fastq.gz and *_2.fastq.gz naming)
      bwa mem $reference_genome ${fastq_files[0]} ${fastq_files[1]} > "$output_dir/$(basename "$sample_dir").sam"
    else
       Single-end read (assuming only one FASTQ file)
       Extract only the sample name from the directory path
      sample_name=$(basename "$sample_dir" /)
      bwa mem $reference_genome ${fastq_files[0]} > "$output_dir/$sample_name.sam"
    fi
  fi
done

""" !/bin/bash
SBATCH --job-name=bwa_single
SBATCH --nodes=1
SBATCH --account=users
SBATCH --qos=cuda
SBATCH --partition=cuda
SBATCH --cpus-per-task=16
SBATCH -o /cta/users/merve.kaftancioglu/alternative_scenario/bwa_as/out-err/%x-%j-slurm.out
SBATCH -e /cta/users/merve.kaftancioglu/alternative_scenario/bwa_as/out-err/%x-%j-slurm.err

 Load BWA module
module load bwa-0.7.17-gcc-9.2.0-xwkclqy

 Define directories (consider adjusting paths)
input_dir="/cta/users/merve.kaftancioglu/alternative_scenario/trimmomatic_as/output"
output_dir="/cta/users/merve.kaftancioglu/alternative_scenario/bwa_as/output"
reference_genome="/cta/users/merve.kaftancioglu/alternative_scenario/bwa_as/reference/human_g1k_v37.fasta"

 Check if BWA index exists
if [ ! -f "$reference_genome.bwt" ]; then
  echo "Error: BWA index not found. Please create the index using: bwa mem index $reference_genome"
  exit 1
fi

 Define sample IDs (modify as needed)
sample_ids=("h3k4me1_ms_cd8_ENCSR940PHE" "h3k4me3_ms_cd4_ENCSR180NCM" "h3k4me3_ms_cd4_ENCSR341QLC" "h3k4me3_ms_cd4_ENCSR482TGI" "h3k4me3_ms_cd4_ENCSR486XJK" "h3k4me3_ms_cd4_ENCSR496LKR" "h3k4me3_ms_cd4_ENCSR603LTN" "h3k4me3_ms_cd4_ENCSR802MXQ" "h3k4me3_ms_cd4_ENCSR878YHM" "h3k4me3_ms_cd4_ENCSR954ZLD" "h3k4me3_ms_cd8_ENCSR231ZZH" "h3k4me3_ms_cd8_ENCSR278QHR" "h3k4me3_ms_cd8_ENCSR516CKJ" "h3k4me3_ms_cd8_ENCSR535YYH" "h3k4me3_ms_cd8_ENCSR741XAE" "h3k4me3_ms_cd8_ENCSR848XJL" "h3k9me3_ms_cd4_ENCSR057IZD" "h3k9me3_ms_cd4_ENCSR550DPT" "h3k9me3_ms_cd4_ENCSR677OEF" "h3k9me3_ms_cd4_ENCSR729ZXH" "h3k9me3_ms_cd4_ENCSR851RJV" "h3k9me3_ms_cd4_ENCSR919DFZ" "h3k9me3_ms_cd4_ENCSR953STZ" "h3k9me3_ms_cd4_ENCSR959VZU" "h3k9me3_ms_cd8_ENCSR101USF" "h3k9me3_ms_cd8_ENCSR354GNT" "h3k9me3_ms_cd8_ENCSR377QYB" "h3k9me3_ms_cd8_ENCSR733LCG" "h3k9me3_ms_cd8_ENCSR980KTW" "h3k27ac_ms_cd4_ENCSR200SSJ")

 Loop through all sample subdirectories
for sample_dir in "$input_dir"/*; do
   Check if directory name matches any sample ID
  if [[ " ${sample_ids[@]} " =~ " $(basename "$sample_dir") " ]]; then
    echo "Processing sample: $(basename "$sample_dir")"

     Create output directory for the sample (if it doesn't exist)
    mkdir -p "$output_dir/$(basename "$sample_dir")"

     Find FASTQ files within the sample directory
     Adjust the pattern below if your file names differ
    fastq_files=( "$sample_dir"/*.fastq.gz )

     Check if any FASTQ files found
    if [[ ${fastq_files[@]} -eq 0 ]]; then
      echo "Warning: No FASTQ files found in $sample_dir"
      continue  Skip to next iteration if no FASTQ files
    fi

     BWA mem command with paired-end or single-end handling
    if [[ ${fastq_files[@]} -eq 2 ]]; then
       Paired-end reads (assuming *_1.fastq.gz and *_2.fastq.gz naming)
      bwa mem $reference_genome ${fastq_files[0]} ${fastq_files[1]} > "$output_dir/$(basename "$sample_dir").sam"
    else
       Single-end read (assuming only one FASTQ file)
       Extract only the sample name from the directory path
      sample_name=$(basename "$sample_dir" /)
      bwa mem $reference_genome ${fastq_files[0]} > "$output_dir/$sample_name.sam"
    fi
  fi
done

""" !/bin/bash
SBATCH --job-name=bwa_single
SBATCH --nodes=1
SBATCH --account=users
SBATCH --qos=cuda
SBATCH --partition=cuda
SBATCH --cpus-per-task=16
SBATCH -o /cta/users/merve.kaftancioglu/alternative_scenario/bwa_as/out-err/%x-%j-slurm.out
SBATCH -e /cta/users/merve.kaftancioglu/alternative_scenario/bwa_as/out-err/%x-%j-slurm.err

 Load BWA module
module load bwa-0.7.17-gcc-9.2.0-xwkclqy

 Define directories (consider adjusting paths)
input_dir="/cta/users/merve.kaftancioglu/alternative_scenario/trimmomatic_as/output"
output_dir="/cta/users/merve.kaftancioglu/alternative_scenario/bwa_as/output"
reference_genome="/cta/users/merve.kaftancioglu/alternative_scenario/bwa_as/reference/human_g1k_v37.fasta"

 Check if BWA index exists
if [ ! -f "$reference_genome.bwt" ]; then
  echo "Error: BWA index not found. Please create the index using: bwa mem index $reference_genome"
  exit 1
fi

 Define sample IDs (modify as needed)
sample_ids=("h3k27ac_ms_cd4_ENCSR322MTA" "h3k27ac_ms_cd4_ENCSR331WMS" "h3k27ac_ms_cd4_ENCSR350UKV" "h3k27ac_ms_cd4_ENCSR474PYR" "h3k27ac_ms_cd4_ENCSR520QDR" "h3k27ac_ms_cd4_ENCSR540XNK" "h3k27ac_ms_cd4_ENCSR705VSO" "h3k27ac_ms_cd4_ENCSR832UMM" "h3k27ac_ms_cd8_ENCSR078ATS" "h3k27ac_ms_cd8_ENCSR348YRH" "h3k27ac_ms_cd8_ENCSR458TOW" "h3k27ac_ms_cd8_ENCSR476IPR" "h3k27ac_ms_cd8_ENCSR787HDF" "h3k27ac_ms_cd8_ENCSR923JIB" "h3k27me3_ms_nk_ENCSR469QVG" "h3k27me3_ms_nk_ENCSR565WDW" "h3k36me3_ms_nk_ENCSR158VSE" "h3k36me3_ms_nk_ENCSR245KON" "h3k36me3_ms_nk_ENCSR530YDY" "h3k36me_ms_cd4_ENCSR532PXR" "h3k27me3_ms_cd4_ENCSR277XYX" "h3k27me3_ms_cd4_ENCSR526TNC" "h3k27me3_ms_cd4_ENCSR592EKF" "h3k27me3_ms_cd4_ENCSR613UFD" "h3k27me3_ms_cd4_ENCSR740SDR" "h3k27me3_ms_cd4_ENCSR779JLY" "h3k27me3_ms_cd4_ENCSR993CTA" "h3k27me3_ms_cd8_ENCSR116FVG" "h3k27me3_ms_cd8_ENCSR122JCM" "h3k27me3_ms_cd8_ENCSR216ZVA")

 Loop through all sample subdirectories
for sample_dir in "$input_dir"/*; do
   Check if directory name matches any sample ID
  if [[ " ${sample_ids[@]} " =~ " $(basename "$sample_dir") " ]]; then
    echo "Processing sample: $(basename "$sample_dir")"

     Create output directory for the sample (if it doesn't exist)
    mkdir -p "$output_dir/$(basename "$sample_dir")"

     Find FASTQ files within the sample directory
     Adjust the pattern below if your file names differ
    fastq_files=( "$sample_dir"/*.fastq.gz )

     Check if any FASTQ files found
    if [[ ${fastq_files[@]} -eq 0 ]]; then
      echo "Warning: No FASTQ files found in $sample_dir"
      continue  Skip to next iteration if no FASTQ files
    fi

     BWA mem command with paired-end or single-end handling
    if [[ ${fastq_files[@]} -eq 2 ]]; then
       Paired-end reads (assuming *_1.fastq.gz and *_2.fastq.gz naming)
      bwa mem $reference_genome ${fastq_files[0]} ${fastq_files[1]} > "$output_dir/$(basename "$sample_dir").sam"
    else
       Single-end read (assuming only one FASTQ file)
       Extract only the sample name from the directory path
      sample_name=$(basename "$sample_dir" /)
      bwa mem $reference_genome ${fastq_files[0]} > "$output_dir/$sample_name.sam"
    fi
  fi
done

""" !/bin/bash
SBATCH --job-name=bwa_single
SBATCH --nodes=1
SBATCH --account=users
SBATCH --qos=cuda
SBATCH --partition=cuda
SBATCH --cpus-per-task=16
SBATCH -o /cta/users/merve.kaftancioglu/alternative_scenario/bwa_as/out-err/%x-%j-slurm.out
SBATCH -e /cta/users/merve.kaftancioglu/alternative_scenario/bwa_as/out-err/%x-%j-slurm.err

 Load BWA module
module load bwa-0.7.17-gcc-9.2.0-xwkclqy

 Define directories (consider adjusting paths)
input_dir="/cta/users/merve.kaftancioglu/alternative_scenario/trimmomatic_as/output"
output_dir="/cta/users/merve.kaftancioglu/alternative_scenario/bwa_as/output"
reference_genome="/cta/users/merve.kaftancioglu/alternative_scenario/bwa_as/reference/human_g1k_v37.fasta"

 Check if BWA index exists
if [ ! -f "$reference_genome.bwt" ]; then
  echo "Error: BWA index not found. Please create the index using: bwa mem index $reference_genome"
  exit 1
fi

 Define sample IDs (modify as needed)
sample_ids=("h3k27me3_ms_cd8_ENCSR284IKS" "h3k27me3_ms_cd8_ENCSR385BOZ" "h3k27me3_ms_cd8_ENCSR521SFR" "h3k36me3_ms_cd4_ENCSR276NGH" "h3k36me3_ms_cd4_ENCSR330CQU" "h3k36me3_ms_cd4_ENCSR482VIB" "h3k36me3_ms_cd4_ENCSR785PKM" "h3k36me3_ms_cd4_ENCSR865FTW" "h3k36me3_ms_cd4_ENCSR898VJE" "h3k36me3_ms_cd8_ENCSR239OWL" "h3k36me3_ms_cd8_ENCSR435MSB" "h3k36me3_ms_cd8_ENCSR631VIW" "h3k36me3_ms_cd8_ENCSR652WFO" "h3k36me3_ms_cd8_ENCSR757FGN" "h3k4me1_normal_b_ENCSR156BXM" "h3k4me3_normal_b_ENCSR791CAF" "h3k9me3_normal_b_ENCSR445LTM" "h3k27ac_normal_b_ENCSR685KZA" "h3k4me1_normal_nk_ENCSR277YKG" "h3k4me3_normal_nk_ENCSR394JFQ" "h3k9me3_normal_nk_ENCSR025UNZ" "h3k27ac_normal_nk_ENCSR977FMZ" "h3k27me3_normal_b_ENCSR589LHR" "h3k36me3_normal_b_ENCSR831AXK")

 Loop through all sample subdirectories
for sample_dir in "$input_dir"/*; do
   Check if directory name matches any sample ID
  if [[ " ${sample_ids[@]} " =~ " $(basename "$sample_dir") " ]]; then
    echo "Processing sample: $(basename "$sample_dir")"

     Create output directory for the sample (if it doesn't exist)
    mkdir -p "$output_dir/$(basename "$sample_dir")"

     Find FASTQ files within the sample directory
     Adjust the pattern below if your file names differ
    fastq_files=( "$sample_dir"/*.fastq.gz )

     Check if any FASTQ files found
    if [[ ${fastq_files[@]} -eq 0 ]]; then
      echo "Warning: No FASTQ files found in $sample_dir"
      continue  Skip to next iteration if no FASTQ files
    fi

     BWA mem command with paired-end or single-end handling
    if [[ ${fastq_files[@]} -eq 2 ]]; then
       Paired-end reads (assuming *_1.fastq.gz and *_2.fastq.gz naming)
      bwa mem $reference_genome ${fastq_files[0]} ${fastq_files[1]} > "$output_dir/$(basename "$sample_dir").sam"
    else
       Single-end read (assuming only one FASTQ file)
       Extract only the sample name from the directory path
      sample_name=$(basename "$sample_dir" /)
      bwa mem $reference_genome ${fastq_files[0]} > "$output_dir/$sample_name.sam"
    fi
  fi
done

""" !/bin/bash
SBATCH --job-name=bwa_single
SBATCH --nodes=1
SBATCH --account=users
SBATCH --qos=cuda
SBATCH --partition=cuda
SBATCH --cpus-per-task=16
SBATCH -o /cta/users/merve.kaftancioglu/alternative_scenario/bwa_as/out-err/%x-%j-slurm.out
SBATCH -e /cta/users/merve.kaftancioglu/alternative_scenario/bwa_as/out-err/%x-%j-slurm.err

 Load BWA module
module load bwa-0.7.17-gcc-9.2.0-xwkclqy

 Define directories (consider adjusting paths)
input_dir="/cta/users/merve.kaftancioglu/alternative_scenario/trimmomatic_as/output"
output_dir="/cta/users/merve.kaftancioglu/alternative_scenario/bwa_as/output"
reference_genome="/cta/users/merve.kaftancioglu/alternative_scenario/bwa_as/reference/human_g1k_v37.fasta"

 Check if BWA index exists
if [ ! -f "$reference_genome.bwt" ]; then
  echo "Error: BWA index not found. Please create the index using: bwa mem index $reference_genome"
  exit 1
fi

 Define sample IDs (modify as needed)
sample_ids=("h3k27me3_normal_cd4_ENCSR002UMT" "h3k4me3_normal_cd4_ENCSR935ELX" "h3k4me1_normal_cd4_ENCSR102SOR" "h3k4me1_normal_cd8_ENCSR217SHH" "h3k4me3_normal_cd4_ENCSR537KJA" "h3k4me3_normal_cd8_ENCSR123ZAT" "h3k9me3_normal_cd4_ENCSR433EWI" "h3k9me3_normal_cd8_ENCSR294HTM" "h3k27ac_normal_cd4_ENCSR819NCZ" "h3k27ac_normal_cd8_ENCSR976RWL" "h3k27me3_normal_nk_ENCSR639NIG" "h3k36me3_normal_nk_ENCSR056GJY" "h3k27me3_normal_cd4_ENCSR068YVZ" "h3k27me3_normal_cd8_ENCSR720QRX" "h3k36me3_normal_cd4_ENCSR864OKB" "h3k36me3_normal_cd8_ENCSR303SQG")

 Loop through all sample subdirectories
for sample_dir in "$input_dir"/*; do
   Check if directory name matches any sample ID
  if [[ " ${sample_ids[@]} " =~ " $(basename "$sample_dir") " ]]; then
    echo "Processing sample: $(basename "$sample_dir")"

     Create output directory for the sample (if it doesn't exist)
    mkdir -p "$output_dir/$(basename "$sample_dir")"

     Find FASTQ files within the sample directory
     Adjust the pattern below if your file names differ
    fastq_files=( "$sample_dir"/*.fastq.gz )

     Check if any FASTQ files found
    if [[ ${fastq_files[@]} -eq 0 ]]; then
      echo "Warning: No FASTQ files found in $sample_dir"
      continue  Skip to next iteration if no FASTQ files
    fi

     BWA mem command with paired-end or single-end handling
    if [[ ${fastq_files[@]} -eq 2 ]]; then
       Paired-end reads (assuming *_1.fastq.gz and *_2.fastq.gz naming)
      bwa mem $reference_genome ${fastq_files[0]} ${fastq_files[1]} > "$output_dir/$(basename "$sample_dir").sam"
    else
       Single-end read (assuming only one FASTQ file)
       Extract only the sample name from the directory path
      sample_name=$(basename "$sample_dir" /)
      bwa mem $reference_genome ${fastq_files[0]} > "$output_dir/$sample_name.sam"
    fi
  fi
done

export PATH=$PATH:/home/merve/apps/hisat2-2.1


""" !/bin/bash
SBATCH --job-name=star_single
SBATCH --nodes=1
SBATCH --account=users
SBATCH --qos=mid_mdbf
SBATCH --partition=mid_mdbf
SBATCH --cpus-per-task=16
SBATCH -o /cta/users/merve.kaftancioglu/alternative_scenario/star_as/out-err/%x-%j-slurm.out
SBATCH -e /cta/users/merve.kaftancioglu/alternative_scenario/star_as/out-err/%x-%j-slurm.err

 Load STAR module (if applicable)
module load star-2.7.0e-gcc-9.2.0-vynasg3   Replace with your module name (if used)

 Check for STAR binary (modify path if needed)
star_path=$(which star 2>/dev/null)

if [[ -z "$star_path" ]]; then
  echo "ERROR: STAR binary not found!"
  echo "Please ensure STAR is installed and accessible on the compute nodes."
  exit 1   Exit script with error code
fi

 Define directories (replace with your actual paths)
genome_dir="/cta/users/merve.kaftancioglu/alternative_scenario/star_as/index"
input_dir="/cta/users/merve.kaftancioglu/alternative_scenario/trimmomatic_as/output"
output_dir="/cta/users/merve.kaftancioglu/alternative_scenario/star_as/output"


 List of sample IDs 
sample_ids=("h3k27ac_ms_b_ENCSR295QZX" "h3k27ac_ms_b_ENCSR538URI" "h3k27ac_ms_b_ENCSR617GMR")

 Process each sample folder
for sample_id in "${sample_ids[@]}"; do
   Input folder 
  input_folder="/cta/users/merve.kaftancioglu/alternative_scenario/trimmomatic_as/output/$sample_id"

   Output folder 
  output_folder="/cta/users/merve.kaftancioglu/alternative_scenario/star_as/output/$sample_id"

   Create output folder if it doesn't exist
  mkdir -p "$output_folder"

   Process all .fastq.gz files in the input folder
  for input_file in "$input_folder"/*.fastq.gz; do
     Get sample name from filename
    sample_name=$(basename "$input_file" .fastq.gz)

   Run STAR for the current sample
  STAR --genomeDir "$genome_dir" \
       --runThreadN 16 \   Adjust threads based on your resource allocation
       --readFilesIn "$input_fastq" \   No colon after flag name
       --outFileNamePrefix "$output_dir/$sample_id/" \   Include trailing slash for proper path
       --outSAMtype BAM SortedByCoordinate \   No colon after flag name
       --outSAMunmapped Within \
       --outSAMattributes Standard \
       --alignEndsType EndToEnd
  done
done

""" !/bin/bash
SBATCH --job-name=fastQC
SBATCH --nodes=1
SBATCH --account=users
SBATCH --qos=cuda
SBATCH --partition=cuda
SBATCH --cpus-per-task=16
SBATCH -o /cta/users/merve.kaftancioglu/alternative_scenario/fastQC_as/out-err/%x-%j-slurm.out
SBATCH -e /cta/users/merve.kaftancioglu/alternative_scenario/fastQC_as/out-err/%x-%j-slurm.err


 Load FastQC module
module load fastqc-0.11.7-gcc-9.2.0-nsjzjof

substrate_ids=("h3k27me3_ms_cd8_ENCSR521SFR" "h3k27me3_ms_nk_ENCSR565WDW" "h3k27me3_normal_b_ENCSR589LHR" "h3k27me3_normal_cd4_ENCSR068YVZ" "h3k27me3_normal_cd8_ENCSR720QRX" "h3k36me3_ms_b_ENCSR119PSR" "h3k36me3_ms_b_ENCSR238WFK" "h3k36me3_ms_b_ENCSR987OPY" "h3k36me3_ms_cd4_ENCSR473TGK" "h3k36me3_ms_cd4_ENCSR482VIB" "h3k36me3_ms_cd4_ENCSR865FTW" "h3k36me3_ms_cd4_ENCSR898VJE" "h3k36me3_ms_cd8_ENCSR239OWL")

for substrate_id in "${substrate_ids[@]}"; do
    echo "Predicting ${substrate_id}..."
    input_directory="/cta/users/merve.kaftancioglu/alternative_scenario/seqs/${substrate_id}"
    output_dir="/cta/users/merve.kaftancioglu/alternative_scenario/fastQC_as/output/${substrate_id}"    
    mkdir -p "${output_dir}"
    fastqc --outdir "$output_directory" "$input_directory"/*.fastq.gz -threads 16
    echo "Finished ${substrate_id}"
done


echo "DONE!""" "



""" !/bin/bash
SBATCH --job-name=fastQC
SBATCH --nodes=1
SBATCH --account=users
SBATCH --qos=mid_mdbf
SBATCH --partition=mid_mdbf
SBATCH --cpus-per-task=16
SBATCH -o /cta/users/merve.kaftancioglu/alternative_scenario/post-trim_fastQC_as/out-err/%x-%j-slurm.out
SBATCH -e /cta/users/merve.kaftancioglu/alternative_scenario/post-trim_fastQC_as/out-err/%x-%j-slurm.err

 Load FastQC module
module load fastqc-0.11.7-gcc-9.2.0-nsjzjof

sample_ids=("h3k27ac_ms_b_ENCSR295QZX" "h3k27ac_ms_b_ENCSR364OIK" "h3k27ac_ms_b_ENCSR538URI" "h3k27ac_ms_b_ENCSR541PLO" "h3k27ac_ms_b_ENCSR617GMR")

for sample_id in "${sample_ids[@]}"; do
  echo "Predicting ${sample_id}..."
  input_directory="/cta/users/merve.kaftancioglu/alternative_scenario/trimmomatic_as/output/${sample_id}"
  output_dir="/cta/users/merve.kaftancioglu/alternative_scenario/post-trim_fastQC_as/output/${sample_id}"
  mkdir -p "${output_dir}"
  fastqc --outdir "$output_dir" "$input_directory"/*.fastq.gz -threads 16
  echo "Finished ${sample_id}"
done

echo "DONE!""" "

""" !/bin/bash
SBATCH --job-name=fastQC
SBATCH --nodes=1
SBATCH --account=users
SBATCH --qos=mid_mdbf
SBATCH --partition=mid_mdbf
SBATCH --cpus-per-task=16
SBATCH -o /cta/users/merve.kaftancioglu/alternative_scenario/post-trim_fastQC_as/out-err/%x-%j-slurm.out
SBATCH -e /cta/users/merve.kaftancioglu/alternative_scenario/post-trim_fastQC_as/out-err/%x-%j-slurm.err

 Load FastQC module
module load fastqc-0.11.7-gcc-9.2.0-nsjzjof

sample_ids=("h3k27ac_ms_cd4_ENCSR200SSJ" "h3k27ac_ms_cd4_ENCSR322MTA" "h3k27ac_ms_cd4_ENCSR331WMS" "h3k27ac_ms_cd4_ENCSR350UKV" "h3k27ac_ms_cd4_ENCSR474PYR")

for sample_id in "${sample_ids[@]}"; do
  echo "Predicting ${sample_id}..."
  input_directory="/cta/users/merve.kaftancioglu/alternative_scenario/trimmomatic_as/output/${sample_id}"
  output_dir="/cta/users/merve.kaftancioglu/alternative_scenario/post-trim_fastQC_as/output/${sample_id}"
  mkdir -p "${output_dir}"
  fastqc --outdir "$output_dir" "$input_directory"/*.fastq.gz -threads 16
  echo "Finished ${sample_id}"
done

echo "DONE!""" "

""" !/bin/bash
SBATCH --job-name=fastQC
SBATCH --nodes=1
SBATCH --account=users
SBATCH --qos=mid_mdbf
SBATCH --partition=mid_mdbf
SBATCH --cpus-per-task=16
SBATCH -o /cta/users/merve.kaftancioglu/alternative_scenario/post-trim_fastQC_as/out-err/%x-%j-slurm.out
SBATCH -e /cta/users/merve.kaftancioglu/alternative_scenario/post-trim_fastQC_as/out-err/%x-%j-slurm.err

 Load FastQC module
module load fastqc-0.11.7-gcc-9.2.0-nsjzjof

sample_ids=("h3k27ac_ms_cd4_ENCSR520QDR" "h3k27ac_ms_cd4_ENCSR540XNK" "h3k27ac_ms_cd4_ENCSR705VSO" "h3k27ac_ms_cd4_ENCSR832UMM" "h3k27ac_ms_cd8_ENCSR078ATS")

for sample_id in "${sample_ids[@]}"; do
  echo "Predicting ${sample_id}..."
  input_directory="/cta/users/merve.kaftancioglu/alternative_scenario/trimmomatic_as/output/${sample_id}"
  output_dir="/cta/users/merve.kaftancioglu/alternative_scenario/post-trim_fastQC_as/output/${sample_id}"
  mkdir -p "${output_dir}"
  fastqc --outdir "$output_dir" "$input_directory"/*.fastq.gz -threads 16
  echo "Finished ${sample_id}"
done

echo "DONE!""" "

""" !/bin/bash
SBATCH --job-name=fastQC
SBATCH --nodes=1
SBATCH --account=users
SBATCH --qos=cuda
SBATCH --partition=cuda
SBATCH --cpus-per-task=16
SBATCH -o /cta/users/merve.kaftancioglu/alternative_scenario/fastQC_as/out-err/%x-%j-slurm.out
SBATCH -e /cta/users/merve.kaftancioglu/alternative_scenario/fastQC_as/out-err/%x-%j-slurm.err

 Load FastQC module
module load fastqc-0.11.7-gcc-9.2.0-nsjzjof

sample_ids=("h3k27ac_ms_cd8_ENCSR348YRH" "h3k27ac_ms_cd8_ENCSR458TOW" "h3k27ac_ms_cd8_ENCSR476IPR" "h3k27ac_ms_cd8_ENCSR787HDF" "h3k27ac_ms_cd8_ENCSR923JIB")

for sample_id in "${sample_ids[@]}"; do
  echo "Predicting ${sample_id}..."
  input_directory="/cta/users/merve.kaftancioglu/alternative_scenario/trimmomatic_as/output/${sample_id}"
  output_dir="/cta/users/merve.kaftancioglu/alternative_scenario/post-trim_fastQC_as/output/${sample_id}"
  mkdir -p "${output_dir}"
  fastqc --outdir "$output_dir" "$input_directory"/*.fastq.gz -threads 16
  echo "Finished ${sample_id}"
done

echo "DONE!""" "

""" !/bin/bash
SBATCH --job-name=fastQC
SBATCH --nodes=1
SBATCH --account=users
SBATCH --qos=mid_mdbf
SBATCH --partition=mid_mdbf
SBATCH --cpus-per-task=16
SBATCH -o /cta/users/merve.kaftancioglu/alternative_scenario/fastQC_as/out-err/%x-%j-slurm.out
SBATCH -e /cta/users/merve.kaftancioglu/alternative_scenario/fastQC_as/out-err/%x-%j-slurm.err

 Load FastQC module
module load fastqc-0.11.7-gcc-9.2.0-nsjzjof

sample_ids=("h3k27ac_ms_nk_ENCSR246TTM" "h3k27ac_ms_nk_ENCSR469CFL" "h3k27ac_ms_nk_ENCSR746AIX" "h3k27ac_normal_b_ENCSR685KZA" "h3k27ac_normal_cd4_ENCSR819NCZ")

for sample_id in "${sample_ids[@]}"; do
  echo "Predicting ${sample_id}..."
  input_directory="/cta/users/merve.kaftancioglu/alternative_scenario/trimmomatic_as/output/${sample_id}"
  output_dir="/cta/users/merve.kaftancioglu/alternative_scenario/post-trim_fastQC_as/output/${sample_id}"
  mkdir -p "${output_dir}"
  fastqc --outdir "$output_dir" "$input_directory"/*.fastq.gz -threads 16
  echo "Finished ${sample_id}"
done

echo "DONE!""" "

""" !/bin/bash
SBATCH --job-name=fastQC
SBATCH --nodes=1
SBATCH --account=users
SBATCH --qos=cuda
SBATCH --partition=cuda
SBATCH --cpus-per-task=16
SBATCH -o /cta/users/merve.kaftancioglu/alternative_scenario/fastQC_as/out-err/%x-%j-slurm.out
SBATCH -e /cta/users/merve.kaftancioglu/alternative_scenario/fastQC_as/out-err/%x-%j-slurm.err

 Load FastQC module
module load fastqc-0.11.7-gcc-9.2.0-nsjzjof

sample_ids=("h3k27ac_normal_cd8_ENCSR976RWL" "h3k27ac_normal_nk_ENCSR977FMZ" "h3k27me3_ms_b_ENCSR009ZRH" "h3k27me3_ms_b_ENCSR182NLA" "h3k27me3_ms_b_ENCSR272YVX")

for sample_id in "${sample_ids[@]}"; do
  echo "Predicting ${sample_id}..."
  input_directory="/cta/users/merve.kaftancioglu/alternative_scenario/trimmomatic_as/output/${sample_id}"
  output_dir="/cta/users/merve.kaftancioglu/alternative_scenario/post-trim_fastQC_as/output/${sample_id}"
  mkdir -p "${output_dir}"
  fastqc --outdir "$output_dir" "$input_directory"/*.fastq.gz -threads 16
  echo "Finished ${sample_id}"
done

echo "DONE!""" "

""" !/bin/bash
SBATCH --job-name=fastQC
SBATCH --nodes=1
SBATCH --account=users
SBATCH --qos=mid_mdbf
SBATCH --partition=mid_mdbf
SBATCH --cpus-per-task=16
SBATCH -o /cta/users/merve.kaftancioglu/alternative_scenario/fastQC_as/out-err/%x-%j-slurm.out
SBATCH -e /cta/users/merve.kaftancioglu/alternative_scenario/fastQC_as/out-err/%x-%j-slurm.err

 Load FastQC module
module load fastqc-0.11.7-gcc-9.2.0-nsjzjof

sample_ids=("h3k27me3_ms_b_ENCSR649FUX" "h3k27me3_ms_b_ENCSR842WWX" "h3k27me3_ms_cd4_ENCSR277XYX" "h3k27me3_ms_cd4_ENCSR526TNC" "h3k27me3_ms_cd4_ENCSR592EKF")

for sample_id in "${sample_ids[@]}"; do
  echo "Predicting ${sample_id}..."
  input_directory="/cta/users/merve.kaftancioglu/alternative_scenario/trimmomatic_as/output/${sample_id}"
  output_dir="/cta/users/merve.kaftancioglu/alternative_scenario/post-trim_fastQC_as/output/${sample_id}"
  mkdir -p "${output_dir}"
  fastqc --outdir "$output_dir" "$input_directory"/*.fastq.gz -threads 16
  echo "Finished ${sample_id}"
done

echo "DONE!""" "

""" !/bin/bash
SBATCH --job-name=fastQC
SBATCH --nodes=1
SBATCH --account=users
SBATCH --qos=mid_mdbf
SBATCH --partition=mid_mdbf
SBATCH --cpus-per-task=16
SBATCH -o /cta/users/merve.kaftancioglu/alternative_scenario/fastQC_as/out-err/%x-%j-slurm.out
SBATCH -e /cta/users/merve.kaftancioglu/alternative_scenario/fastQC_as/out-err/%x-%j-slurm.err

 Load FastQC module
module load fastqc-0.11.7-gcc-9.2.0-nsjzjof

sample_ids=("h3k27me3_ms_cd4_ENCSR613UFD" "h3k27me3_ms_cd4_ENCSR740SDR" "h3k27me3_ms_cd4_ENCSR779JLY" "h3k27me3_ms_cd4_ENCSR993CTA" "h3k27me3_ms_cd8_ENCSR116FVG")

for sample_id in "${sample_ids[@]}"; do
  echo "Predicting ${sample_id}..."
  input_directory="/cta/users/merve.kaftancioglu/alternative_scenario/trimmomatic_as/output/${sample_id}"
  output_dir="/cta/users/merve.kaftancioglu/alternative_scenario/post-trim_fastQC_as/output/${sample_id}"
  mkdir -p "${output_dir}"
  fastqc --outdir "$output_dir" "$input_directory"/*.fastq.gz -threads 16
  echo "Finished ${sample_id}"
done

echo "DONE!""" "

""" !/bin/bash
SBATCH --job-name=fastQC
SBATCH --nodes=1
SBATCH --account=users
SBATCH --qos=mid_mdbf
SBATCH --partition=mid_mdbf
SBATCH --cpus-per-task=16
SBATCH -o /cta/users/merve.kaftancioglu/alternative_scenario/fastQC_as/out-err/%x-%j-slurm.out
SBATCH -e /cta/users/merve.kaftancioglu/alternative_scenario/fastQC_as/out-err/%x-%j-slurm.err

 Load FastQC module
module load fastqc-0.11.7-gcc-9.2.0-nsjzjof

sample_ids=("h3k27me3_ms_cd8_ENCSR122JCM" "h3k27me3_ms_cd8_ENCSR216ZVA" "h3k27me3_ms_cd8_ENCSR284IKS" "h3k27me3_ms_cd8_ENCSR385BOZ" "h3k27me3_ms_cd8_ENCSR521SFR")

for sample_id in "${sample_ids[@]}"; do
  echo "Predicting ${sample_id}..."
  input_directory="/cta/users/merve.kaftancioglu/alternative_scenario/trimmomatic_as/output/${sample_id}"
  output_dir="/cta/users/merve.kaftancioglu/alternative_scenario/post-trim_fastQC_as/output/${sample_id}"
  mkdir -p "${output_dir}"
  fastqc --outdir "$output_dir" "$input_directory"/*.fastq.gz -threads 16
  echo "Finished ${sample_id}"
done

echo "DONE!""" "

""" !/bin/bash
SBATCH --job-name=fastQC
SBATCH --nodes=1
SBATCH --account=users
SBATCH --qos=mid_mdbf
SBATCH --partition=mid_mdbf
SBATCH --cpus-per-task=16
SBATCH -o /cta/users/merve.kaftancioglu/alternative_scenario/fastQC_as/out-err/%x-%j-slurm.out
SBATCH -e /cta/users/merve.kaftancioglu/alternative_scenario/fastQC_as/out-err/%x-%j-slurm.err

 Load FastQC module
module load fastqc-0.11.7-gcc-9.2.0-nsjzjof

sample_ids=("h3k27me3_ms_nk_ENCSR469QVG" "h3k27me3_ms_nk_ENCSR565WDW" "h3k27me3_normal_b_ENCSR589LHR" "h3k27me3_normal_cd4_ENCSR068YVZ" "h3k27me3_normal_cd8_ENCSR720QRX")

for sample_id in "${sample_ids[@]}"; do
  echo "Predicting ${sample_id}..."
  input_directory="/cta/users/merve.kaftancioglu/alternative_scenario/trimmomatic_as/output/${sample_id}"
  output_dir="/cta/users/merve.kaftancioglu/alternative_scenario/post-trim_fastQC_as/output/${sample_id}"
  mkdir -p "${output_dir}"
  fastqc --outdir "$output_dir" "$input_directory"/*.fastq.gz -threads 16
  echo "Finished ${sample_id}"
done

echo "DONE!""" "

""" !/bin/bash
SBATCH --job-name=fastQC
SBATCH --nodes=1
SBATCH --account=users
SBATCH --qos=mid_mdbf
SBATCH --partition=mid_mdbf
SBATCH --cpus-per-task=16
SBATCH -o /cta/users/merve.kaftancioglu/alternative_scenario/fastQC_as/out-err/%x-%j-slurm.out
SBATCH -e /cta/users/merve.kaftancioglu/alternative_scenario/fastQC_as/out-err/%x-%j-slurm.err

 Load FastQC module
module load fastqc-0.11.7-gcc-9.2.0-nsjzjof

sample_ids=("h3k27me3_normal_nk_ENCSR639NIG" "h3k36me3_ms_b_ENCSR089VQL" "h3k36me3_ms_b_ENCSR119PSR" "h3k36me3_ms_b_ENCSR238WFK" "h3k36me3_ms_b_ENCSR987OPY")

for sample_id in "${sample_ids[@]}"; do
  echo "Predicting ${sample_id}..."
  input_directory="/cta/users/merve.kaftancioglu/alternative_scenario/trimmomatic_as/output/${sample_id}"
  output_dir="/cta/users/merve.kaftancioglu/alternative_scenario/post-trim_fastQC_as/output/${sample_id}"
  mkdir -p "${output_dir}"
  fastqc --outdir "$output_dir" "$input_directory"/*.fastq.gz -threads 16
  echo "Finished ${sample_id}"
done

echo "DONE!""" "

""" !/bin/bash
SBATCH --job-name=fastQC
SBATCH --nodes=1
SBATCH --account=users
SBATCH --qos=mid_mdbf
SBATCH --partition=mid_mdbf
SBATCH --cpus-per-task=16
SBATCH -o /cta/users/merve.kaftancioglu/alternative_scenario/fastQC_as/out-err/%x-%j-slurm.out
SBATCH -e /cta/users/merve.kaftancioglu/alternative_scenario/fastQC_as/out-err/%x-%j-slurm.err

 Load FastQC module
module load fastqc-0.11.7-gcc-9.2.0-nsjzjof

sample_ids=("h3k36me3_ms_cd4_ENCSR276NGH" "h3k36me3_ms_cd4_ENCSR330CQU" "h3k36me3_ms_cd4_ENCSR471WZE" "h3k36me3_ms_cd4_ENCSR473TGK" "h3k36me3_ms_cd4_ENCSR482VIB")

for sample_id in "${sample_ids[@]}"; do
  echo "Predicting ${sample_id}..."
  input_directory="/cta/users/merve.kaftancioglu/alternative_scenario/trimmomatic_as/output/${sample_id}"
  output_dir="/cta/users/merve.kaftancioglu/alternative_scenario/post-trim_fastQC_as/output/${sample_id}"
  mkdir -p "${output_dir}"
  fastqc --outdir "$output_dir" "$input_directory"/*.fastq.gz -threads 16
  echo "Finished ${sample_id}"
done

echo "DONE!""" "

""" !/bin/bash
SBATCH --job-name=fastQC
SBATCH --nodes=1
SBATCH --account=users
SBATCH --qos=mid_mdbf
SBATCH --partition=mid_mdbf
SBATCH --cpus-per-task=16
SBATCH -o /cta/users/merve.kaftancioglu/alternative_scenario/fastQC_as/out-err/%x-%j-slurm.out
SBATCH -e /cta/users/merve.kaftancioglu/alternative_scenario/fastQC_as/out-err/%x-%j-slurm.err

 Load FastQC module
module load fastqc-0.11.7-gcc-9.2.0-nsjzjof

sample_ids=("h3k36me3_ms_cd4_ENCSR532PXR" "h3k36me3_ms_cd4_ENCSR785PKM" "h3k36me3_ms_cd4_ENCSR865FTW" "h3k36me3_ms_cd4_ENCSR898VJE" "h3k36me3_ms_cd8_ENCSR239OWL")

for sample_id in "${sample_ids[@]}"; do
  echo "Predicting ${sample_id}..."
  input_directory="/cta/users/merve.kaftancioglu/alternative_scenario/trimmomatic_as/output/${sample_id}"
  output_dir="/cta/users/merve.kaftancioglu/alternative_scenario/post-trim_fastQC_as/output/${sample_id}"
  mkdir -p "${output_dir}"
  fastqc --outdir "$output_dir" "$input_directory"/*.fastq.gz -threads 16
  echo "Finished ${sample_id}"
done

echo "DONE!""" "

""" !/bin/bash
SBATCH --job-name=fastQC
SBATCH --nodes=1
SBATCH --account=users
SBATCH --qos=mid_mdbf
SBATCH --partition=mid_mdbf
SBATCH --cpus-per-task=16
SBATCH -o /cta/users/merve.kaftancioglu/alternative_scenario/fastQC_as/out-err/%x-%j-slurm.out
SBATCH -e /cta/users/merve.kaftancioglu/alternative_scenario/fastQC_as/out-err/%x-%j-slurm.err

 Load FastQC module
module load fastqc-0.11.7-gcc-9.2.0-nsjzjof

sample_ids=("h3k36me3_ms_cd8_ENCSR435MSB" "h3k36me3_ms_cd8_ENCSR631VIW" "h3k36me3_ms_cd8_ENCSR652WFO" "h3k36me3_ms_cd8_ENCSR757FGN" "h3k36me3_ms_nk_ENCSR158VSE")

for sample_id in "${sample_ids[@]}"; do
  echo "Predicting ${sample_id}..."
  input_directory="/cta/users/merve.kaftancioglu/alternative_scenario/trimmomatic_as/output/${sample_id}"
  output_dir="/cta/users/merve.kaftancioglu/alternative_scenario/post-trim_fastQC_as/output/${sample_id}"
  mkdir -p "${output_dir}"
  fastqc --outdir "$output_dir" "$input_directory"/*.fastq.gz -threads 16
  echo "Finished ${sample_id}"
done

echo "DONE!""" "

""" !/bin/bash
SBATCH --job-name=fastQC
SBATCH --nodes=1
SBATCH --account=users
SBATCH --qos=mid_mdbf
SBATCH --partition=mid_mdbf
SBATCH --cpus-per-task=16
SBATCH -o /cta/users/merve.kaftancioglu/alternative_scenario/fastQC_as/out-err/%x-%j-slurm.out
SBATCH -e /cta/users/merve.kaftancioglu/alternative_scenario/fastQC_as/out-err/%x-%j-slurm.err

 Load FastQC module
module load fastqc-0.11.7-gcc-9.2.0-nsjzjof

sample_ids=("h3k36me3_ms_nk_ENCSR245KON" "h3k36me3_ms_nk_ENCSR530YDY" "h3k36me3_normal_b_ENCSR831AXK" "h3k36me3_normal_cd4_ENCSR640OKB" "h3k36me3_normal_cd8_ENCSR303SQG")

for sample_id in "${sample_ids[@]}"; do
  echo "Predicting ${sample_id}..."
  input_directory="/cta/users/merve.kaftancioglu/alternative_scenario/trimmomatic_as/output/${sample_id}"
  output_dir="/cta/users/merve.kaftancioglu/alternative_scenario/post-trim_fastQC_as/output/${sample_id}"
  mkdir -p "${output_dir}"
  fastqc --outdir "$output_dir" "$input_directory"/*.fastq.gz -threads 16
  echo "Finished ${sample_id}"
done

echo "DONE!""" "

""" !/bin/bash
SBATCH --job-name=fastQC
SBATCH --nodes=1
SBATCH --account=users
SBATCH --qos=mid_mdbf
SBATCH --partition=mid_mdbf
SBATCH --cpus-per-task=16
SBATCH -o /cta/users/merve.kaftancioglu/alternative_scenario/fastQC_as/out-err/%x-%j-slurm.out
SBATCH -e /cta/users/merve.kaftancioglu/alternative_scenario/fastQC_as/out-err/%x-%j-slurm.err

 Load FastQC module
module load fastqc-0.11.7-gcc-9.2.0-nsjzjof

sample_ids=("h3k36me3_normal_nk_ENCSR056GJY" "h3k4me1_ms_b_ENCSR084FXT" "h3k4me1_ms_b_ENCSR480ITK" "h3k4me1_ms_b_ENCSR759DHN" "h3k4me1_ms_b_ENCSR900IAM")

for sample_id in "${sample_ids[@]}"; do
  echo "Predicting ${sample_id}..."
  input_directory="/cta/users/merve.kaftancioglu/alternative_scenario/trimmomatic_as/output/${sample_id}"
  output_dir="/cta/users/merve.kaftancioglu/alternative_scenario/post-trim_fastQC_as/output/${sample_id}"
  mkdir -p "${output_dir}"
  fastqc --outdir "$output_dir" "$input_directory"/*.fastq.gz -threads 16
  echo "Finished ${sample_id}"
done

echo "DONE!""" "

""" !/bin/bash
SBATCH --job-name=fastQC
SBATCH --nodes=1
SBATCH --account=users
SBATCH --qos=mid_mdbf
SBATCH --partition=mid_mdbf
SBATCH --cpus-per-task=16
SBATCH -o /cta/users/merve.kaftancioglu/alternative_scenario/fastQC_as/out-err/%x-%j-slurm.out
SBATCH -e /cta/users/merve.kaftancioglu/alternative_scenario/fastQC_as/out-err/%x-%j-slurm.err

 Load FastQC module
module load fastqc-0.11.7-gcc-9.2.0-nsjzjof

sample_ids=("h3k4me1_ms_b_ENCSR931JYC" "h3k4me1_ms_cd4_ENCSR036YVP" "h3k4me1_ms_cd4_ENCSR043JIA" "h3k4me1_ms_cd4_ENCSR093VUP" "h3k4me1_ms_cd4_ENCSR124IGC")

for sample_id in "${sample_ids[@]}"; do
  echo "Predicting ${sample_id}..."
  input_directory="/cta/users/merve.kaftancioglu/alternative_scenario/trimmomatic_as/output/${sample_id}"
  output_dir="/cta/users/merve.kaftancioglu/alternative_scenario/post-trim_fastQC_as/output/${sample_id}"
  mkdir -p "${output_dir}"
  fastqc --outdir "$output_dir" "$input_directory"/*.fastq.gz -threads 16
  echo "Finished ${sample_id}"
done

echo "DONE!""" "

""" !/bin/bash
SBATCH --job-name=fastQC
SBATCH --nodes=1
SBATCH --account=users
SBATCH --qos=mid_mdbf
SBATCH --partition=mid_mdbf
SBATCH --cpus-per-task=16
SBATCH -o /cta/users/merve.kaftancioglu/alternative_scenario/fastQC_as/out-err/%x-%j-slurm.out
SBATCH -e /cta/users/merve.kaftancioglu/alternative_scenario/fastQC_as/out-err/%x-%j-slurm.err

 Load FastQC module
module load fastqc-0.11.7-gcc-9.2.0-nsjzjof

sample_ids=("h3k4me1_ms_cd4_ENCSR315MYP" "h3k4me1_ms_cd4_ENCSR458QEY" "h3k4me1_ms_cd4_ENCSR485FBT" "h3k4me1_ms_cd4_ENCSR641ZFV" "h3k4me1_ms_cd4_ENCSR687CQX")

for sample_id in "${sample_ids[@]}"; do
  echo "Predicting ${sample_id}..."
  input_directory="/cta/users/merve.kaftancioglu/alternative_scenario/trimmomatic_as/output/${sample_id}"
  output_dir="/cta/users/merve.kaftancioglu/alternative_scenario/post-trim_fastQC_as/output/${sample_id}"
  mkdir -p "${output_dir}"
  fastqc --outdir "$output_dir" "$input_directory"/*.fastq.gz -threads 16
  echo "Finished ${sample_id}"
done

echo "DONE!""" "

""" !/bin/bash
SBATCH --job-name=fastQC
SBATCH --nodes=1
SBATCH --account=users
SBATCH --qos=mid_mdbf
SBATCH --partition=mid_mdbf
SBATCH --cpus-per-task=16
SBATCH -o /cta/users/merve.kaftancioglu/alternative_scenario/fastQC_as/out-err/%x-%j-slurm.out
SBATCH -e /cta/users/merve.kaftancioglu/alternative_scenario/fastQC_as/out-err/%x-%j-slurm.err

 Load FastQC module
module load fastqc-0.11.7-gcc-9.2.0-nsjzjof

sample_ids=("h3k4me1_ms_cd8_ENCSR231XAP" "h3k4me1_ms_cd8_ENCSR572XTB" "h3k4me1_ms_cd8_ENCSR788TEF" "h3k4me1_ms_cd8_ENCSR815YZL" "h3k4me1_ms_cd8_ENCSR861FEC")

for sample_id in "${sample_ids[@]}"; do
  echo "Predicting ${sample_id}..."
  input_directory="/cta/users/merve.kaftancioglu/alternative_scenario/trimmomatic_as/output/${sample_id}"
  output_dir="/cta/users/merve.kaftancioglu/alternative_scenario/post-trim_fastQC_as/output/${sample_id}"
  mkdir -p "${output_dir}"
  fastqc --outdir "$output_dir" "$input_directory"/*.fastq.gz -threads 16
  echo "Finished ${sample_id}"
done

echo "DONE!""" "

""" !/bin/bash
SBATCH --job-name=fastQC
SBATCH --nodes=1
SBATCH --account=users
SBATCH --qos=mid_mdbf
SBATCH --partition=mid_mdbf
SBATCH --cpus-per-task=16
SBATCH -o /cta/users/merve.kaftancioglu/alternative_scenario/fastQC_as/out-err/%x-%j-slurm.out
SBATCH -e /cta/users/merve.kaftancioglu/alternative_scenario/fastQC_as/out-err/%x-%j-slurm.err

 Load FastQC module
module load fastqc-0.11.7-gcc-9.2.0-nsjzjof

sample_ids=("h3k4me1_ms_cd8_ENCSR940PHE" "h3k4me1_ms_nk_ENCSR482UGV" "h3k4me1_ms_nk_ENCSR718LCI" "h3k4me1_ms_nk_ENCSR806HUT" "h3k4me1_normal_b_ENCSR156BXM")

for sample_id in "${sample_ids[@]}"; do
  echo "Predicting ${sample_id}..."
  input_directory="/cta/users/merve.kaftancioglu/alternative_scenario/trimmomatic_as/output/${sample_id}"
  output_dir="/cta/users/merve.kaftancioglu/alternative_scenario/post-trim_fastQC_as/output/${sample_id}"
  mkdir -p "${output_dir}"
  fastqc --outdir "$output_dir" "$input_directory"/*.fastq.gz -threads 16
  echo "Finished ${sample_id}"
done

echo "DONE!""" "

""" !/bin/bash
SBATCH --job-name=fastQC
SBATCH --nodes=1
SBATCH --account=users
SBATCH --qos=mid_mdbf
SBATCH --partition=mid_mdbf
SBATCH --cpus-per-task=16
SBATCH -o /cta/users/merve.kaftancioglu/alternative_scenario/fastQC_as/out-err/%x-%j-slurm.out
SBATCH -e /cta/users/merve.kaftancioglu/alternative_scenario/fastQC_as/out-err/%x-%j-slurm.err

 Load FastQC module
module load fastqc-0.11.7-gcc-9.2.0-nsjzjof

sample_ids=("h3k4me1_normal_cd4_ENCSR102SOR" "h3k4me1_normal_cd8_ENCSR217SHH" "h3k4me1_normal_nk_ENCSR277YKG" "h3k4me3_ms_b_ENCSR260CRI" "h3k4me3_ms_b_ENCSR461QMZ")

for sample_id in "${sample_ids[@]}"; do
  echo "Predicting ${sample_id}..."
  input_directory="/cta/users/merve.kaftancioglu/alternative_scenario/trimmomatic_as/output/${sample_id}"
  output_dir="/cta/users/merve.kaftancioglu/alternative_scenario/post-trim_fastQC_as/output/${sample_id}"
  mkdir -p "${output_dir}"
  fastqc --outdir "$output_dir" "$input_directory"/*.fastq.gz -threads 16
  echo "Finished ${sample_id}"
done

echo "DONE!""" "

""" !/bin/bash
SBATCH --job-name=fastQC
SBATCH --nodes=1
SBATCH --account=users
SBATCH --qos=mid_mdbf
SBATCH --partition=mid_mdbf
SBATCH --cpus-per-task=16
SBATCH -o /cta/users/merve.kaftancioglu/alternative_scenario/fastQC_as/out-err/%x-%j-slurm.out
SBATCH -e /cta/users/merve.kaftancioglu/alternative_scenario/fastQC_as/out-err/%x-%j-slurm.err

 Load FastQC module
module load fastqc-0.11.7-gcc-9.2.0-nsjzjof

sample_ids=("h3k4me3_ms_b_ENCSR530AVK" "h3k4me3_ms_b_ENCSR713ZYF" "h3k4me3_ms_b_ENCSR923EGO" "h3k4me3_ms_cd4_ENCSR180NCM" "h3k4me3_ms_cd4_ENCSR341QLC")

for sample_id in "${sample_ids[@]}"; do
  echo "Predicting ${sample_id}..."
  input_directory="/cta/users/merve.kaftancioglu/alternative_scenario/trimmomatic_as/output/${sample_id}"
  output_dir="/cta/users/merve.kaftancioglu/alternative_scenario/post-trim_fastQC_as/output/${sample_id}"
  mkdir -p "${output_dir}"
  fastqc --outdir "$output_dir" "$input_directory"/*.fastq.gz -threads 16
  echo "Finished ${sample_id}"
done

echo "DONE!""" "

""" !/bin/bash
SBATCH --job-name=fastQC
SBATCH --nodes=1
SBATCH --account=users
SBATCH --qos=mid_mdbf
SBATCH --partition=mid_mdbf
SBATCH --cpus-per-task=16
SBATCH -o /cta/users/merve.kaftancioglu/alternative_scenario/fastQC_as/out-err/%x-%j-slurm.out
SBATCH -e /cta/users/merve.kaftancioglu/alternative_scenario/fastQC_as/out-err/%x-%j-slurm.err

 Load FastQC module
module load fastqc-0.11.7-gcc-9.2.0-nsjzjof

sample_ids=("h3k4me3_ms_cd4_ENCSR482TGI" "h3k4me3_ms_cd4_ENCSR486XJK" "h3k4me3_ms_cd4_ENCSR496LKR" "h3k4me3_ms_cd4_ENCSR603LTN" "h3k4me3_ms_cd4_ENCSR802MXQ")

for sample_id in "${sample_ids[@]}"; do
  echo "Predicting ${sample_id}..."
  input_directory="/cta/users/merve.kaftancioglu/alternative_scenario/trimmomatic_as/output/${sample_id}"
  output_dir="/cta/users/merve.kaftancioglu/alternative_scenario/post-trim_fastQC_as/output/${sample_id}"
  mkdir -p "${output_dir}"
  fastqc --outdir "$output_dir" "$input_directory"/*.fastq.gz -threads 16
  echo "Finished ${sample_id}"
done

echo "DONE!""" "

""" !/bin/bash
SBATCH --job-name=fastQC
SBATCH --nodes=1
SBATCH --account=users
SBATCH --qos=mid_mdbf
SBATCH --partition=mid_mdbf
SBATCH --cpus-per-task=16
SBATCH -o /cta/users/merve.kaftancioglu/alternative_scenario/fastQC_as/out-err/%x-%j-slurm.out
SBATCH -e /cta/users/merve.kaftancioglu/alternative_scenario/fastQC_as/out-err/%x-%j-slurm.err

 Load FastQC module
module load fastqc-0.11.7-gcc-9.2.0-nsjzjof

sample_ids=("h3k4me3_ms_cd4_ENCSR878YHM" "h3k4me3_ms_cd4_ENCSR954ZLD" "h3k4me3_ms_cd8_ENCSR231ZZH" "h3k4me3_ms_cd8_ENCSR278QHR" "h3k4me3_ms_cd8_ENCSR516CKJ")

for sample_id in "${sample_ids[@]}"; do
  echo "Predicting ${sample_id}..."
  input_directory="/cta/users/merve.kaftancioglu/alternative_scenario/trimmomatic_as/output/${sample_id}"
  output_dir="/cta/users/merve.kaftancioglu/alternative_scenario/post-trim_fastQC_as/output/${sample_id}"
  mkdir -p "${output_dir}"
  fastqc --outdir "$output_dir" "$input_directory"/*.fastq.gz -threads 16
  echo "Finished ${sample_id}"
done

echo "DONE!""" "

""" !/bin/bash
SBATCH --job-name=fastQC
SBATCH --nodes=1
SBATCH --account=users
SBATCH --qos=mid_mdbf
SBATCH --partition=mid_mdbf
SBATCH --cpus-per-task=16
SBATCH -o /cta/users/merve.kaftancioglu/alternative_scenario/fastQC_as/out-err/%x-%j-slurm.out
SBATCH -e /cta/users/merve.kaftancioglu/alternative_scenario/fastQC_as/out-err/%x-%j-slurm.err

 Load FastQC module
module load fastqc-0.11.7-gcc-9.2.0-nsjzjof

sample_ids=("h3k4me3_ms_cd8_ENCSR535YYH" "h3k4me3_ms_cd8_ENCSR741XAE" "h3k4me3_ms_cd8_ENCSR848XJL" "h3k4me3_ms_nk_ENCSR234BAD" "h3k4me3_ms_nk_ENCSR703TWY")

for sample_id in "${sample_ids[@]}"; do
  echo "Predicting ${sample_id}..."
  input_directory="/cta/users/merve.kaftancioglu/alternative_scenario/trimmomatic_as/output/${sample_id}"
  output_dir="/cta/users/merve.kaftancioglu/alternative_scenario/post-trim_fastQC_as/output/${sample_id}"
  mkdir -p "${output_dir}"
  fastqc --outdir "$output_dir" "$input_directory"/*.fastq.gz -threads 16
  echo "Finished ${sample_id}"
done

echo "DONE!""" "

""" !/bin/bash
SBATCH --job-name=fastQC
SBATCH --nodes=1
SBATCH --account=users
SBATCH --qos=mid_mdbf
SBATCH --partition=mid_mdbf
SBATCH --cpus-per-task=16
SBATCH -o /cta/users/merve.kaftancioglu/alternative_scenario/fastQC_as/out-err/%x-%j-slurm.out
SBATCH -e /cta/users/merve.kaftancioglu/alternative_scenario/fastQC_as/out-err/%x-%j-slurm.err

 Load FastQC module
module load fastqc-0.11.7-gcc-9.2.0-nsjzjof

sample_ids=("h3k4me3_ms_nk_ENCSR755ZGS" "h3k4me3_normal_b_ENCSR791CAF" "h3k4me3_normal_cd4_ENCSR537KJA" "h3k4me3_normal_cd8_ENCSR123ZAT" "h3k4me3_normal_nk_ENCSR394JFQ")

for sample_id in "${sample_ids[@]}"; do
  echo "Predicting ${sample_id}..."
  input_directory="/cta/users/merve.kaftancioglu/alternative_scenario/trimmomatic_as/output/${sample_id}"
  output_dir="/cta/users/merve.kaftancioglu/alternative_scenario/post-trim_fastQC_as/output/${sample_id}"
  mkdir -p "${output_dir}"
  fastqc --outdir "$output_dir" "$input_directory"/*.fastq.gz -threads 16
  echo "Finished ${sample_id}"
done

echo "DONE!""" "

""" !/bin/bash
SBATCH --job-name=fastQC
SBATCH --nodes=1
SBATCH --account=users
SBATCH --qos=mid_mdbf
SBATCH --partition=mid_mdbf
SBATCH --cpus-per-task=16
SBATCH -o /cta/users/merve.kaftancioglu/alternative_scenario/fastQC_as/out-err/%x-%j-slurm.out
SBATCH -e /cta/users/merve.kaftancioglu/alternative_scenario/fastQC_as/out-err/%x-%j-slurm.err

 Load FastQC module
module load fastqc-0.11.7-gcc-9.2.0-nsjzjof

sample_ids=("h3k9me3_ms_b_ENCSR445DNH" "h3k9me3_ms_b_ENCSR486SZN" "h3k9me3_ms_b_ENCSR983MKX" "h3k9me3_ms_cd4_ENCSR057IZD" "h3k9me3_ms_cd4_ENCSR550DPT")

for sample_id in "${sample_ids[@]}"; do
  echo "Predicting ${sample_id}..."
  input_directory="/cta/users/merve.kaftancioglu/alternative_scenario/trimmomatic_as/output/${sample_id}"
  output_dir="/cta/users/merve.kaftancioglu/alternative_scenario/post-trim_fastQC_as/output/${sample_id}"
  mkdir -p "${output_dir}"
  fastqc --outdir "$output_dir" "$input_directory"/*.fastq.gz -threads 16
  echo "Finished ${sample_id}"
done

echo "DONE!""" "

""" !/bin/bash
SBATCH --job-name=fastQC
SBATCH --nodes=1
SBATCH --account=users
SBATCH --qos=mid_mdbf
SBATCH --partition=mid_mdbf
SBATCH --cpus-per-task=16
SBATCH -o /cta/users/merve.kaftancioglu/alternative_scenario/fastQC_as/out-err/%x-%j-slurm.out
SBATCH -e /cta/users/merve.kaftancioglu/alternative_scenario/fastQC_as/out-err/%x-%j-slurm.err

 Load FastQC module
module load fastqc-0.11.7-gcc-9.2.0-nsjzjof

sample_ids=("h3k9me3_ms_cd4_ENCSR677OEF" "h3k9me3_ms_cd4_ENCSR729ZXH" "h3k9me3_ms_cd4_ENCSR851RJV" "h3k9me3_ms_cd4_ENCSR919DFZ" "h3k9me3_ms_cd4_ENCSR953STZ")

for sample_id in "${sample_ids[@]}"; do
  echo "Predicting ${sample_id}..."
  input_directory="/cta/users/merve.kaftancioglu/alternative_scenario/trimmomatic_as/output/${sample_id}"
  output_dir="/cta/users/merve.kaftancioglu/alternative_scenario/post-trim_fastQC_as/output/${sample_id}"
  mkdir -p "${output_dir}"
  fastqc --outdir "$output_dir" "$input_directory"/*.fastq.gz -threads 16
  echo "Finished ${sample_id}"
done

echo "DONE!""" "

""" !/bin/bash
SBATCH --job-name=fastQC
SBATCH --nodes=1
SBATCH --account=users
SBATCH --qos=mid_mdbf
SBATCH --partition=mid_mdbf
SBATCH --cpus-per-task=16
SBATCH -o /cta/users/merve.kaftancioglu/alternative_scenario/fastQC_as/out-err/%x-%j-slurm.out
SBATCH -e /cta/users/merve.kaftancioglu/alternative_scenario/fastQC_as/out-err/%x-%j-slurm.err

 Load FastQC module
module load fastqc-0.11.7-gcc-9.2.0-nsjzjof

sample_ids=("h3k9me3_ms_cd4_ENCSR959VZU" "h3k9me3_ms_cd8_ENCSR101USF" "h3k9me3_ms_cd8_ENCSR354GNT" "h3k9me3_ms_cd8_ENCSR377QYB" "h3k9me3_ms_cd8_ENCSR733LCG")

for sample_id in "${sample_ids[@]}"; do
  echo "Predicting ${sample_id}..."
  input_directory="/cta/users/merve.kaftancioglu/alternative_scenario/trimmomatic_as/output/${sample_id}"
  output_dir="/cta/users/merve.kaftancioglu/alternative_scenario/post-trim_fastQC_as/output/${sample_id}"
  mkdir -p "${output_dir}"
  fastqc --outdir "$output_dir" "$input_directory"/*.fastq.gz -threads 16
  echo "Finished ${sample_id}"
done

echo "DONE!""" "

""" !/bin/bash
SBATCH --job-name=fastQC
SBATCH --nodes=1
SBATCH --account=users
SBATCH --qos=mid_mdbf
SBATCH --partition=mid_mdbf
SBATCH --cpus-per-task=16
SBATCH -o /cta/users/merve.kaftancioglu/alternative_scenario/fastQC_as/out-err/%x-%j-slurm.out
SBATCH -e /cta/users/merve.kaftancioglu/alternative_scenario/fastQC_as/out-err/%x-%j-slurm.err

 Load FastQC module
module load fastqc-0.11.7-gcc-9.2.0-nsjzjof

sample_ids=("h3k9me3_ms_cd8_ENCSR980KTW" "h3k9me3_ms_nk_ENCSR061ATV" "h3k9me3_ms_nk_ENCSR611YIA" "h3k9me3_ms_nk_ENCSR667HUI" "h3k9me3_normal_b_ENCSR445LTM")

for sample_id in "${sample_ids[@]}"; do
  echo "Predicting ${sample_id}..."
  input_directory="/cta/users/merve.kaftancioglu/alternative_scenario/trimmomatic_as/output/${sample_id}"
  output_dir="/cta/users/merve.kaftancioglu/alternative_scenario/post-trim_fastQC_as/output/${sample_id}"
  mkdir -p "${output_dir}"
  fastqc --outdir "$output_dir" "$input_directory"/*.fastq.gz -threads 16
  echo "Finished ${sample_id}"
done

echo "DONE!""" "

""" !/bin/bash
SBATCH --job-name=fastQC
SBATCH --nodes=1
SBATCH --account=users
SBATCH --qos=mid_mdbf
SBATCH --partition=mid_mdbf
SBATCH --cpus-per-task=16
SBATCH -o /cta/users/merve.kaftancioglu/alternative_scenario/fastQC_as/out-err/%x-%j-slurm.out
SBATCH -e /cta/users/merve.kaftancioglu/alternative_scenario/fastQC_as/out-err/%x-%j-slurm.err

 Load FastQC module
module load fastqc-0.11.7-gcc-9.2.0-nsjzjof

sample_ids=("h3k9me3_normal_cd4_ENCSR433EWI" "h3k9me3_normal_cd8_ENCSR294HTM" "h3k9me3_normal_nk_ENCSR025UNZ")

for sample_id in "${sample_ids[@]}"; do
  echo "Predicting ${sample_id}..."
  input_directory="/cta/users/merve.kaftancioglu/alternative_scenario/trimmomatic_as/output/${sample_id}"
  output_dir="/cta/users/merve.kaftancioglu/alternative_scenario/post-trim_fastQC_as/output/${sample_id}"
  mkdir -p "${output_dir}"
  fastqc --outdir "$output_dir" "$input_directory"/*.fastq.gz -threads 16
  echo "Finished ${sample_id}"
done

echo "DONE!""" "

""" !/bin/bash
SBATCH --job-name=fastQC
SBATCH --nodes=1
SBATCH --account=users
SBATCH --qos=mid_mdbf
SBATCH --partition=mid_mdbf
SBATCH --cpus-per-task=16
SBATCH -o /cta/users/merve.kaftancioglu/alternative_scenario/post-trim_fastQC_as/out-err/%x-%j-slurm.out
SBATCH -e /cta/users/merve.kaftancioglu/alternative_scenario/post-trim_fastQC_as/out-err/%x-%j-slurm.err

 Load FastQC module
module load fastqc-0.11.7-gcc-9.2.0-nsjzjof

sample_ids=("h3k27ac_ms_b_ENCSR295QZX" "h3k27ac_ms_b_ENCSR364OIK" "h3k27ac_ms_b_ENCSR538URI" "h3k27ac_ms_b_ENCSR541PLO" "h3k27ac_ms_b_ENCSR617GMR")

for sample_id in "${sample_ids[@]}"; do
  echo "Predicting ${sample_id}..."
  input_directory="/cta/users/merve.kaftancioglu/alternative_scenario/trimmomatic_as/output/${sample_id}"
  output_dir="/cta/users/merve.kaftancioglu/alternative_scenario/post-trim_fastQC_as/output/${sample_id}"
  mkdir -p "${output_dir}"
   Fix: Expand the variable before using it
  output_dir_expanded="${output_directory}"
  fastqc --outdir "$output_dir_expanded" "$input_directory"/*.fastq.gz -threads 16
  echo "Finished ${sample_id}"
done


echo "DONE!""" "

""" !/bin/bash
SBATCH --job-name=fastQC
SBATCH --nodes=1
SBATCH --account=users
SBATCH --qos=mid_mdbf
SBATCH --partition=mid_mdbf
SBATCH --cpus-per-task=16
SBATCH -o /cta/users/merve.kaftancioglu/alternative_scenario/post-trim_fastQC_as/out-err/%x-%j-slurm.out
SBATCH -e /cta/users/merve.kaftancioglu/alternative_scenario/post-trim_fastQC_as/out-err/%x-%j-slurm.err

 Load FastQC module
module load fastqc-0.11.7-gcc-9.2.0-nsjzjof

sample_ids=("h3k27ac_ms_cd4_ENCSR200SSJ" "h3k27ac_ms_cd4_ENCSR322MTA" "h3k27ac_ms_cd4_ENCSR331WMS" "h3k27ac_ms_cd4_ENCSR350UKV" "h3k27ac_ms_cd4_ENCSR474PYR")

for sample_id in "${sample_ids[@]}"; do
  echo "Predicting ${sample_id}..."
  input_directory="/cta/users/merve.kaftancioglu/alternative_scenario/trimmomatic_as/output/${sample_id}"
  output_dir="/cta/users/merve.kaftancioglu/alternative_scenario/post-trim_fastQC_as/output/${sample_id}"
  mkdir -p "${output_dir}"
  fastqc --outdir "$output_dir" "$input_directory"/*.fastq.gz -threads 16
  echo "Finished ${sample_id}"
done

echo "DONE!""" "


""" !/bin/bash
SBATCH --job-name=fastQC
SBATCH --nodes=1
SBATCH --account=users
SBATCH --qos=mid_mdbf
SBATCH --partition=mid_mdbf
SBATCH --cpus-per-task=16
SBATCH -o /cta/users/merve.kaftancioglu/alternative_scenario/fastQC_as/out-err/%x-%j-slurm.out
SBATCH -e /cta/users/merve.kaftancioglu/alternative_scenario/fastQC_as/out-err/%x-%j-slurm.err


 Load FastQC module
module load fastqc-0.11.7-gcc-9.2.0-nsjzjof

substrate_ids=("h3k36me3_ms_cd8_ENCSR652WFO" "h3k36me3_ms_cd8_ENCSR757FGN" "h3k36me3_ms_nk_ENCSR158VSE" "h3k36me3_normal_cd8_ENCSR303SQG" "h3k4me1_ms_b_ENCSR084FXT" "h3k4me1_ms_cd4_ENCSR315MYP" "h3k4me1_ms_cd8_ENCSR815YZL" "h3k4me1_ms_cd8_ENCSR861FEC" "h3k4me1_ms_nk_ENCSR718LCI" "h3k4me1_ms_nk_ENCSR806HUT" "h3k4me1_normal_b_ENCSR156BXM" "h3k4me3_ms_b_ENCSR923EGO" "h3k4me3_ms_cd4_ENCSR180NCM" "h3k4me3_ms_cd4_ENCSR341QLC" "h3k4me3_ms_cd4_ENCSR603LTN" "h3k4me3_ms_cd4_ENCSR802MXQ" "h3k4me3_ms_cd8_ENCSR278QHR" "h3k4me3_ms_cd8_ENCSR516CKJ" "h3k4me3_ms_nk_ENCSR234BAD" "h3k4me3_ms_nk_ENCSR703TWY" "h3k4me3_normal_cd8_ENCSR123ZAT" "h3k4me3_normal_nk_ENCSR394JFQ" "h3k9me3_ms_cd4_ENCSR057IZD" "h3k9me3_ms_cd4_ENCSR550DPT" "h3k9me3_ms_nk_ENCSR611YIA" "h3k9me3_ms_nk_ENCSR667HUI" "h3k9me3_normal_b_ENCSR445LTM" "h3k9me3_normal_nk_ENCSR025UNZ")

for substrate_id in "${substrate_ids[@]}"; do
    echo "Predicting ${substrate_id}..."
    input_directory="/cta/users/merve.kaftancioglu/alternative_scenario/seqs/${substrate_id}"
    output_dir="/cta/users/merve.kaftancioglu/alternative_scenario/seqs/${substrate_id}"    
    mkdir -p "${output_dir}"
    fastqc --outdir "$output_directory" "$input_directory"/*.fastq.gz -threads 16
    echo "Finished ${substrate_id}"
done


echo "DONE!""" "



""" !/bin/bash
SBATCH --job-name=bwa_single
SBATCH --nodes=1
SBATCH --account=users
SBATCH --qos=cuda
SBATCH --partition=cuda
SBATCH --cpus-per-task=16
SBATCH -o /cta/users/merve.kaftancioglu/alternative_scenario/bwa_as/out-err/%x-%j-slurm.out
SBATCH -e /cta/users/merve.kaftancioglu/alternative_scenario/bwa_as/out-err/%x-%j-slurm.err

 Load BWA module
module load bwa-0.7.17-gcc-9.2.0-xwkclqy

 Define directories (consider adjusting paths)
input_dir="/cta/users/merve.kaftancioglu/alternative_scenario/trimmomatic_as/output"
output_dir="/cta/users/merve.kaftancioglu/alternative_scenario/bwa_as/output"
reference_genome="/cta/users/merve.kaftancioglu/alternative_scenario/bwa_as/reference/human_g1k_v37.fasta"

 Check if BWA index exists
if [ ! -f "$reference_genome.bwt" ]; then
  echo "Error: BWA index not found. Please create the index using: bwa mem index $reference_genome"
  exit 1
fi

 Define sample IDs (modify as needed)
sample_ids=("h3k36me_ms_cd4_ENCSR532PXR")

 Loop through all sample subdirectories
for sample_dir in "$input_dir"/*; do
   Check if directory name matches any sample ID
  if [[ " ${sample_ids[@]} " =~ " $(basename "$sample_dir") " ]]; then
    echo "Processing sample: $(basename "$sample_dir")"

     Create output directory for the sample (if it doesn't exist)
    mkdir -p "$output_dir/$(basename "$sample_dir")"

     Find FASTQ files within the sample directory
     Adjust the pattern below if your file names differ
    fastq_files=( "$sample_dir"/*.fastq.gz )

     Check if any FASTQ files found
    if [[ ${fastq_files[@]} -eq 0 ]]; then
      echo "Warning: No FASTQ files found in $sample_dir"
      continue  Skip to next iteration if no FASTQ files
    fi

     BWA mem command with paired-end or single-end handling
    if [[ ${fastq_files[@]} -eq 2 ]]; then
       Paired-end reads (assuming *_1.fastq.gz and *_2.fastq.gz naming)
      bwa mem $reference_genome ${fastq_files[0]} ${fastq_files[1]} > "$output_dir/$(basename "$sample_dir").sam"
    else
       Single-end read (assuming only one FASTQ file)
       Extract only the sample name from the directory path
      sample_name=$(basename "$sample_dir" /)
      bwa mem $reference_genome ${fastq_files[0]} > "$output_dir/$sample_name.sam"
    fi
  fi
done


""" !/bin/bash
SBATCH --job-name=sam_to_bam
SBATCH --nodes=1
SBATCH --account=users
SBATCH --qos=cuda
SBATCH --partition=cuda
SBATCH --cpus-per-task=16
SBATCH -o /cta/users/merve.kaftancioglu/alternative_scenario/samtools_as/out-err/%x-%j-slurm.out
SBATCH -e /cta/users/merve.kaftancioglu/alternative_scenario/samtools_as/out-err/%x-%j-slurm.err

 Load samtools module
module load samtools-1.9-gcc-9.2.0-w7pulwi

 Define directories (consider adjusting paths)
input_dir="/cta/users/merve.kaftancioglu/alternative_scenario/bwa_as/output"
output_dir="/cta/users/merve.kaftancioglu/alternative_scenario/samtools_as/output"

 Loop through all subdirectories in the input directory
for subdir in $input_dir/*; do
   Check if it's a directory (avoid hidden folders)
  if [ -d "$subdir" ]; then
     Get the sample name from the subdirectory name
    sample_name=$(basename "$subdir")

     Look for a .sam file within the subdirectory
    sam_file="$subdir/$sample_name.sam"

     Check if the SAM file exists
    if [ -f "$sam_file" ]; then
       Convert SAM to BAM using samtools view
      samtools view -bS "$sam_file" > "$output_dir/$sample_name.bam"

       Check for conversion errors (optional)
      if [ $? -ne 0 ]; then
        echo "Error converting $sam_file to BAM. See slurm error logs for details."
      fi
    fi
  fi
done

echo "DONE!""" "

""" !/bin/bash
SBATCH --job-name=star_alignment
SBATCH --nodes=1
SBATCH --account=users
SBATCH --qos=mid_mdbf
SBATCH --partition=mid_mdbf
SBATCH --cpus-per-task=16
SBATCH -o /cta/users/merve.kaftancioglu/alternative_scenario/bowtie2_as/out-err/%x-%j-slurm.out
SBATCH -e /cta/users/merve.kaftancioglu/alternative_scenario/bowtie2_as/out-err/%x-%j-slurm.err

 Load STAR module
module load star-2.7.0e-gcc-9.2.0-vynasg3

 Define directories
genome_dir="/cta/users/merve.kaftancioglu/alternative_scenario/star_as/index"   Update with your genome directory
input_dir="/cta/users/merve.kaftancioglu/alternative_scenario/trimmomatic_as/output"   Update with your input directory
output_dir="/cta/users/merve.kaftancioglu/alternative_scenario/star_as/output"   Update with your output directory

 Run STAR for each folder within the input directory
for sample_dir in "$input_dir"/*; do
   Check if it's a directory (avoid hidden files or non-directories)
  if [[ -d "$sample_dir" ]]; then
     Extract sample ID from directory name (modify if needed)
    sample_id=$(basename "$sample_dir")

     Define input and output paths with sample ID
    input_fastq="$sample_dir"/*.fastq.gz
    output_bam="$output_dir/$sample_id/$sample_id.bam"

     Run STAR for the current sample
    STAR --genomeDir "$genome_dir" \
         --runThreadN 16 \   Adjust threads based on your resource allocation
         --readFilesIn "$input_fastq" \   No colon after flag name
         --outFileNamePrefix "$output_dir/$sample_id/" \   Include trailing slash for proper path
         --outSAMtype BAM SortedByCoordinate \   No colon after flag name
         --outSAMunmapped Within \
         --outSAMattributes Standard \
         --alignEndsType EndToEnd
  fi
done

""" !/bin/bash
SBATCH --job-name=move
SBATCH --nodes=1
SBATCH --account=users
SBATCH --qos=cuda
SBATCH --partition=cuda
SBATCH --cpus-per-task=16
SBATCH -o /cta/users/merve.kaftancioglu/alternative_scenario/fastQC_as/out-err/%x-%j-slurm.out
SBATCH -e /cta/users/merve.kaftancioglu/alternative_scenario/fastQC_as/out-err/%x-%j-slurm.err

 Dry run (set to 'cp' for actual move)
DRY_RUN="cp"  

 Replace with the actual source directory path
 Use curly brackets only if $sample_id contains special characters or spaces
SOURCE_DIR="/cta/users/merve.kaftancioglu/alternative_scenario/trimmomatic_as/output/{$sample_id}"

 Destination directory 
DEST_DIR="/cta/users/merve.kaftancioglu/alternative_scenario/post-trim_fastQC_as/output/{$sample_id}"

 Create destination directory if it doesn't exist
mkdir -p "$DEST_DIR"

 Find zip files (ensure there are some in the source directory)
zip_files=$(find "$SOURCE_DIR" -type f -name "*.zip")

 Check if any zip files were found
if [ -z "$zip_files" ]; then
  echo "No zip files found in source directory: $SOURCE_DIR"
  exit 0
fi

 Loop through zip files and copy them
for file in $zip_files
do
  $DRY_RUN cp "$file" "$DEST_DIR"
done

 ~/.bashrc: executed by bash(1) for non-login shells.
 see /usr/share/doc/bash/examples/startup-files (in the package bash-doc)
 for examples

 If not running interactively, don't do anything
case $- in
    *i*) ;;
      *) return;;
esac

 don't put duplicate lines or lines starting with space in the history.
 See bash(1) for more options
HISTCONTROL=ignoreboth

 append to the history file, don't overwrite it
shopt -s histappend

 for setting history length see HISTSIZE and HISTFILESIZE in bash(1)
HISTSIZE=1000
HISTFILESIZE=2000

 check the window size after each command and, if necessary,
 update the values of LINES and COLUMNS.
shopt -s checkwinsize

 If set, the pattern "**" used in a pathname expansion context will
 match all files and zero or more directories and subdirectories.
shopt -s globstar

 make less more friendly for non-text input files, see lesspipe(1)
[ -x /usr/bin/lesspipe ] && eval "$(SHELL=/bin/sh lesspipe)"

 set variable identifying the chroot you work in (used in the prompt below)
if [ -z "${debian_chroot:-}" ] && [ -r /etc/debian_chroot ]; then
    debian_chroot=$(cat /etc/debian_chroot)
fi

 set a fancy prompt (non-color, unless we know we "want" color)
case "$TERM" in
    xterm-color|*-256color) color_prompt=yes;;
esac

 uncomment for a colored prompt, if the terminal has the capability; turned
 off by default to not distract the user: the focus in a terminal window
 should be on the output of commands, not on the prompt
force_color_prompt=yes

if [ -n "$force_color_prompt" ]; then
    if [ -x /usr/bin/tput ] && tput setaf 1 >&/dev/null; then
	 We have color support; assume it's compliant with Ecma-48
	 (ISO/IEC-6429). (Lack of such support is extremely rare, and such
	 a case would tend to support setf rather than setaf.)
	color_prompt=yes
    else
	color_prompt=
    fi
fi

if [ "$color_prompt" = yes ]; then
    PS1='${debian_chroot:+($debian_chroot)}\[\033[01;32m\]\u@\h\[\033[00m\]:\[\033[01;34m\]\w\[\033[00m\]\$ '
else
    PS1='${debian_chroot:+($debian_chroot)}\u@\h:\w\$ '
fi
unset color_prompt force_color_prompt

 If this is an xterm set the title to user@host:dir
case "$TERM" in
xterm*|rxvt*)
    PS1="\[\e]0;${debian_chroot:+($debian_chroot)}\u@\h: \w\a\]$PS1"
    ;;
*)
    ;;
esac

 enable color support of ls and also add handy aliases
if [ -x /usr/bin/dircolors ]; then
    test -r ~/.dircolors && eval "$(dircolors -b ~/.dircolors)" || eval "$(dircolors -b)"
    alias ls='ls --color=auto'
    alias dir='dir --color=auto'
    alias vdir='vdir --color=auto'

    alias grep='grep --color=auto'
    alias fgrep='fgrep --color=auto'
    alias egrep='egrep --color=auto'
fi

 colored GCC warnings and errors
export GCC_COLORS='error=01;31:warning=01;35:note=01;36:caret=01;32:locus=01:quote=01'

 some more ls aliases
alias ll='ls -alF'
alias la='ls -A'
alias l='ls -CF'

 Add an "alert" alias for long running commands.  Use like so:
   sleep 10; alert
alias alert='notify-send --urgency=low -i "$([ $? = 0 ] && echo terminal || echo error)" "$(history|tail -n1|sed -e '\''s/^\s*[0-9]\+\s*//;s/[;&|]\s*alert$//'\'')"'

 Alias definitions.
 You may want to put all your additions into a separate file like
 ~/.bash_aliases, instead of adding them here directly.
 See /usr/share/doc/bash-doc/examples in the bash-doc package.

if [ -f ~/.bash_aliases ]; then
    . ~/.bash_aliases
fi

 enable programmable completion features (you don't need to enable
 this, if it's already enabled in /etc/bash.bashrc and /etc/profile
 sources /etc/bash.bashrc).
if ! shopt -oq posix; then
  if [ -f /usr/share/bash-completion/bash_completion ]; then
    . /usr/share/bash-completion/bash_completion
  elif [ -f /etc/bash_completion ]; then
    . /etc/bash_completion
  fi
fi

 >>> conda initialize >>>
 !! Contents within this block are managed by 'conda init' !!
__conda_setup="$('/home/merve/anaconda3/bin/conda' 'shell.bash' 'hook' 2> /dev/null)"
if [ $? -eq 0 ]; then
    eval "$__conda_setup"
else
    if [ -f "/home/merve/anaconda3/etc/profile.d/conda.sh" ]; then
        . "/home/merve/anaconda3/etc/profile.d/conda.sh"
    else
        export PATH="/home/merve/anaconda3/bin:$PATH"
    fi
fi
unset __conda_setup
 <<< conda initialize <<<


""" !/bin/bash
SBATCH --job-name=star_single
SBATCH --nodes=1
SBATCH --account=users
SBATCH --qos=mid_mdbf
SBATCH --partition=mid_mdbf
SBATCH --cpus-per-task=16
SBATCH -o /cta/users/merve.kaftancioglu/alternative_scenario/star_as/out-err/%x-%j-slurm.out
SBATCH -e /cta/users/merve.kaftancioglu/alternative_scenario/star_as/out-err/%x-%j-slurm.err

 Define paths 
genome_fasta="/cta/users/merve.kaftancioglu/alternative_scenario/star_as/reference/human_g1k_v37.fasta"
star_index_dir="/cta/users/merve.kaftancioglu/alternative_scenario/star_as/index"

 Check if genome fasta file exists
if [ ! -f "$genome_fasta" ]; then
  echo "Error: Genome fasta file not found at: $genome_fasta"
  exit 1
fi

 Check if star_index_dir exists (and create if not)
if [ ! -d "$star_index_dir" ]; then
  echo "Creating star_genome_index directory: $star_index_dir"
  mkdir -p "$star_index_dir"
fi

 Define number of threads (adjust based on your system)
threads=8

 STAR command
star_cmd="STAR \
  --runThreadN $threads \
  --runMode genomeGenerate \
  --genomeDir $star_index_dir \
  --genomeFastaFiles $genome_fasta"

 Run STAR and handle potential errors
echo "Building STAR genome index..."
$star_cmd

if [ $? -ne 0 ]; then
  echo "Error: STAR failed to generate genome index."
  exit 1
fi

echo "STAR genome index created successfully in: $star_index_dir"

""" !/bin/bash
SBATCH --job-name=bowtie2_alignment
SBATCH --nodes=1
SBATCH --account=users
SBATCH --qos=mid_mdbf
SBATCH --partition=mid_mdbf
SBATCH --cpus-per-task=16
SBATCH -o /cta/users/merve.kaftancioglu/alternative_scenario/bowtie2_as/out-err/%x-%j-slurm.out
SBATCH -e /cta/users/merve.kaftancioglu/alternative_scenario/bowtie2_as/out-err/%x-%j-slurm.err

 Load Bowtie2 module
module load bowtie2-2.3.5.1-gcc-9.2.0-ae4dozk

 Sample IDs 
sample_ids=("h3k27ac_ms_b_ENCSR295QZX" "h3k27ac_ms_b_ENCSR364OIK" "h3k27ac_ms_b_ENCSR538URI")

 Directories
input_dir="/cta/users/merve.kaftancioglu/alternative_scenario/trimmomatic_as/output"
output_dir="/cta/users/merve.kaftancioglu/alternative_scenario/bowtie2_as/output"
index_dir="/cta/users/merve.kaftancioglu/alternative_scenario/bowtie2_as/index"  

 Create Bowtie2 index 
bowtie2-build "$reference_genome" "$index_dir"

 Loop through samples
for sample_id in "${sample_ids[@]}"; do
  echo "Processing ${sample_id}..."

   Create output directory 
  mkdir -p "${output_dir}/${sample_id}"

   Combine input and output paths with sample ID
  input_fastq="${input_dir}/${sample_id}"/*.fastq.gz
  output_bam="${output_dir}/${sample_id}/${sample_id}.bam"

   Run Bowtie2 with trimming
  bowtie2 --threads 16 -p 16 --trim3 30 -x "$index_dir" -U "$input_fastq" | samtools view -bS - > "$output_bam"
  
   Convert BAM to SAM
  input_bam="$output_bam"
  output_sam="${output_dir}/${sample_id}/${sample_id}.sam"
  samtools view -h -o "$output_sam" "$input_b

  echo "Finished processing ${sample_id}"
done

echo "All samples processed!"

""" !/bin/bash
SBATCH --job-name=star_single
SBATCH --nodes=1
SBATCH --account=users
SBATCH --qos=mid_mdbf
SBATCH --partition=mid_mdbf
SBATCH --cpus-per-task=16
SBATCH -o /cta/users/merve.kaftancioglu/alternative_scenario/bwa_as/out-err/%x-%j-slurm.out
SBATCH -e /cta/users/merve.kaftancioglu/alternative_scenario/bwa_as/out-err/%x-%j-slurm.err

 Load BWA module
module load bwa-0.7.17-gcc-9.2.0-xwkclqy

 Define directories
input_dir="/cta/users/merve.kaftancioglu/alternative_scenario/trimmomatic_as/output"   Folder containing sample subdirectories
output_dir="/cta/users/merve.kaftancioglu/alternative_scenario/bwa_as/output"   Output directory

 Path to reference genome FASTA file
reference_genome="/cta/users/merve.kaftancioglu/alternative_scenario/bwa_as/reference/human_g1k_v37.fasta"

 Create BWA index 
bwa index $reference_genome

 Loop through all sample subdirectories
for sample_dir in "$input_dir"/*; do
   Extract sample ID from directory name (assuming it matches)
  sample_id=$(basename "$sample_dir")

   Check if it's a directory (skip non-directories)
  if [[ -d "$sample_dir" ]]; then
    echo "Processing sample: $sample_id"

     Find FASTQ files within the sample directory
     Adjust the pattern below if your file names differ
    fastq_files=( "$sample_dir"/*.fastq.gz )

     Check if any FASTQ files found
    if [[ ${fastq_files[@]} -eq 0 ]]; then
      echo "Warning: No FASTQ files found in $sample_dir"
      continue   Skip to next iteration if no FASTQ files
    fi

     BWA command with paired-end or single-end handling
    if [[ ${fastq_files[@]} -eq 2 ]]; then
       Paired-end reads (assuming *_1.fastq.gz and *_2.fastq.gz naming)
      bwa mem $reference_genome ${fastq_files[0]} ${fastq_files[1]} > "$output_dir/$sample_id.sam"
    else
       Single-end read (assuming only one FASTQ file)
      bwa mem $reference_genome ${fastq_files[0]} > "$output_dir/$sample_id.sam"
    fi
  fi
done

 BWA Bwasw 
bwa bwasw index_prefix input_reads.fasta -t $SLURM_NTASKS_PER_NODE > bwa_bwasw_alignments.sam

""" !/bin/bash
SBATCH --job-name=bwa_single
SBATCH --nodes=1
SBATCH --account=users
SBATCH --qos=cuda
SBATCH --partition=cuda
SBATCH --cpus-per-task=16
SBATCH -o /cta/users/merve.kaftancioglu/alternative_scenario/bwa_as/out-err/%x-%j-slurm.out
SBATCH -e /cta/users/merve.kaftancioglu/alternative_scenario/bwa_as/out-err/%x-%j-slurm.err

 Load BWA module
module load bwa-0.7.17-gcc-9.2.0-xwkclqy

 Define directories
input_dir="/cta/users/merve.kaftancioglu/alternative_scenario/trimmomatic_as/output"   Folder containing sample subdirectories
output_dir="/cta/users/merve.kaftancioglu/alternative_scenario/bwa_as/output"   Output directory

 Path to reference genome FASTA file
reference_genome="/cta/users/merve.kaftancioglu/alternative_scenario/bwa_as/reference/human_g1k_v37.fasta"

 Assuming you already have the BWA index files created (remove the commented line below if you need to create them)
 bwa mem index $reference_genome   Only uncomment this line if you need to recreate the index

 Loop through all sample subdirectories
for sample_dir in "$input_dir"/*; do
   Extract sample ID from directory name (assuming it matches)
  sample_id=$(basename "$sample_dir")

   Check if it's a directory (skip non-directories)
  if [[ -d "$sample_dir" ]]; then
    echo "Processing sample: $sample_id"

     Create output directory for the sample (if it doesn't exist)
    mkdir -p "$output_dir/$sample_id"   Create directory with sample ID in the output directory

     Find FASTQ files within the sample directory
     Adjust the pattern below if your file names differ
    fastq_files=( "$sample_dir"/*.fastq.gz )

     Check if any FASTQ files found
    if [[ ${fastq_files[@]} -eq 0 ]]; then
      echo "Warning: No FASTQ files found in $sample_dir"
      continue   Skip to next iteration if no FASTQ files
    fi

     BWA mem command with paired-end or single-end handling
    if [[ ${fastq_files[@]} -eq 2 ]]; then
       Paired-end reads (assuming *_1.fastq.gz and *_2.fastq.gz naming)
      bwa mem $reference_genome ${fastq_files[0]} ${fastq_files[1]} > "$output_dir/$sample_id.sam"
    else
       Single-end read (assuming only one FASTQ file)
      bwa mem $reference_genome ${fastq_files[0]} > "$output_dir/$sample_dir.sam"   Use sample directory name for single-end output
    fi
  fi
done


""" !/bin/bash
SBATCH --job-name=bwa_single
SBATCH --nodes=1
SBATCH --account=users
SBATCH --qos=cuda
SBATCH --partition=cuda
SBATCH --cpus-per-task=16
SBATCH -o /cta/users/merve.kaftancioglu/alternative_scenario/bwa_as/out-err/%x-%j-slurm.out
SBATCH -e /cta/users/merve.kaftancioglu/alternative_scenario/bwa_as/out-err/%x-%j-slurm.err

 Load BWA module
module load bwa-0.7.17-gcc-9.2.0-xwkclqy

 Define directories (consider adjusting paths)
input_dir="/cta/users/merve.kaftancioglu/alternative_scenario/trimmomatic_as/output"
output_dir="/cta/users/merve.kaftancioglu/alternative_scenario/bwa_as/output"
reference_genome="/cta/users/merve.kaftancioglu/alternative_scenario/bwa_as/reference/human_g1k_v37.fasta"

 Check if BWA index exists
if [ ! -f "$reference_genome.bwt" ]; then
  echo "Error: BWA index not found. Please create the index using: bwa mem index $reference_genome"
  exit 1
fi

 Define sample IDs (modify as needed)
sample_ids=("h3k27ac_ms_b_ENCSR295QZX" "h3k27ac_ms_b_ENCSR538URI" "h3k27ac_ms_b_ENCSR617GMR" "h3k4me1_ms_b_ENCSR084FXT" "h3k4me1_ms_b_ENCSR480ITK" "h3k4me1_ms_b_ENCSR759DHN" "h3k4me1_ms_b_ENCSR900IAM" "h3k4me1_ms_b_ENCSR931JYC" "h3k4me3_ms_b_ENCSR260CRI" "h3k4me3_ms_b_ENCSR461QMZ" "h3k4me3_ms_b_ENCSR530AVK" "h3k4me3_ms_b_ENCSR713ZYF" "h3k4me3_ms_b_ENCSR923EGO" "h3k9me3_ms_b_ENCSR445DNH" "h3k9me3_ms_b_ENCSR486SZN" "h3k9me3_ms_b_ENCSR983MKX" "h3k27ac_ms_b_ENCSR364OIK" "h3k27ac_ms_b_ENCSR541PLO" "h3k36me3_ms_cd4_ENCSR471WZE" "h3k36me3_ms_cd4_ENCSR473TGK" "h3k4me1_ms_nk_ENCSR482UGV" "h3k4me1_ms_nk_ENCSR718LCI" "h3k4me1_ms_nk_ENCSR806HUT" "h3k4me3_ms_nk_ENCSR234BAD" "h3k4me3_ms_nk_ENCSR703TWY" "h3k4me3_ms_nk_ENCSR75")

 Loop through all sample subdirectories
for sample_dir in "$input_dir"/*; do
   Check if directory name matches any sample ID
  if [[ " ${sample_ids[@]} " =~ " $(basename "$sample_dir") " ]]; then
    echo "Processing sample: $(basename "$sample_dir")"

     Create output directory for the sample (if it doesn't exist)
    mkdir -p "<span class="math-inline">output\_dir/</span>(basename "$sample_dir")"

     Find FASTQ files within the sample directory
     Adjust the pattern below if your file names differ
    fastq_files=( "$sample_dir"/*.fastq.gz )

     Check if any FASTQ files found
    if [[ ${fastq_files[@]} -eq 0 ]]; then
      echo "Warning: No FASTQ files found in $sample_dir"
      continue  Skip to next iteration if no FASTQ files
    fi

     BWA mem command with paired-end or single-end handling
    if [[ ${fastq_files[@]} -eq 2 ]]; then
       Paired-end reads (assuming *_1.fastq.gz and *_2.fastq.gz naming)
      bwa mem $reference_genome ${fastq_files[0]} ${fastq_files[1]} > "$output_dir/$sample_id.sam"
    else
       Single-end read (assuming only one FASTQ file)
       Extract only the sample name from the directory path
      sample_name=$(basename "$sample_dir" /)
      bwa mem $reference_genome ${fastq_files[0]} > "$output_dir/$sample_name.sam"
    fi

     Convert SAM to BAM format
    samtools view -S -b "$output_dir/$sample_id.sam" > "$output_dir/$sample_id.bam"
    
     Remove the intermediate SAM file (optional, saves space)
    rm "$output_dir/$sample_id.sam"   Uncomment this line if you want to remove SAM files
  fi
done

echo "BWA alignment completed! Both SAM and BAM files generated."

""" !/bin/bash
SBATCH --job-name=bwa_single
SBATCH --nodes=1
SBATCH --account=users
SBATCH --qos=cuda
SBATCH --partition=cuda
SBATCH --cpus-per-task=16
SBATCH -o /cta/users/merve.kaftancioglu/alternative_scenario/bwa_as/out-err/%x-%j-slurm.out
SBATCH -e /cta/users/merve.kaftancioglu/alternative_scenario/bwa_as/out-err/%x-%j-slurm.err

 Load BWA module
module load bwa-0.7.17-gcc-9.2.0-xwkclqy

 Define directories (consider adjusting paths)
input_dir="/cta/users/merve.kaftancioglu/alternative_scenario/trimmomatic_as/output"
output_dir="/cta/users/merve.kaftancioglu/alternative_scenario/bwa_as/output"
reference_genome="/cta/users/merve.kaftancioglu/alternative_scenario/bwa_as/reference/human_g1k_v37.fasta"

 Check if BWA index exists
if [ ! -f "$reference_genome.bwt" ]; then
  echo "Error: BWA index not found. Please create the index using: bwa mem index $reference_genome"
  exit 1
fi

 Define sample IDs (modify as needed)
sample_ids=("h3k27ac_ms_b_ENCSR295QZX" "h3k27ac_ms_b_ENCSR538URI" "h3k27ac_ms_b_ENCSR617GMR" "h3k4me1_ms_b_ENCSR084FXT" "h3k4me1_ms_b_ENCSR480ITK" "h3k4me1_ms_b_ENCSR759DHN" "h3k4me1_ms_b_ENCSR900IAM" "h3k4me1_ms_b_ENCSR931JYC" "h3k4me3_ms_b_ENCSR260CRI" "h3k4me3_ms_b_ENCSR461QMZ" "h3k4me3_ms_b_ENCSR530AVK" "h3k4me3_ms_b_ENCSR713ZYF" "h3k4me3_ms_b_ENCSR923EGO" "h3k9me3_ms_b_ENCSR445DNH" "h3k9me3_ms_b_ENCSR486SZN" "h3k9me3_ms_b_ENCSR983MKX" "h3k27ac_ms_b_ENCSR364OIK" "h3k27ac_ms_b_ENCSR541PLO" "h3k36me3_ms_cd4_ENCSR471WZE" "h3k36me3_ms_cd4_ENCSR473TGK" "h3k4me1_ms_nk_ENCSR482UGV" "h3k4me1_ms_nk_ENCSR718LCI" "h3k4me1_ms_nk_ENCSR806HUT" "h3k4me3_ms_nk_ENCSR234BAD" "h3k4me3_ms_nk_ENCSR703TWY" "h3k4me3_ms_nk_ENCSR75")

 Loop through all sample subdirectories
for sample_dir in "$input_dir"/*; do
   Check if directory name matches any sample ID
  if [[ " ${sample_ids[@]} " =~ " $(basename "$sample_dir") " ]]; then
    echo "Processing sample: $(basename "$sample_dir")"

     Create output directory for the sample (if it doesn't exist)
    mkdir -p "$output_dir/$(basename "$sample_dir")"

     Find FASTQ files within the sample directory
     Adjust the pattern below if your file names differ
    fastq_files=( "$sample_dir"/*.fastq.gz )

     Check if any FASTQ files found
    if [[ ${fastq_files[@]} -eq 0 ]]; then
      echo "Warning: No FASTQ files found in $sample_dir"
      continue  Skip to next iteration if no FASTQ files
    fi

     BWA mem command with paired-end or single-end handling
    if [[ ${fastq_files[@]} -eq 2 ]]; then
       Paired-end reads (assuming *_1.fastq.gz and *_2.fastq.gz naming)
      sample_id=$(basename "$sample_dir")
      bwa mem "$reference_genome" "${fastq_files[0]}" "${fastq_files[1]}" > "$output_dir/$sample_id.sam"
    else
       Single-end read (assuming only one FASTQ file)
       Extract only the sample name from the directory path
      sample_name=$(basename "$sample_dir")
      bwa mem "$reference_genome" "${fastq_files[0]}" > "$output_dir/$sample_name.sam"
    fi

     Convert SAM to BAM format
    samtools view -S -b "$output_dir/$sample_id.sam" > "$output_dir/$sample_id.bam"
    
     Remove the intermediate SAM file (optional, saves space)
    rm "$output_dir/$sample_id.sam"   Uncomment this line if you want to remove SAM files
  fi
done

echo "BWA alignment completed! Both SAM and BAM files generated."

""" !/bin/bash
SBATCH --job-name=fastQC
SBATCH --nodes=1
SBATCH --account=users
SBATCH --qos=short_mdbf
SBATCH --partition=short_mdbf
SBATCH --cpus-per-task=8
SBATCH -o /cta/users/merve.kaftancioglu/fastQC_M2/job_outputs/%x-%j-slurm.out
SBATCH -e /cta/users/merve.kaftancioglu/fastQC_M2/job_outputs/%x-%j-slurm.err


input_directory="/cta/users/cta/users/merve.kaftancioglu/alternative_scenario/seqs/h3k27me3_ms"
output_directory="/cta/users/merve.kaftancioglu/fastQC_M2/job_outputs"

 Load FastQC module
module load fastqc-0.11.7-gcc-9.2.0-nsjzjof

 Run FastQC on .fastq files in input directory and write output to output directory
fastqc --outdir "$output_directory" "$input_directory"/*.fastq -threads 8

""" !/bin/bash
SBATCH --job-name=chipseq_pipeline
SBATCH --nodes=1
SBATCH --account=users
SBATCH --qos=cuda
SBATCH --partition=cuda
SBATCH --cpus-per-task=16
SBATCH -o /cta/users/merve.kaftancioglu/alternative_scenario/processed_seqs/out-err/%x-%j-slurm.out
SBATCH -e /cta/users/merve.kaftancioglu/alternative_scenario/processed_seqs/out-err/%x-%j-slurm.err

cd /cta/users/merve.kaftancioglu/alternative_scenario/nf-core/chipseq-master

nextflow run main.nf -profile singularity --outdir /output

""" !/bin/bash
SBATCH --job-name=trimmomatic
SBATCH --nodes=1
SBATCH --account=users
SBATCH --qos=short_mdbf
SBATCH --partition=short_mdbf
SBATCH --cpus-per-task=16
SBATCH -o /cta/users/merve.kaftancioglu/alternative_scenario/fastQC_as/out-err/%x-%j-slurm.out
SBATCH -e /cta/users/merve.kaftancioglu/alternative_scenario/fastQC_as/out-err/%x-%j-slurm.err

 Load FastQC module
module load fastqc-0.11.7-gcc-9.2.0-nsjzjof

substrate_ids=("h3k27ac_ms_b_ENCSR295QZX" "h3k27ac_ms_b_ENCSR364OIK" "h3k27ac_ms_b_ENCSR538URI" "h3k27ac_ms_b_ENCSR541PLO" "h3k27ac_ms_b_ENCSR617GMR" "h3k27ac_ms_cd4_ENCSR200SSJ" "h3k27ac_ms_cd4_ENCSR322MTA" "h3k27ac_ms_cd4_ENCSR331WMS" "h3k27ac_ms_cd4_ENCSR350UKV" "h3k27ac_ms_cd4_ENCSR474PYR" "h3k27ac_ms_cd4_ENCSR520QDR" "h3k27ac_ms_cd4_ENCSR540XNK" "h3k27ac_ms_cd4_ENCSR705VSO" "h3k27ac_ms_cd4_ENCSR832UMM" "h3k27ac_ms_cd8_ENCSR078ATS" "h3k27ac_ms_cd8_ENCSR348YRH" "h3k27ac_ms_cd8_ENCSR458TOW" "h3k27ac_ms_cd8_ENCSR476IPR" "h3k27ac_ms_cd8_ENCSR787HDF" "h3k27ac_ms_cd8_ENCSR923JIB" "h3k27ac_ms_nk_ENCSR246TTM" "h3k27ac_ms_nk_ENCSR469CFL" "h3k27ac_ms_nk_ENCSR746AIX" "h3k27ac_normal_b_ENCSR685KZA" "h3k27ac_normal_cd4_ENCSR819NCZ" "h3k27ac_normal_cd8_ENCSR976RWL" "h3k27ac_normal_nk_ENCSR977FMZ" "h3k27me3_ms_b_ENCSR009ZRH" "h3k27me3_ms_b_ENCSR182NLA" "h3k27me3_ms_b_ENCSR272YVX" "h3k27me3_ms_b_ENCSR649FUX" "h3k27me3_ms_b_ENCSR842WWX" "h3k27me3_ms_cd4_ENCSR277XYX" "h3k27me3_ms_cd4_ENCSR526TNC" "h3k27me3_ms_cd4_ENCSR592EKF" "h3k27me3_ms_cd4_ENCSR613UFD" "h3k27me3_ms_cd4_ENCSR740SDR" "h3k27me3_ms_cd4_ENCSR779JLY" "h3k27me3_ms_cd4_ENCSR993CTA" "h3k27me3_ms_cd8_ENCSR116FVG" "h3k27me3_ms_cd8_ENCSR122JCM" "h3k27me3_ms_cd8_ENCSR216ZVA" "h3k27me3_ms_cd8_ENCSR284IKS" "h3k27me3_ms_cd8_ENCSR385BOZ" "h3k27me3_ms_cd8_ENCSR521SFR" "h3k27me3_ms_nk_ENCSR469QVG" "h3k27me3_ms_nk_ENCSR565WDW" "h3k27me3_normal_b_ENCSR589LHR" "h3k27me3_normal_cd4_ENCSR068YVZ" "h3k27me3_normal_cd8_ENCSR720QRX" "h3k27me3_normal_nk_ENCSR639NIG" "h3k36me3_ms_b_ENCSR089VQL" "h3k36me3_ms_b_ENCSR119PSR" "h3k36me3_ms_b_ENCSR238WFK" "h3k36me3_ms_b_ENCSR987OPY" "h3k36me3_ms_cd4_ENCSR276NGH" "h3k36me3_ms_cd4_ENCSR330CQU" "h3k36me3_ms_cd4_ENCSR471WZE" "h3k36me3_ms_cd4_ENCSR473TGK" "h3k36me3_ms_cd4_ENCSR482VIB" "h3k36me3_ms_cd4_ENCSR532PXR" "h3k36me3_ms_cd4_ENCSR785PKM" "h3k36me3_ms_cd4_ENCSR865FTW" "h3k36me3_ms_cd4_ENCSR898VJE" "h3k36me3_ms_cd8_ENCSR239OWL" "h3k36me3_ms_cd8_ENCSR435MSB" "h3k36me3_ms_cd8_ENCSR631VIW" "h3k36me3_ms_cd8_ENCSR652WFO" "h3k36me3_ms_cd8_ENCSR757FGN" "h3k36me3_ms_nk_ENCSR158VSE" "h3k36me3_ms_nk_ENCSR245KON" "h3k36me3_ms_nk_ENCSR530YDY" "h3k36me3_normal_b_ENCSR831AXK" "h3k36me3_normal_cd4_ENCSR640OKB" "h3k36me3_normal_cd8_ENCSR303SQG" "h3k36me3_normal_nk_ENCSR056GJY" "h3k4me1_ms_b_ENCSR084FXT" "h3k4me1_ms_b_ENCSR480ITK" "h3k4me1_ms_b_ENCSR759DHN" "h3k4me1_ms_b_ENCSR900IAM" "h3k4me1_ms_b_ENCSR931JYC" "h3k4me1_ms_cd4_ENCSR036YVP" "h3k4me1_ms_cd4_ENCSR043JIA" "h3k4me1_ms_cd4_ENCSR093VUP" "h3k4me1_ms_cd4_ENCSR124IGC" "h3k4me1_ms_cd4_ENCSR315MYP" "h3k4me1_ms_cd4_ENCSR458QEY" "h3k4me1_ms_cd4_ENCSR485FBT" "h3k4me1_ms_cd4_ENCSR641ZFV" "h3k4me1_ms_cd4_ENCSR687CQX" "h3k4me1_ms_cd8_ENCSR231XAP" "h3k4me1_ms_cd8_ENCSR572XTB" "h3k4me1_ms_cd8_ENCSR788TEF" "h3k4me1_ms_cd8_ENCSR815YZL" "h3k4me1_ms_cd8_ENCSR861FEC" "h3k4me1_ms_cd8_ENCSR940PHE" "h3k4me1_ms_nk_ENCSR482UGV" "h3k4me1_ms_nk_ENCSR718LCI" "h3k4me1_ms_nk_ENCSR806HUT" "h3k4me1_normal_b_ENCSR156BXM" "h3k4me1_normal_cd4_ENCSR102SOR" "h3k4me1_normal_cd8_ENCSR217SHH" "h3k4me1_normal_nk_ENCSR277YKG" "h3k4me3_ms_b_ENCSR260CRI" "h3k4me3_ms_b_ENCSR461QMZ" "h3k4me3_ms_b_ENCSR530AVK" "h3k4me3_ms_b_ENCSR713ZYF" "h3k4me3_ms_b_ENCSR923EGO" "h3k4me3_ms_cd4_ENCSR180NCM" "h3k4me3_ms_cd4_ENCSR341QLC" "h3k4me3_ms_cd4_ENCSR482TGI" "h3k4me3_ms_cd4_ENCSR486XJK" "h3k4me3_ms_cd4_ENCSR496LKR" "h3k4me3_ms_cd4_ENCSR603LTN" "h3k4me3_ms_cd4_ENCSR802MXQ" "h3k4me3_ms_cd4_ENCSR878YHM" "h3k4me3_ms_cd4_ENCSR954ZLD" "h3k4me3_ms_cd8_ENCSR231ZZH" "h3k4me3_ms_cd8_ENCSR278QHR" "h3k4me3_ms_cd8_ENCSR516CKJ" "h3k4me3_ms_cd8_ENCSR535YYH" "h3k4me3_ms_cd8_ENCSR741XAE" "h3k4me3_ms_cd8_ENCSR848XJL" "h3k4me3_ms_nk_ENCSR234BAD" "h3k4me3_ms_nk_ENCSR703TWY" "h3k4me3_ms_nk_ENCSR755ZGS" "h3k4me3_normal_b_ENCSR791CAF" "h3k4me3_normal_cd4_ENCSR537KJA" "h3k4me3_normal_cd8_ENCSR123ZAT" "h3k4me3_normal_nk_ENCSR394JFQ" "h3k9me3_ms_b_ENCSR445DNH" "h3k9me3_ms_b_ENCSR486SZN" "h3k9me3_ms_b_ENCSR983MKX" "h3k9me3_ms_cd4_ENCSR057IZD" "h3k9me3_ms_cd4_ENCSR550DPT" "h3k9me3_ms_cd4_ENCSR677OEF" "h3k9me3_ms_cd4_ENCSR729ZXH" "h3k9me3_ms_cd4_ENCSR851RJV" "h3k9me3_ms_cd4_ENCSR919DFZ" "h3k9me3_ms_cd4_ENCSR953STZ" "h3k9me3_ms_cd4_ENCSR959VZU" "h3k9me3_ms_cd8_ENCSR101USF" "h3k9me3_ms_cd8_ENCSR354GNT" "h3k9me3_ms_cd8_ENCSR377QYB" "h3k9me3_ms_cd8_ENCSR733LCG" "h3k9me3_ms_cd8_ENCSR980KTW" "h3k9me3_ms_nk_ENCSR061ATV" "h3k9me3_ms_nk_ENCSR611YIA" "h3k9me3_ms_nk_ENCSR667HUI" "h3k9me3_normal_b_ENCSR445LTM" "h3k9me3_normal_cd4_ENCSR433EWI" "h3k9me3_normal_cd8_ENCSR294HTM" "h3k9me3_normal_nk_ENCSR025UNZ")

for substrate_id in "${substrate_ids[@]}"; do
    echo "Predicting ${substrate_id}..."
    input_directory="/cta/users/merve.kaftancioglu/alternative_scenario/seqs/${substrate_id}"
    output_dir="/cta/users/merve.kaftancioglu/alternative_scenario/fastQC_as/output/${substrate_id}"    
    mkdir -p "${output_dir}"
    trimmomatic PE -threads 16 -phred33 -trimlog "${output_dir}/trimmomatic_log.txt" "$input_directory"/R1.fastq.gz "$input_directory"/R2.fastq.gz "${output_dir}/output_forward_paired.fastq.gz" "${output_dir}/output_forward_unpaired.fastq.gz" "${output_dir}/output_reverse_paired.fastq.gz" "${output_dir}/output_reverse_unpaired.fastq.gz" ILLUMINACLIP:/path/to/adapters.fa:2:30:10 LEADING:3 TRAILING:3 SLIDINGWINDOW:4:15 MINLEN:36
    echo "Finished ${substrate_id}"
done


""" !/bin/bash
SBATCH --job-name=multiQC
SBATCH --nodes=1
SBATCH --account=users
SBATCH --qos=cuda
SBATCH --partition=cuda
SBATCH --cpus-per-task=16
SBATCH -o /cta/users/merve.kaftancioglu/alternative_scenario/post-trim_fastQC_as/output/out-err/%x-%j-slurm.out
SBATCH -e /cta/users/merve.kaftancioglu/alternative_scenario/post-trim_fastQC_as/output/out-err/%x-%j-slurm.err

 Specify the main directory containing the folders to merge
main_directory="/cta/users/merve.kaftancioglu/alternative_scenario/post-trim_fastQC_as/output_merged"

 Specify the directory where you want to merge the files
merge_directory="/cta/users/merve.kaftancioglu/alternative_scenario/post-trim_fastQC_as/output_merged"

 Copy all files from subdirectories to the merge directory
find "$main_directory" -type f -exec cp {} "$merge_directory" \;

echo "All files merged into: $merge_directory"

""" !/bin/bash
SBATCH --job-name=star_single
SBATCH --nodes=1
SBATCH --account=users
SBATCH --qos=mid_mdbf
SBATCH --partition=mid_mdbf
SBATCH --cpus-per-task=16
SBATCH -o /cta/users/merve.kaftancioglu/alternative_scenario/star_as/out-err/%x-%j-slurm.out
SBATCH -e /cta/users/merve.kaftancioglu/alternative_scenario/star_as/out-err/%x-%j-slurm.err

 Load STAR module
module load star-2.7.0e-gcc-9.2.0-vynasg3

 Define paths (replace with your actual paths)
genome_fasta="/cta/users/merve.kaftancioglu/alternative_scenario/star_as/reference/human_g1k_v37.fasta"
star_index_dir="/cta/users/merve.kaftancioglu/star_as/index"

 Check if genome fasta file exists
if [ ! -f "$genome_fasta" ]; then
  echo "Error: Genome fasta file not found at: $genome_fasta"
  exit 1
fi

 Check if star_index_dir exists (and create if not)
if [ ! -d "$star_index_dir" ]; then
  echo "Creating star_genome_index directory: $star_index_dir"
  mkdir -p "$star_index_dir"
fi

 Define number of threads (adjust based on your system)
threads=8

 STAR command
star_cmd="STAR \
  --runThreadN $threads \
  --runMode genomeGenerate \
  --genomeDir $star_index_dir \
  --genomeFastaFiles $genome_fasta"

 Run STAR and handle potential errors
echo "Building STAR genome index..."
$star_cmd

if [ $? -ne 0 ]; then
  echo "Error: STAR failed to generate genome index."
  exit 1
fi

echo "STAR genome index created successfully in: $star_index_dir"

""" !/bin/bash
SBATCH --job-name=nf-core
SBATCH --nodes=1
SBATCH --account=users
SBATCH --qos=cuda
SBATCH --partition=cuda
SBATCH --cpus-per-task=16
SBATCH -o /cta/users/merve.kaftancioglu/alternative_scenario/processed_seqs/out-err/%x-%j-slurm.out
SBATCH -e /cta/users/merve.kaftancioglu/alternative_scenario/processed_seqs/out-err/%x-%j-slurm.err

 Load necessary modules (if applicable for your cluster)
module load docker

 Specify the Docker image (replace with nf-core/chipseq version you downloaded)
docker_image=nf-core/chipseq:latest

 Navigate to the directory containing your Slurm script and pipeline files
cd /cta/users/merve.kaftancioglu/alternative_scenario/processed_seqs/all_data

 Create symbolic links for all fastq files in the current directory
for file in *.fastq; do ln -s "$file" all_fastq_files; done

 Run the nextflow command within the container
docker run --rm -v /cta/users/merve.kaftancioglu/alternative_scenario/trimmomatic_as/output

 Create a directory to hold symbolic links (soft links) to your fastq files (optional)
 This step is optional, you can comment it out if you don't want to create a separate directory
 mkdir all_fastq_files

 Loop through all subdirectories under processed_seqs
for sample_dir in */ ; do
   Navigate to the current subdirectory
  cd "$sample_dir"

   Loop through all fastq.gz files in the current subdirectory
  for file in *.fastq.gz; do
     Create a symbolic link (optional, uncomment if using the directory)
     ln -s "$file" ../all_fastq_files   Uncomment if using all_fastq_files directory
    
     Construct the full path to the fastq file (alternative to symbolic links)
    full_fastq_path=$(pwd)/"$file"
    
     You can echo the full path for verification (optional)
     echo "Found fastq: $full_fastq_path"
  done
  
   Move back to the processed_seqs directory after processing each subdirectory
  cd ..
done

""" !/bin/bash
SBATCH --job-name=picard
SBATCH --nodes=1
SBATCH --account=users
SBATCH --qos=cuda
SBATCH --partition=cuda
SBATCH --cpus-per-task=16
SBATCH -o /cta/users/merve.kaftancioglu/alternative_scenario/picard_as/out-err/%x-%j-slurm.out
SBATCH -e /cta/users/merve.kaftancioglu/alternative_scenario/picard_as/out-err/%x-%j-slurm.err

 Load BWA module

 Load required modules
module load picard samtools

 Define input and output file variables
MAIN_DIR="/cta/users/merve.kaftancioglu/alternative_scenario/samtools_as/output"   Replace with path to main directory
REALIGNED_BAMS_DIR="$MAIN_DIR/realigned_bams"   Output directory for realigned BAMs
TEMP_DIR="$MAIN_DIR/temp"   Optional temporary directory (create if needed)

 **Step 1: Mark duplicates on each individual BAM**

 Loop through subdirectories
for subdir in $MAIN_DIR/*; do
   Check if subdirectory is a directory (avoid hidden files)
  if [[ -d "$subdir" ]]; then
     Extract subdirectory name
    subdir_name=$(basename "$subdir")

     Find BAM files with wildcard
    bam_file="$subdir/$subdir_name*.bam"

     Check if BAM file exists
    if [[ -f "$bam_file" ]]; then
       Mark duplicates and place output in subdirectory
      mark_dups_cmd="picard MarkDuplicates \
        I=$bam_file \
        O=$subdir/$subdir_name.marked_dups.bam \
        M=$subdir/$subdir_name.marked_dups_metrics.txt \
        CREATE_INDEX=true"

      echo "Marking duplicates for $subdir_name..."
      $mark_dups_cmd
    else
      echo "No BAM file found in $subdir."
    fi
  fi
done

 **Step 2: Merge marked duplicate BAMs (optional)**

 List of temporary marked duplicate BAMs
marked_dups_bams="$TEMP_DIR/$SAMPLE_NAME*.marked_dups.bam"

 Merge marked duplicate BAMs
merge_bam_cmd="picard MergeSamFiles \
  O=$REALIGNED_BAMS_DIR/$SAMPLE_NAME.merged.marked_dups.bam \
  $marked_dups_bams"

echo "Merging marked duplicate BAMs..."
$merge_bam_cmd

 **Step 3: Re-mark duplicates on the merged BAM**

 Re-mark duplicates on the merged BAM
re_mark_dups_cmd="picard MarkDuplicates \
  I=$REALIGNED_BAMS_DIR/$SAMPLE_NAME.merged.marked_dups.bam \
  O=$REALIGNED_BAMS_DIR/$SAMPLE_NAME.realigned.bam \
  M=$REALIGNED_BAMS_DIR/$SAMPLE_NAME.realigned_metrics.txt \
  CREATE_INDEX=true"

echo "Re-marking duplicates on merged BAM..."
$re_mark_dups_cmd

 **Optional Step: Clean up temporary files (if used)**
 rm -rf $TEMP_DIR   Uncomment to remove temporary files

echo "DONE!""" "

 **Step 3: Re-mark duplicates on merged BAMs (optional)**

 If you performed merging, adjust commands to work with merged BAMs.

 **Optional Cleanup (if using temporary directory)**
 rm -rf $TEMP_DIR   Uncomment to remove temporary files

echo "DONE!""" "


""" !/bin/bash
SBATCH --job-name=bowtie2_alignment
SBATCH --nodes=1
SBATCH --account=users
SBATCH --qos=mid_mdbf
SBATCH --partition=mid_mdbf
SBATCH --cpus-per-task=16
SBATCH -o /cta/users/merve.kaftancioglu/alternative_scenario/bowtie2_as/out-err/%x-%j-slurm.out
SBATCH -e /cta/users/merve.kaftancioglu/alternative_scenario/bowtie2_as/out-err/%x-%j-slurm.err

 Load Bowtie2 module
module load bowtie2-2.3.5.1-gcc-9.2.0-ae4dozk

 Define directories
input_dir="/cta/users/merve.kaftancioglu/alternative_scenario/trimmomatic_as/output"
output_dir="/cta/users/merve.kaftancioglu/alternative_scenario/bowtie2_as/output"

 Path to reference genome FASTA file
reference_genome="/cta/users/merve.kaftancioglu/alternative_scenario/bowtie2_as/reference/human_g1k_v37.fasta"

 Create Bowtie2 index (if not already done)
bowtie2-build "$reference_genome" /cta/users/merve.kaftancioglu/alternative_scenario/bowtie2_as/index

 Loop through all FASTQ files in the input directory using find
for input_fastq in $(find "${input_dir}" -type f -name "*.fastq.gz" ! -name "*.*"); do

  sample_id="${input_fastq*/}"  
  sample_id="${sample_id//.fastq.gz}" 

  echo "Processing ${sample_id}..."

   Create output directory (if it doesn't exist)
  mkdir -p "${output_dir}/${sample_id}"

   Output SAM file path
  output_sam="${output_dir}/${sample_id}/${sample_id}.sam"

   Run Bowtie2 with the current FASTQ file, specifying output format (-U) as SAM
  bowtie2 --threads 16 -p 16 --trim3 30 -x /cta/users/merve.kaftancioglu/alternative_scenario/bowtie2_as/index -U $input_fastq -U $output_sam
   Note the change: -U $output_sam specifies the output SAM file
done

""" !/bin/bash
SBATCH --job-name=bowtie2_alignment
SBATCH --nodes=1
SBATCH --account=users
SBATCH --qos=mid_mdbf
SBATCH --partition=mid_mdbf
SBATCH --cpus-per-task=16
SBATCH -o /cta/users/merve.kaftancioglu/alternative_scenario/bowtie2_as/out-err/%x-%j-slurm.out
SBATCH -e /cta/users/merve.kaftancioglu/alternative_scenario/bowtie2_as/out-err/%x-%j-slurm.err

 Load Bowtie2 module
module load bowtie2-2.3.5.1-gcc-9.2.0-ae4dozk

 Define directories
input_dir="/cta/users/merve.kaftancioglu/alternative_scenario/trimmomatic_as/output"
output_dir="/cta/users/merve.kaftancioglu/alternative_scenario/bowtie2_as/output"

 Path to reference genome FASTA file
reference_genome="/cta/users/merve.kaftancioglu/alternative_scenario/bowtie2_as/reference/human_g1k_v37.fasta"

 Create Bowtie2 index (if not already done)
bowtie2-build "$reference_genome" /cta/users/merve.kaftancioglu/alternative_scenario/bowtie2_as/index

 Loop through all FASTQ files in the input directory
for input_fastq in "${input_dir}"/*.fastq.gz; do
   Extract sample ID from filename (assuming filenames start with sample ID)
   Modify this part to match your specific filename format if needed
  sample_id="${input_fastq*/}"   Double  removes everything before the last slash (/)
  sample_id="${sample_id//.fastq.gz}"   Remove .fastq.gz extension

  echo "Processing ${sample_id}..."

   Create output directory (if it doesn't exist)
  mkdir -p "${output_dir}/${sample_id}"

   Output BAM file path
  output_bam="${output_dir}/${sample_id}/${sample_id}.bam"

   Run Bowtie2 with the current FASTQ file
  bowtie2 -x /cta/users/merve.kaftancioglu/alternative_scenario/bowtie2_as/index -U "$input_fastq" -S "$output_bam"
done

""" !/bin/bash
SBATCH --job-name=bowtie2_alignment
SBATCH --nodes=1
SBATCH --account=users
SBATCH --qos=mid_mdbf
SBATCH --partition=mid_mdbf
SBATCH --cpus-per-task=16
SBATCH -o /cta/users/merve.kaftancioglu/alternative_scenario/bowtie2_as/out-err/%x-%j-slurm.out
SBATCH -e /cta/users/merve.kaftancioglu/alternative_scenario/bowtie2_as/out-err/%x-%j-slurm.err

 Load Bowtie2 module
module load bowtie2-2.3.5.1-gcc-9.2.0-ae4dozk

 Define directories
input_dir="/cta/users/merve.kaftancioglu/alternative_scenario/trimmomatic_as/output"
output_dir="/cta/users/merve.kaftancioglu/alternative_scenario/bowtie2_as/output"

 Path to reference genome FASTA file
reference_genome="/cta/users/merve.kaftancioglu/alternative_scenario/bowtie2_as/reference/human_g1k_v37.fasta"

 Create Bowtie2 index (if not already done)
bowtie2-build "$reference_genome" /cta/users/merve.kaftancioglu/alternative_scenario/bowtie2_as/index

 Loop through all FASTQ files in the input directory using find
for input_fastq in $(find "${input_dir}" -type f -name "*.fastq.gz" ! -name "*.*"); do

  sample_id="${input_fastq*/}"  
  sample_id="${sample_id//.fastq.gz}" 

  echo "Processing ${sample_id}..."

   Create output directory (if it doesn't exist)
  mkdir -p "${output_dir}/${sample_id}"

   Output BAM file path
  output_bam="${output_dir}/${sample_id}/${sample_id}.bam"

   Run Bowtie2 with the current FASTQ file
  bowtie2 --threads 16 -p 16 --trim3 30 -x reference_genome -U $input_fastq -S $output_bam
   bowtie2 -x /cta/users/merve.kaftancioglu/alternative_scenario/bowtie2_as/index -U "$input_fastq" -S "$output_bam"
done

""" !/bin/bash
SBATCH --job-name=bowtie2_alignment
SBATCH --nodes=1
SBATCH --account=users
SBATCH --qos=mid_mdbf
SBATCH --partition=mid_mdbf
SBATCH --cpus-per-task=16
SBATCH -o /cta/users/merve.kaftancioglu/alternative_scenario/bowtie2_as/out-err/%x-%j-slurm.out
SBATCH -e /cta/users/merve.kaftancioglu/alternative_scenario/bowtie2_as/out-err/%x-%j-slurm.err

 Load Bowtie2 module
module load bowtie2-2.3.5.1-gcc-9.2.0-ae4dozk

 Define directories
input_dir="/cta/users/merve.kaftancioglu/alternative_scenario/trimmomatic_as/output"
output_dir="/cta/users/merve.kaftancioglu/alternative_scenario/bowtie2_as/output"

 Path to reference genome FASTA file
reference_genome="/cta/users/merve.kaftancioglu/alternative_scenario/bowtie2_as/reference/human_g1k_v37.fasta"

 Create Bowtie2 index (if not already done)
bowtie2-build "$reference_genome" /cta/users/merve.kaftancioglu/alternative_scenario/bowtie2_as/index

 Loop through all FASTQ files in the input directory using find
for input_fastq in $(find "${input_dir}" -type f -name "*.fastq.gz" ! -name "*.*"); do

  sample_id="${input_fastq*/}"  
  sample_id="${sample_id//.fastq.gz}" 

  echo "Processing ${sample_id}..."

   Create output directory (if it doesn't exist)
  mkdir -p "${output_dir}/${sample_id}"

   Output SAM file path
  output_sam="${output_dir}/${sample_id}/${sample_id}.sam"

   Run Bowtie2 with the current FASTQ file, specifying output format (-S) as SAM
  bowtie2 --threads 16 -p 16 --trim3 30 -x /cta/users/merve.kaftancioglu/alternative_scenario/bowtie2_as/index -U $input_fastq -S $output_sam
done

""" !/bin/bash
SBATCH --job-name=bowtie2_longread_alignment
SBATCH --nodes=1
SBATCH --account=users
SBATCH --qos=mid_mdbf
SBATCH --partition=mid_mdbf
SBATCH --cpus-per-task=16
SBATCH -o /cta/users/merve.kaftancioglu/alternative_scenario/bowtie2_as/out-err/%x-%j-slurm.out
SBATCH -e /cta/users/merve.kaftancioglu/alternative_scenario/bowtie2_as/out-err/%x-%j-slurm.err

 Load Bowtie2 module
module load bowtie2-2.3.5.1-gcc-9.2.0-ae4dozk

 Define directories (consider adjusting paths)
input_dir="/cta/users/merve.kaftancioglu/alternative_scenario/trimmomatic_as/output"
output_dir="/cta/users/merve.kaftancioglu/alternative_scenario/bowtie2_as/longread_output"   Separate output directory for long reads

 Path to reference genome FASTA file
reference_genome="/cta/users/merve.kaftancioglu/alternative_scenario/bowtie2_as/reference/human_g1k_v37.fasta"

 Create Bowtie2 index (if not already done)
bowtie2-build "$reference_genome" /cta/users/merve.kaftancioglu/alternative_scenario/bowtie2_as/index

 Loop through all FASTQ files in the input directory
for input_fastq in $(find "${input_dir}" -type f -name "*.fastq.gz" ! -name "*.*"); do

  sample_id="${input_fastq*/}"  
  sample_id="${sample_id//.fastq.gz}" 

  echo "Processing ${sample_id}..."

   Create output directory (if it doesn't exist)
  mkdir -p "${output_dir}/${sample_id}"

   Output SAM file path
  output_sam="${output_dir}/${sample_id}/${sample_id}.sam"

   Run Bowtie2 with the current FASTQ file, specifying output format (-S) as SAM and "--very-sensitive" flag
  bowtie2 --threads 16 -p 16 --trim3 30 --very-sensitive -x /cta/users/merve.kaftancioglu/alternative_scenario/bowtie2_as/index -U $input_fastq -S $output_sam

   Optional: Check for successful completion (add this block if desired)
  if [ $? -eq 0 ]; then
    echo "Alignment for ${sample_id} completed successfully."
  else
    echo "Error occurred during alignment for ${sample_id}. Check slurm error logs for details."
  fi
done

@ECHO OFF

REM Get user inputs
SET /P pipeline_name=Enter the nf-core pipeline name (should be chipseq): 
IF NOT "%pipeline_name%"=="chipseq" (
  ECHO Error: Please enter "chipseq" as the pipeline name.
  GOTO EXIT
)

SET /P data_format=Enter data format (single-end or paired-end): 
SET /P input_path=/cta/users/merve.kaftancioglu/alternative_scenario/seqs: 
SET /P genome_path= /cta/users/merve.kaftancioglu/alternative_scenario/nf-core/reference/human_g1k_v37.fasta:
 SET /P blacklist_path=Enter the path to your blacklist BED file (optional): 
SET /P output_path=/cta/users/merve.kaftancioglu/alternative_scenario/nf-core/chipseq-2.0.0/output: 
SET /P profile=singularity:3.8.4  

REM Build the command with user inputs
SET command=nextflow run nf-core/chipseq -r 2.0.0 

REM Add options based on data format (only for paired-end)
IF "%data_format%"=="paired-end" (
  SET command=%command% --paired
)

SET command=%command% --input %input_path% --genome %genome_path% --outdir %output_path% -profile %profile%

REM Add optional blacklist path (uncomment if needed)
IF NOT "%blacklist_path%"=="" (
  SET command=%command% --blacklist %blacklist_path%
)

REM Echo the final command for verification
ECHO The command to be executed:
ECHO %command%

REM Pause before execution (optional)
PAUSE

REM Run the pipeline command
%command%

ECHO Pipeline execution finished!

PAUSE

:EXIT



""" !/bin/bash
SBATCH --job-name=bwa_single
SBATCH --nodes=1
SBATCH --account=users
SBATCH --qos=cuda
SBATCH --partition=cuda
SBATCH --cpus-per-task=16
SBATCH -o /cta/users/merve.kaftancioglu/alternative_scenario/processed_seqs/out-err/%x-%j-slurm.out
SBATCH -e /cta/users/merve.kaftancioglu/alternative_scenario/processed_seqs/out-err/%x-%j-slurm.err

 Load BWA module
module load rsync-3.1.3-gcc-9.2.0-mrwcim2

 Define directories (consider adjusting paths)
input_dir="/cta/users/merve.kaftancioglu/alternative_scenario/processed_seqs/all_data"
output_dir="/cta/users/merve.kaftancioglu/alternative_scenario/bwa_as/output"
reference_genome="/cta/users/merve.kaftancioglu/alternative_scenario/bwa_as/reference/human_g1k_v37.fasta"

 Check if BWA index exists
if [ ! -f "$reference_genome.bwt" ]; then
  echo "Error: BWA index not found. Please create the index using: bwa mem index $reference_genome"
  exit 1
fi

 Define sample IDs (modify as needed)
sample_ids=("h3k27ac_ms_b_ENCSR295QZX" "h3k27ac_ms_b_ENCSR538URI" "h3k27ac_ms_b_ENCSR617GMR" "h3k4me1_ms_b_ENCSR084FXT" "h3k4me1_ms_b_ENCSR480ITK" "h3k4me1_ms_b_ENCSR759DHN" "h3k4me1_ms_b_ENCSR900IAM" "h3k4me1_ms_b_ENCSR931JYC" "h3k4me3_ms_b_ENCSR260CRI" "h3k4me3_ms_b_ENCSR461QMZ" "h3k4me3_ms_b_ENCSR530AVK" "h3k4me3_ms_b_ENCSR713ZYF" "h3k4me3_ms_b_ENCSR923EGO" "h3k9me3_ms_b_ENCSR445DNH" "h3k9me3_ms_b_ENCSR486SZN" "h3k9me3_ms_b_ENCSR983MKX" "h3k27ac_ms_b_ENCSR364OIK" "h3k27ac_ms_b_ENCSR541PLO" "h3k36me3_ms_cd4_ENCSR471WZE" "h3k36me3_ms_cd4_ENCSR473TGK" "h3k4me1_ms_nk_ENCSR482UGV" "h3k4me1_ms_nk_ENCSR718LCI" "h3k4me1_ms_nk_ENCSR806HUT" "h3k4me3_ms_nk_ENCSR234BAD" "h3k4me3_ms_nk_ENCSR703TWY" "h3k4me3_ms_nk_ENCSR755ZGS" "h3k9me3_ms_nk_ENCSR061ATV" "h3k9me3_ms_nk_ENCSR611YIA" "h3k9me3_ms_nk_ENCSR667HUI" "h3k27ac_ms_nk_ENCSR246TTM")

 Loop through all sample subdirectories
for sample_dir in "$input_dir"/*; do
   Check if directory name matches any sample ID
  if [[ " ${sample_ids[@]} " =~ " $(basename "$sample_dir") " ]]; then
    echo "Processing sample: $(basename "$sample_dir")"

     Create output directory for the sample (if it doesn't exist)
    mkdir -p "$output_dir/$(basename "$sample_dir")"

     Find FASTQ files within the sample directory
     Adjust the pattern below if your file names differ
    fastq_files=( "$sample_dir"/*.fastq.gz )

     Check if any FASTQ files found
    if [[ ${fastq_files[@]} -eq 0 ]]; then
      echo "Warning: No FASTQ files found in $sample_dir"
      continue  Skip to next iteration if no FASTQ files
    fi

     BWA mem command with paired-end or single-end handling
    if [[ ${fastq_files[@]} -eq 2 ]]; then
       Paired-end reads (assuming *_1.fastq.gz and *_2.fastq.gz naming)
      bwa mem $reference_genome ${fastq_files[0]} ${fastq_files[1]} > "$output_dir/$(basename "$sample_dir").sam"
    else
       Single-end read (assuming only one FASTQ file)
       Extract only the sample name from the directory path
      sample_name=$(basename "$sample_dir" /)
      bwa mem $reference_genome ${fastq_files[0]} > "$output_dir/$sample_name.sam"
    fi
  fi
done

""" !/bin/bash
SBATCH --job-name=trimmomatic_single
SBATCH --nodes=1
SBATCH --account=users
SBATCH --qos=mid_mdbf
SBATCH --partition=mid_mdbf
SBATCH --cpus-per-task=16
SBATCH -o /cta/users/merve.kaftancioglu/alternative_scenario/trimmomatic_as/out-err/%x-%j-slurm.out
SBATCH -e /cta/users/merve.kaftancioglu/alternative_scenario/trimmomatic_as/out-err/%x-%j-slurm.err

 Load FastQC module
module load trimmomatic-0.38-gcc-9.2.0-3evdnxc

 List of sample IDs (folders)
sample_ids=("h3k27ac_ms_b_ENCSR541PLO" "h3k27ac_ms_b_ENCSR617GMR" "h3k27ac_ms_cd4_ENCSR200SSJ" "h3k27ac_ms_cd4_ENCSR322MTA" "h3k27ac_ms_cd4_ENCSR331WMS" "h3k27ac_ms_cd4_ENCSR350UKV" "h3k27ac_ms_cd4_ENCSR474PYR" "h3k27ac_ms_cd4_ENCSR520QDR" "h3k27ac_ms_cd4_ENCSR540XNK" "h3k27ac_ms_cd4_ENCSR705VSO" "h3k27ac_ms_cd4_ENCSR832UMM" "h3k27ac_ms_cd8_ENCSR078ATS" "h3k27ac_ms_cd8_ENCSR348YRH" "h3k27ac_ms_cd8_ENCSR458TOW" "h3k27ac_ms_cd8_ENCSR476IPR" "h3k27ac_ms_cd8_ENCSR787HDF" "h3k27ac_ms_cd8_ENCSR923JIB" "h3k27ac_ms_nk_ENCSR246TTM" "h3k27ac_ms_nk_ENCSR469CFL" "h3k27ac_ms_nk_ENCSR746AIX" "h3k27ac_normal_b_ENCSR685KZA" "h3k27ac_normal_cd4_ENCSR819NCZ" "h3k27ac_normal_cd8_ENCSR976RWL" "h3k27ac_normal_nk_ENCSR977FMZ" "h3k27me3_ms_b_ENCSR009ZRH" "h3k27me3_ms_b_ENCSR182NLA" "h3k27me3_ms_b_ENCSR272YVX" "h3k27me3_ms_b_ENCSR649FUX" "h3k27me3_ms_b_ENCSR842WWX" "h3k27me3_ms_cd4_ENCSR277XYX" "h3k27me3_ms_cd4_ENCSR526TNC" "h3k27me3_ms_cd4_ENCSR592EKF" "h3k27me3_ms_cd4_ENCSR613UFD" "h3k27me3_ms_cd4_ENCSR740SDR" "h3k27me3_ms_cd4_ENCSR779JLY" "h3k27me3_ms_cd4_ENCSR993CTA" "h3k27me3_ms_cd8_ENCSR116FVG" "h3k27me3_ms_cd8_ENCSR122JCM" "h3k27me3_ms_cd8_ENCSR216ZVA" "h3k27me3_ms_cd8_ENCSR284IKS" "h3k27me3_ms_cd8_ENCSR385BOZ" "h3k27me3_ms_cd8_ENCSR521SFR" "h3k27me3_ms_nk_ENCSR469QVG" "h3k27me3_ms_nk_ENCSR565WDW" "h3k27me3_normal_b_ENCSR589LHR" "h3k27me3_normal_cd4_ENCSR068YVZ" "h3k27me3_normal_cd8_ENCSR720QRX" "h3k27me3_normal_nk_ENCSR639NIG" "h3k36me3_ms_b_ENCSR089VQL" "h3k36me3_ms_b_ENCSR119PSR" "h3k36me3_ms_b_ENCSR238WFK" "h3k36me3_ms_b_ENCSR987OPY" "h3k36me3_ms_cd4_ENCSR276NGH" "h3k36me3_ms_cd4_ENCSR330CQU" "h3k36me3_ms_cd4_ENCSR471WZE" "h3k36me3_ms_cd4_ENCSR473TGK" "h3k36me3_ms_cd4_ENCSR482VIB" "h3k36me3_ms_cd4_ENCSR532PXR" "h3k36me3_ms_cd4_ENCSR785PKM" "h3k36me3_ms_cd4_ENCSR865FTW" "h3k36me3_ms_cd4_ENCSR898VJE" "h3k36me3_ms_cd8_ENCSR239OWL" "h3k36me3_ms_cd8_ENCSR435MSB" "h3k36me3_ms_cd8_ENCSR631VIW" "h3k36me3_ms_cd8_ENCSR652WFO" "h3k36me3_ms_cd8_ENCSR757FGN" "h3k36me3_ms_nk_ENCSR158VSE" "h3k36me3_ms_nk_ENCSR245KON" "h3k36me3_ms_nk_ENCSR530YDY" "h3k36me3_normal_b_ENCSR831AXK" "h3k36me3_normal_cd4_ENCSR640OKB" "h3k36me3_normal_cd8_ENCSR303SQG" "h3k36me3_normal_nk_ENCSR056GJY" "h3k4me1_ms_b_ENCSR084FXT" "h3k4me1_ms_b_ENCSR480ITK" "h3k4me1_ms_b_ENCSR759DHN" "h3k4me1_ms_b_ENCSR900IAM" "h3k4me1_ms_b_ENCSR931JYC" "h3k4me1_ms_cd4_ENCSR036YVP" "h3k4me1_ms_cd4_ENCSR043JIA" "h3k4me1_ms_cd4_ENCSR093VUP" "h3k4me1_ms_cd4_ENCSR124IGC" "h3k4me1_ms_cd4_ENCSR315MYP" "h3k4me1_ms_cd4_ENCSR458QEY" "h3k4me1_ms_cd4_ENCSR485FBT" "h3k4me1_ms_cd4_ENCSR641ZFV" "h3k4me1_ms_cd4_ENCSR687CQX" "h3k4me1_ms_cd8_ENCSR231XAP" "h3k4me1_ms_cd8_ENCSR572XTB" "h3k4me1_ms_cd8_ENCSR788TEF" "h3k4me1_ms_cd8_ENCSR815YZL" "h3k4me1_ms_cd8_ENCSR861FEC" "h3k4me1_ms_cd8_ENCSR940PHE" "h3k4me1_ms_nk_ENCSR482UGV" "h3k4me1_ms_nk_ENCSR718LCI" "h3k4me1_ms_nk_ENCSR806HUT" "h3k4me1_normal_b_ENCSR156BXM" "h3k4me1_normal_cd4_ENCSR102SOR" "h3k4me1_normal_cd8_ENCSR217SHH" "h3k4me1_normal_nk_ENCSR277YKG" "h3k4me3_ms_b_ENCSR260CRI" "h3k4me3_ms_b_ENCSR461QMZ" "h3k4me3_ms_b_ENCSR530AVK" "h3k4me3_ms_b_ENCSR713ZYF" "h3k4me3_ms_b_ENCSR923EGO" "h3k4me3_ms_cd4_ENCSR180NCM" "h3k4me3_ms_cd4_ENCSR341QLC" "h3k4me3_ms_cd4_ENCSR482TGI" "h3k4me3_ms_cd4_ENCSR486XJK" "h3k4me3_ms_cd4_ENCSR496LKR" "h3k4me3_ms_cd4_ENCSR603LTN" "h3k4me3_ms_cd4_ENCSR802MXQ" "h3k4me3_ms_cd4_ENCSR878YHM" "h3k4me3_ms_cd4_ENCSR954ZLD" "h3k4me3_ms_cd8_ENCSR231ZZH" "h3k4me3_ms_cd8_ENCSR278QHR" "h3k4me3_ms_cd8_ENCSR516CKJ" "h3k4me3_ms_cd8_ENCSR535YYH" "h3k4me3_ms_cd8_ENCSR741XAE" "h3k4me3_ms_cd8_ENCSR848XJL" "h3k4me3_ms_nk_ENCSR234BAD" "h3k4me3_ms_nk_ENCSR703TWY" "h3k4me3_ms_nk_ENCSR755ZGS" "h3k4me3_normal_b_ENCSR791CAF" "h3k4me3_normal_cd4_ENCSR537KJA" "h3k4me3_normal_cd8_ENCSR123ZAT" "h3k4me3_normal_nk_ENCSR394JFQ" "h3k9me3_ms_b_ENCSR445DNH" "h3k9me3_ms_b_ENCSR486SZN" "h3k9me3_ms_b_ENCSR983MKX" "h3k9me3_ms_cd4_ENCSR057IZD" "h3k9me3_ms_cd4_ENCSR550DPT" "h3k9me3_ms_cd4_ENCSR677OEF" "h3k9me3_ms_cd4_ENCSR729ZXH" "h3k9me3_ms_cd4_ENCSR851RJV" "h3k9me3_ms_cd4_ENCSR919DFZ" "h3k9me3_ms_cd4_ENCSR953STZ" "h3k9me3_ms_cd4_ENCSR959VZU" "h3k9me3_ms_cd8_ENCSR101USF" "h3k9me3_ms_cd8_ENCSR354GNT" "h3k9me3_ms_cd8_ENCSR377QYB" "h3k9me3_ms_cd8_ENCSR733LCG" "h3k9me3_ms_cd8_ENCSR980KTW" "h3k9me3_ms_nk_ENCSR061ATV" "h3k9me3_ms_nk_ENCSR611YIA" "h3k9me3_ms_nk_ENCSR667HUI" "h3k9me3_normal_b_ENCSR445LTM" "h3k9me3_normal_cd4_ENCSR433EWI" "h3k9me3_normal_cd8_ENCSR294HTM" "h3k9me3_normal_nk_ENCSR025UNZ")

 Process each sample folder
for sample_id in "${sample_ids[@]}"; do
    input_folder="/cta/users/merve.kaftancioglu/alternative_scenario/seqs/$sample_id"
    output_folder="/cta/users/merve.kaftancioglu/alternative_scenario/trimmomatic_as/output/$sample_id"

     Create output folder if it doesn't exist
    mkdir -p "$output_folder"

     Process all .fastq.gz files in the specified folder
    for input_file in "$input_folder"/*.fastq.gz; do
        output_file="$output_folder/output_paired_$(basename "$input_file")"
        trimmomatic SE -threads 16 -phred33 -trimlog "$output_folder/trimmomatic_log.txt" "$input_file" "$output_file" HEADCROP:14
    done
done

""" !/bin/bash
SBATCH --job-name=trimmomatic_single
SBATCH --nodes=1
SBATCH --account=users
SBATCH --qos=mid_mdbf
SBATCH --partition=mid_mdbf
SBATCH --cpus-per-task=16
SBATCH -o /cta/users/merve.kaftancioglu/alternative_scenario/trimmomatic_as/out-err/%x-%j-slurm.out
SBATCH -e /cta/users/merve.kaftancioglu/alternative_scenario/trimmomatic_as/out-err/%x-%j-slurm.err

 Load FastQC module
module load trimmomatic-0.38-gcc-9.2.0-3evdnxc

 List of sample IDs (folders)
sample_ids=("h3k27ac_ms_cd4_ENCSR331WMS" "h3k27ac_ms_cd4_ENCSR350UKV" "h3k27ac_ms_cd4_ENCSR474PYR")

 Process each sample folder
for sample_id in "${sample_ids[@]}"; do
    input_folder="/cta/users/merve.kaftancioglu/alternative_scenario/seqs/$sample_id"
    output_folder="/cta/users/merve.kaftancioglu/alternative_scenario/trimmomatic_as/output/$sample_id"

     Create output folder if it doesn't exist
    mkdir -p "$output_folder"

     Process all .fastq.gz files in the specified folder
    for input_file in "$input_folder"/*.fastq.gz; do
        output_file="$output_folder/output_paired_$(basename "$input_file")"
        trimmomatic SE -threads 16 -phred33 -trimlog "$output_folder/trimmomatic_log.txt" "$input_file" "$output_file" HEADCROP:14
    done
done

""" !/bin/bash
SBATCH --job-name=trimmomatic_single
SBATCH --nodes=1
SBATCH --account=users
SBATCH --qos=mid_mdbf
SBATCH --partition=mid_mdbf
SBATCH --cpus-per-task=16
SBATCH -o /cta/users/merve.kaftancioglu/alternative_scenario/trimmomatic_as/out-err/%x-%j-slurm.out
SBATCH -e /cta/users/merve.kaftancioglu/alternative_scenario/trimmomatic_as/out-err/%x-%j-slurm.err

 Load FastQC module
module load trimmomatic-0.38-gcc-9.2.0-3evdnxc

 List of sample IDs (folders)
sample_ids=("h3k27ac_ms_cd4_ENCSR520QDR" "h3k27ac_ms_cd4_ENCSR540XNK" "h3k27ac_ms_cd4_ENCSR705VSO")

 Process each sample folder
for sample_id in "${sample_ids[@]}"; do
    input_folder="/cta/users/merve.kaftancioglu/alternative_scenario/seqs/$sample_id"
    output_folder="/cta/users/merve.kaftancioglu/alternative_scenario/trimmomatic_as/output/$sample_id"

     Create output folder if it doesn't exist
    mkdir -p "$output_folder"

     Process all .fastq.gz files in the specified folder
    for input_file in "$input_folder"/*.fastq.gz; do
        output_file="$output_folder/output_paired_$(basename "$input_file")"
        trimmomatic SE -threads 16 -phred33 -trimlog "$output_folder/trimmomatic_log.txt" "$input_file" "$output_file" HEADCROP:14
    done
done

""" !/bin/bash
SBATCH --job-name=trimmomatic_single
SBATCH --nodes=1
SBATCH --account=users
SBATCH --qos=mid_mdbf
SBATCH --partition=mid_mdbf
SBATCH --cpus-per-task=16
SBATCH -o /cta/users/merve.kaftancioglu/alternative_scenario/trimmomatic_as/out-err/%x-%j-slurm.out
SBATCH -e /cta/users/merve.kaftancioglu/alternative_scenario/trimmomatic_as/out-err/%x-%j-slurm.err

 Load FastQC module
module load trimmomatic-0.38-gcc-9.2.0-3evdnxc

 List of sample IDs (folders)
sample_ids=("h3k27ac_ms_cd8_ENCSR458TOW" "h3k27ac_ms_cd8_ENCSR476IPR" "h3k27ac_ms_cd8_ENCSR787HDF")

 Process each sample folder
for sample_id in "${sample_ids[@]}"; do
    input_folder="/cta/users/merve.kaftancioglu/alternative_scenario/seqs/$sample_id"
    output_folder="/cta/users/merve.kaftancioglu/alternative_scenario/trimmomatic_as/output/$sample_id"

     Create output folder if it doesn't exist
    mkdir -p "$output_folder"

     Process all .fastq.gz files in the specified folder
    for input_file in "$input_folder"/*.fastq.gz; do
        output_file="$output_folder/output_paired_$(basename "$input_file")"
        trimmomatic SE -threads 16 -phred33 -trimlog "$output_folder/trimmomatic_log.txt" "$input_file" "$output_file" HEADCROP:14
    done
done

""" !/bin/bash
SBATCH --job-name=trimmomatic_single
SBATCH --nodes=1
SBATCH --account=users
SBATCH --qos=mid_mdbf
SBATCH --partition=mid_mdbf
SBATCH --cpus-per-task=16
SBATCH -o /cta/users/merve.kaftancioglu/alternative_scenario/trimmomatic_as/out-err/%x-%j-slurm.out
SBATCH -e /cta/users/merve.kaftancioglu/alternative_scenario/trimmomatic_as/out-err/%x-%j-slurm.err

 Load FastQC module
module load trimmomatic-0.38-gcc-9.2.0-3evdnxc

 List of sample IDs (folders)
sample_ids=("h3k27ac_ms_cd8_ENCSR923JIB" "h3k27ac_ms_nk_ENCSR246TTM" "h3k27ac_ms_nk_ENCSR469CFL")

 Process each sample folder
for sample_id in "${sample_ids[@]}"; do
    input_folder="/cta/users/merve.kaftancioglu/alternative_scenario/seqs/$sample_id"
    output_folder="/cta/users/merve.kaftancioglu/alternative_scenario/trimmomatic_as/output/$sample_id"

     Create output folder if it doesn't exist
    mkdir -p "$output_folder"

     Process all .fastq.gz files in the specified folder
    for input_file in "$input_folder"/*.fastq.gz; do
        output_file="$output_folder/output_paired_$(basename "$input_file")"
        trimmomatic SE -threads 16 -phred33 -trimlog "$output_folder/trimmomatic_log.txt" "$input_file" "$output_file" HEADCROP:14
    done
done

""" !/bin/bash
SBATCH --job-name=trimmomatic_single
SBATCH --nodes=1
SBATCH --account=users
SBATCH --qos=mid_mdbf
SBATCH --partition=mid_mdbf
SBATCH --cpus-per-task=16
SBATCH -o /cta/users/merve.kaftancioglu/alternative_scenario/trimmomatic_as/out-err/%x-%j-slurm.out
SBATCH -e /cta/users/merve.kaftancioglu/alternative_scenario/trimmomatic_as/out-err/%x-%j-slurm.err

 Load FastQC module
module load trimmomatic-0.38-gcc-9.2.0-3evdnxc

 List of sample IDs (folders)
sample_ids=("h3k27ac_ms_nk_ENCSR746AIX" "h3k27ac_normal_b_ENCSR685KZA" "h3k27ac_normal_cd4_ENCSR819NCZ")

 Process each sample folder
for sample_id in "${sample_ids[@]}"; do
    input_folder="/cta/users/merve.kaftancioglu/alternative_scenario/seqs/$sample_id"
    output_folder="/cta/users/merve.kaftancioglu/alternative_scenario/trimmomatic_as/output/$sample_id"

     Create output folder if it doesn't exist
    mkdir -p "$output_folder"

     Process all .fastq.gz files in the specified folder
    for input_file in "$input_folder"/*.fastq.gz; do
        output_file="$output_folder/output_paired_$(basename "$input_file")"
        trimmomatic SE -threads 16 -phred33 -trimlog "$output_folder/trimmomatic_log.txt" "$input_file" "$output_file" HEADCROP:14
    done
done

""" !/bin/bash
SBATCH --job-name=trimmomatic_single
SBATCH --nodes=1
SBATCH --account=users
SBATCH --qos=mid_mdbf
SBATCH --partition=mid_mdbf
SBATCH --cpus-per-task=16
SBATCH -o /cta/users/merve.kaftancioglu/alternative_scenario/trimmomatic_as/out-err/%x-%j-slurm.out
SBATCH -e /cta/users/merve.kaftancioglu/alternative_scenario/trimmomatic_as/out-err/%x-%j-slurm.err

 Load FastQC module
module load trimmomatic-0.38-gcc-9.2.0-3evdnxc

 List of sample IDs (folders)
sample_ids=("h3k27ac_normal_cd8_ENCSR976RWL" "h3k27ac_normal_nk_ENCSR977FMZ" "h3k27me3_ms_b_ENCSR009ZRH")

 Process each sample folder
for sample_id in "${sample_ids[@]}"; do
    input_folder="/cta/users/merve.kaftancioglu/alternative_scenario/seqs/$sample_id"
    output_folder="/cta/users/merve.kaftancioglu/alternative_scenario/trimmomatic_as/output/$sample_id"

     Create output folder if it doesn't exist
    mkdir -p "$output_folder"

     Process all .fastq.gz files in the specified folder
    for input_file in "$input_folder"/*.fastq.gz; do
        output_file="$output_folder/output_paired_$(basename "$input_file")"
        trimmomatic SE -threads 16 -phred33 -trimlog "$output_folder/trimmomatic_log.txt" "$input_file" "$output_file" HEADCROP:14
    done
done

""" !/bin/bash
SBATCH --job-name=trimmomatic_single
SBATCH --nodes=1
SBATCH --account=users
SBATCH --qos=mid_mdbf
SBATCH --partition=mid_mdbf
SBATCH --cpus-per-task=16
SBATCH -o /cta/users/merve.kaftancioglu/alternative_scenario/trimmomatic_as/out-err/%x-%j-slurm.out
SBATCH -e /cta/users/merve.kaftancioglu/alternative_scenario/trimmomatic_as/out-err/%x-%j-slurm.err

 Load FastQC module
module load trimmomatic-0.38-gcc-9.2.0-3evdnxc

 List of sample IDs (folders)
sample_ids=("h3k27me3_ms_b_ENCSR182NLA" "h3k27me3_ms_b_ENCSR272YVX" "h3k27me3_ms_b_ENCSR649FUX")

 Process each sample folder
for sample_id in "${sample_ids[@]}"; do
    input_folder="/cta/users/merve.kaftancioglu/alternative_scenario/seqs/$sample_id"
    output_folder="/cta/users/merve.kaftancioglu/alternative_scenario/trimmomatic_as/output/$sample_id"

     Create output folder if it doesn't exist
    mkdir -p "$output_folder"

     Process all .fastq.gz files in the specified folder
    for input_file in "$input_folder"/*.fastq.gz; do
        output_file="$output_folder/output_paired_$(basename "$input_file")"
        trimmomatic SE -threads 16 -phred33 -trimlog "$output_folder/trimmomatic_log.txt" "$input_file" "$output_file" HEADCROP:14
    done
done

""" !/bin/bash
SBATCH --job-name=fastQC
SBATCH --nodes=1
SBATCH --account=users
SBATCH --qos=cuda
SBATCH --partition=cuda
SBATCH --cpus-per-task=16
SBATCH -o /cta/users/merve.kaftancioglu/alternative_scenario/fastQC_as/out-err/%x-%j-slurm.out
SBATCH -e /cta/users/merve.kaftancioglu/alternative_scenario/fastQC_as/out-err/%x-%j-slurm.err


 Load FastQC module
module load fastqc-0.11.7-gcc-9.2.0-nsjzjof

substrate_ids=("h3k27ac_ms_b_ENCSR295QZX" "h3k27ac_ms_b_ENCSR364OIK" "h3k27ac_ms_b_ENCSR538URI" "h3k27ac_ms_b_ENCSR541PLO" "h3k27ac_ms_b_ENCSR617GMR" "h3k27ac_ms_cd4_ENCSR200SSJ" "h3k27ac_ms_cd4_ENCSR322MTA" "h3k27ac_ms_cd4_ENCSR331WMS" "h3k27ac_ms_cd4_ENCSR350UKV" "h3k27ac_ms_cd4_ENCSR474PYR" "h3k27ac_ms_cd4_ENCSR520QDR" "h3k27ac_ms_cd4_ENCSR540XNK" "h3k27ac_ms_cd4_ENCSR705VSO" "h3k27ac_ms_cd4_ENCSR832UMM" "h3k27ac_ms_cd8_ENCSR078ATS" "h3k27ac_ms_cd8_ENCSR348YRH" "h3k27ac_ms_cd8_ENCSR458TOW" "h3k27ac_ms_cd8_ENCSR476IPR" "h3k27ac_ms_cd8_ENCSR787HDF" "h3k27ac_ms_cd8_ENCSR923JIB" "h3k27ac_ms_nk_ENCSR246TTM" "h3k27ac_ms_nk_ENCSR469CFL" "h3k27ac_ms_nk_ENCSR746AIX" "h3k27ac_normal_b_ENCSR685KZA" "h3k27ac_normal_cd4_ENCSR819NCZ" "h3k27ac_normal_cd8_ENCSR976RWL" "h3k27ac_normal_nk_ENCSR977FMZ" "h3k27me3_ms_b_ENCSR009ZRH" "h3k27me3_ms_b_ENCSR182NLA" "h3k27me3_ms_b_ENCSR272YVX" "h3k27me3_ms_b_ENCSR649FUX" "h3k27me3_ms_b_ENCSR842WWX" "h3k27me3_ms_cd4_ENCSR277XYX" "h3k27me3_ms_cd4_ENCSR526TNC" "h3k27me3_ms_cd4_ENCSR592EKF" "h3k27me3_ms_cd4_ENCSR613UFD" "h3k27me3_ms_cd4_ENCSR740SDR" "h3k27me3_ms_cd4_ENCSR779JLY" "h3k27me3_ms_cd4_ENCSR993CTA" "h3k27me3_ms_cd8_ENCSR116FVG" "h3k27me3_ms_cd8_ENCSR122JCM" "h3k27me3_ms_cd8_ENCSR216ZVA" "h3k27me3_ms_cd8_ENCSR284IKS" "h3k27me3_ms_cd8_ENCSR385BOZ" "h3k27me3_ms_cd8_ENCSR521SFR" "h3k27me3_ms_nk_ENCSR469QVG" "h3k27me3_ms_nk_ENCSR565WDW" "h3k27me3_normal_b_ENCSR589LHR" "h3k27me3_normal_cd4_ENCSR068YVZ" "h3k27me3_normal_cd8_ENCSR720QRX" "h3k27me3_normal_nk_ENCSR639NIG" "h3k36me3_ms_b_ENCSR089VQL" "h3k36me3_ms_b_ENCSR119PSR" "h3k36me3_ms_b_ENCSR238WFK" "h3k36me3_ms_b_ENCSR987OPY" "h3k36me3_ms_cd4_ENCSR276NGH" "h3k36me3_ms_cd4_ENCSR330CQU" "h3k36me3_ms_cd4_ENCSR471WZE" "h3k36me3_ms_cd4_ENCSR473TGK" "h3k36me3_ms_cd4_ENCSR482VIB" "h3k36me3_ms_cd4_ENCSR532PXR" "h3k36me3_ms_cd4_ENCSR785PKM" "h3k36me3_ms_cd4_ENCSR865FTW" "h3k36me3_ms_cd4_ENCSR898VJE" "h3k36me3_ms_cd8_ENCSR239OWL" "h3k36me3_ms_cd8_ENCSR435MSB" "h3k36me3_ms_cd8_ENCSR631VIW" "h3k36me3_ms_cd8_ENCSR652WFO" "h3k36me3_ms_cd8_ENCSR757FGN" "h3k36me3_ms_nk_ENCSR158VSE" "h3k36me3_ms_nk_ENCSR245KON" "h3k36me3_ms_nk_ENCSR530YDY" "h3k36me3_normal_b_ENCSR831AXK" "h3k36me3_normal_cd4_ENCSR640OKB" "h3k36me3_normal_cd8_ENCSR303SQG" "h3k36me3_normal_nk_ENCSR056GJY" "h3k4me1_ms_b_ENCSR084FXT" "h3k4me1_ms_b_ENCSR480ITK" "h3k4me1_ms_b_ENCSR759DHN" "h3k4me1_ms_b_ENCSR900IAM" "h3k4me1_ms_b_ENCSR931JYC" "h3k4me1_ms_cd4_ENCSR036YVP" "h3k4me1_ms_cd4_ENCSR043JIA" "h3k4me1_ms_cd4_ENCSR093VUP" "h3k4me1_ms_cd4_ENCSR124IGC" "h3k4me1_ms_cd4_ENCSR315MYP" "h3k4me1_ms_cd4_ENCSR458QEY" "h3k4me1_ms_cd4_ENCSR485FBT" "h3k4me1_ms_cd4_ENCSR641ZFV" "h3k4me1_ms_cd4_ENCSR687CQX" "h3k4me1_ms_cd8_ENCSR231XAP" "h3k4me1_ms_cd8_ENCSR572XTB" "h3k4me1_ms_cd8_ENCSR788TEF" "h3k4me1_ms_cd8_ENCSR815YZL" "h3k4me1_ms_cd8_ENCSR861FEC" "h3k4me1_ms_cd8_ENCSR940PHE" "h3k4me1_ms_nk_ENCSR482UGV" "h3k4me1_ms_nk_ENCSR718LCI" "h3k4me1_ms_nk_ENCSR806HUT" "h3k4me1_normal_b_ENCSR156BXM" "h3k4me1_normal_cd4_ENCSR102SOR" "h3k4me1_normal_cd8_ENCSR217SHH" "h3k4me1_normal_nk_ENCSR277YKG" "h3k4me3_ms_b_ENCSR260CRI" "h3k4me3_ms_b_ENCSR461QMZ" "h3k4me3_ms_b_ENCSR530AVK" "h3k4me3_ms_b_ENCSR713ZYF" "h3k4me3_ms_b_ENCSR923EGO" "h3k4me3_ms_cd4_ENCSR180NCM" "h3k4me3_ms_cd4_ENCSR341QLC" "h3k4me3_ms_cd4_ENCSR482TGI" "h3k4me3_ms_cd4_ENCSR486XJK" "h3k4me3_ms_cd4_ENCSR496LKR" "h3k4me3_ms_cd4_ENCSR603LTN" "h3k4me3_ms_cd4_ENCSR802MXQ" "h3k4me3_ms_cd4_ENCSR878YHM" "h3k4me3_ms_cd4_ENCSR954ZLD" "h3k4me3_ms_cd8_ENCSR231ZZH" "h3k4me3_ms_cd8_ENCSR278QHR" "h3k4me3_ms_cd8_ENCSR516CKJ" "h3k4me3_ms_cd8_ENCSR535YYH" "h3k4me3_ms_cd8_ENCSR741XAE" "h3k4me3_ms_cd8_ENCSR848XJL" "h3k4me3_ms_nk_ENCSR234BAD" "h3k4me3_ms_nk_ENCSR703TWY" "h3k4me3_ms_nk_ENCSR755ZGS" "h3k4me3_normal_b_ENCSR791CAF" "h3k4me3_normal_cd4_ENCSR537KJA" "h3k4me3_normal_cd8_ENCSR123ZAT" "h3k4me3_normal_nk_ENCSR394JFQ" "h3k9me3_ms_b_ENCSR445DNH" "h3k9me3_ms_b_ENCSR486SZN" "h3k9me3_ms_b_ENCSR983MKX" "h3k9me3_ms_cd4_ENCSR057IZD" "h3k9me3_ms_cd4_ENCSR550DPT" "h3k9me3_ms_cd4_ENCSR677OEF" "h3k9me3_ms_cd4_ENCSR729ZXH" "h3k9me3_ms_cd4_ENCSR851RJV" "h3k9me3_ms_cd4_ENCSR919DFZ" "h3k9me3_ms_cd4_ENCSR953STZ" "h3k9me3_ms_cd4_ENCSR959VZU" "h3k9me3_ms_cd8_ENCSR101USF" "h3k9me3_ms_cd8_ENCSR354GNT" "h3k9me3_ms_cd8_ENCSR377QYB" "h3k9me3_ms_cd8_ENCSR733LCG" "h3k9me3_ms_cd8_ENCSR980KTW" "h3k9me3_ms_nk_ENCSR061ATV" "h3k9me3_ms_nk_ENCSR611YIA" "h3k9me3_ms_nk_ENCSR667HUI" "h3k9me3_normal_b_ENCSR445LTM" "h3k9me3_normal_cd4_ENCSR433EWI" "h3k9me3_normal_cd8_ENCSR294HTM" "h3k9me3_normal_nk_ENCSR025UNZ")

for substrate_id in "${substrate_ids[@]}"; do
    echo "Predicting ${substrate_id}..."
    input_directory="/cta/users/merve.kaftancioglu/alternative_scenario/seqs/${substrate_id}"
    output_dir="/cta/users/merve.kaftancioglu/alternative_scenario/fastQC_as/output/${substrate_id}"    
    mkdir -p "${output_dir}"
    fastqc --outdir "$output_directory" "$input_directory"/*.fastq.gz -threads 16
    echo "Finished ${substrate_id}"
done


echo "DONE!""" "



""" !/bin/bash
SBATCH --job-name=conversion
SBATCH --nodes=1
SBATCH --account=users
SBATCH --qos=cuda
SBATCH --partition=cuda
SBATCH --cpus-per-task=16
SBATCH -o /cta/users/merve.kaftancioglu/alternative_scenario/bowtie2_as/output/out-err/%x-%j-slurm.out
SBATCH -e /cta/users/merve.kaftancioglu/alternative_scenario/bowtie2_as/output/out-err/%x-%j-slurm.err

 Set directory containing fastq files
fastq_dir="/cta/users/merve.kaftancioglu/alternative_scenario/seqs/h3k27ac_normal_cd4_ENCSR819NCZ"

 Set output directory (optional, use same directory by default)
fasta_dir="/cta/users/merve.kaftancioglu/alternative_scenario/seqs/h3k27ac_normal_cd4_ENCSR819NCZ_fasta"   Uncomment if you want a separate directory

 Loop through all fastq files in the directory
for fastq_file in "$fastq_dir"/*.fastq.gz; do
   Extract filename without extension
  fasta_filename="${fastq_file%.fastq.gz}.fasta"

   Convert fastq.gz to fasta using seqtk
  seqtk seq -a "$fastq_file" > "$fasta_filename"

  echo "Converted $fastq_file to $fasta_filename"
done

echo "All fastq files converted to fasta!"

""" !/bin/bash
SBATCH --job-name=trimmomatic_single
SBATCH --nodes=1
SBATCH --account=users
SBATCH --qos=mid_mdbf
SBATCH --partition=mid_mdbf
SBATCH --cpus-per-task=16
SBATCH -o /cta/users/merve.kaftancioglu/alternative_scenario/trimmomatic_as/out-err/%x-%j-slurm.out
SBATCH -e /cta/users/merve.kaftancioglu/alternative_scenario/trimmomatic_as/out-err/%x-%j-slurm.err

 Load FastQC module
module load trimmomatic-0.38-gcc-9.2.0-3evdnxc

 List of sample IDs (folders)
sample_ids=("h3k27me3_ms_b_ENCSR842WWX" "h3k27me3_ms_cd4_ENCSR277XYX" "h3k27me3_ms_cd4_ENCSR526TNC")

 Process each sample folder
for sample_id in "${sample_ids[@]}"; do
    input_folder="/cta/users/merve.kaftancioglu/alternative_scenario/seqs/$sample_id"
    output_folder="/cta/users/merve.kaftancioglu/alternative_scenario/trimmomatic_as/output/$sample_id"

     Create output folder if it doesn't exist
    mkdir -p "$output_folder"

     Process all .fastq.gz files in the specified folder
    for input_file in "$input_folder"/*.fastq.gz; do
        output_file="$output_folder/output_paired_$(basename "$input_file")"
        trimmomatic SE -threads 16 -phred33 -trimlog "$output_folder/trimmomatic_log.txt" "$input_file" "$output_file" HEADCROP:14
    done
done

""" !/bin/bash
SBATCH --job-name=trimmomatic_single
SBATCH --nodes=1
SBATCH --account=users
SBATCH --qos=mid_mdbf
SBATCH --partition=mid_mdbf
SBATCH --cpus-per-task=16
SBATCH -o /cta/users/merve.kaftancioglu/alternative_scenario/trimmomatic_as/out-err/%x-%j-slurm.out
SBATCH -e /cta/users/merve.kaftancioglu/alternative_scenario/trimmomatic_as/out-err/%x-%j-slurm.err

 Load FastQC module
module load trimmomatic-0.38-gcc-9.2.0-3evdnxc

 List of sample IDs (folders)
sample_ids=("h3k27me3_ms_cd4_ENCSR592EKF" "h3k27me3_ms_cd4_ENCSR613UFD" "h3k27me3_ms_cd4_ENCSR740SDR")

 Process each sample folder
for sample_id in "${sample_ids[@]}"; do
    input_folder="/cta/users/merve.kaftancioglu/alternative_scenario/seqs/$sample_id"
    output_folder="/cta/users/merve.kaftancioglu/alternative_scenario/trimmomatic_as/output/$sample_id"

     Create output folder if it doesn't exist
    mkdir -p "$output_folder"

     Process all .fastq.gz files in the specified folder
    for input_file in "$input_folder"/*.fastq.gz; do
        output_file="$output_folder/output_paired_$(basename "$input_file")"
        trimmomatic SE -threads 16 -phred33 -trimlog "$output_folder/trimmomatic_log.txt" "$input_file" "$output_file" HEADCROP:14
    done
done

""" !/bin/bash
SBATCH --job-name=trimmomatic_single
SBATCH --nodes=1
SBATCH --account=users
SBATCH --qos=mid_mdbf
SBATCH --partition=mid_mdbf
SBATCH --cpus-per-task=16
SBATCH -o /cta/users/merve.kaftancioglu/alternative_scenario/trimmomatic_as/out-err/%x-%j-slurm.out
SBATCH -e /cta/users/merve.kaftancioglu/alternative_scenario/trimmomatic_as/out-err/%x-%j-slurm.err

 Load FastQC module
module load trimmomatic-0.38-gcc-9.2.0-3evdnxc

 List of sample IDs (folders)
sample_ids=("h3k27me3_ms_cd4_ENCSR779JLY" "h3k27me3_ms_cd4_ENCSR993CTA" "h3k27me3_ms_cd8_ENCSR116FVG")

 Process each sample folder
for sample_id in "${sample_ids[@]}"; do
    input_folder="/cta/users/merve.kaftancioglu/alternative_scenario/seqs/$sample_id"
    output_folder="/cta/users/merve.kaftancioglu/alternative_scenario/trimmomatic_as/output/$sample_id"

     Create output folder if it doesn't exist
    mkdir -p "$output_folder"

     Process all .fastq.gz files in the specified folder
    for input_file in "$input_folder"/*.fastq.gz; do
        output_file="$output_folder/output_paired_$(basename "$input_file")"
        trimmomatic SE -threads 16 -phred33 -trimlog "$output_folder/trimmomatic_log.txt" "$input_file" "$output_file" HEADCROP:14
    done
done

""" !/bin/bash
SBATCH --job-name=trimmomatic_single
SBATCH --nodes=1
SBATCH --account=users
SBATCH --qos=mid_mdbf
SBATCH --partition=mid_mdbf
SBATCH --cpus-per-task=16
SBATCH -o /cta/users/merve.kaftancioglu/alternative_scenario/trimmomatic_as/out-err/%x-%j-slurm.out
SBATCH -e /cta/users/merve.kaftancioglu/alternative_scenario/trimmomatic_as/out-err/%x-%j-slurm.err

 Load FastQC module
module load trimmomatic-0.38-gcc-9.2.0-3evdnxc

 List of sample IDs (folders)
sample_ids=("h3k27me3_ms_cd8_ENCSR122JCM" "h3k27me3_ms_cd8_ENCSR216ZVA" "h3k27me3_ms_cd8_ENCSR284IKS")

 Process each sample folder
for sample_id in "${sample_ids[@]}"; do
    input_folder="/cta/users/merve.kaftancioglu/alternative_scenario/seqs/$sample_id"
    output_folder="/cta/users/merve.kaftancioglu/alternative_scenario/trimmomatic_as/output/$sample_id"

     Create output folder if it doesn't exist
    mkdir -p "$output_folder"

     Process all .fastq.gz files in the specified folder
    for input_file in "$input_folder"/*.fastq.gz; do
        output_file="$output_folder/output_paired_$(basename "$input_file")"
        trimmomatic SE -threads 16 -phred33 -trimlog "$output_folder/trimmomatic_log.txt" "$input_file" "$output_file" HEADCROP:14
    done
done

""" !/bin/bash
SBATCH --job-name=trimmomatic_single
SBATCH --nodes=1
SBATCH --account=users
SBATCH --qos=mid_mdbf
SBATCH --partition=mid_mdbf
SBATCH --cpus-per-task=16
SBATCH -o /cta/users/merve.kaftancioglu/alternative_scenario/trimmomatic_as/out-err/%x-%j-slurm.out
SBATCH -e /cta/users/merve.kaftancioglu/alternative_scenario/trimmomatic_as/out-err/%x-%j-slurm.err

 Load FastQC module
module load trimmomatic-0.38-gcc-9.2.0-3evdnxc

 List of sample IDs (folders)
sample_ids=("h3k27me3_ms_cd8_ENCSR385BOZ" "h3k27me3_ms_cd8_ENCSR521SFR" "h3k27me3_ms_nk_ENCSR469QVG")

 Process each sample folder
for sample_id in "${sample_ids[@]}"; do
    input_folder="/cta/users/merve.kaftancioglu/alternative_scenario/seqs/$sample_id"
    output_folder="/cta/users/merve.kaftancioglu/alternative_scenario/trimmomatic_as/output/$sample_id"

     Create output folder if it doesn't exist
    mkdir -p "$output_folder"

     Process all .fastq.gz files in the specified folder
    for input_file in "$input_folder"/*.fastq.gz; do
        output_file="$output_folder/output_paired_$(basename "$input_file")"
        trimmomatic SE -threads 16 -phred33 -trimlog "$output_folder/trimmomatic_log.txt" "$input_file" "$output_file" HEADCROP:14
    done
done

""" !/bin/bash
SBATCH --job-name=trimmomatic_single
SBATCH --nodes=1
SBATCH --account=users
SBATCH --qos=mid_mdbf
SBATCH --partition=mid_mdbf
SBATCH --cpus-per-task=16
SBATCH -o /cta/users/merve.kaftancioglu/alternative_scenario/trimmomatic_as/out-err/%x-%j-slurm.out
SBATCH -e /cta/users/merve.kaftancioglu/alternative_scenario/trimmomatic_as/out-err/%x-%j-slurm.err

 Load FastQC module
module load trimmomatic-0.38-gcc-9.2.0-3evdnxc

 List of sample IDs (folders)
sample_ids=("h3k27me3_ms_nk_ENCSR565WDW" "h3k27me3_normal_b_ENCSR589LHR" "h3k27me3_normal_cd4_ENCSR068YVZ")

 Process each sample folder
for sample_id in "${sample_ids[@]}"; do
    input_folder="/cta/users/merve.kaftancioglu/alternative_scenario/seqs/$sample_id"
    output_folder="/cta/users/merve.kaftancioglu/alternative_scenario/trimmomatic_as/output/$sample_id"

     Create output folder if it doesn't exist
    mkdir -p "$output_folder"

     Process all .fastq.gz files in the specified folder
    for input_file in "$input_folder"/*.fastq.gz; do
        output_file="$output_folder/output_paired_$(basename "$input_file")"
        trimmomatic SE -threads 16 -phred33 -trimlog "$output_folder/trimmomatic_log.txt" "$input_file" "$output_file" HEADCROP:14
    done
done

""" !/bin/bash
SBATCH --job-name=trimmomatic_single
SBATCH --nodes=1
SBATCH --account=users
SBATCH --qos=mid_mdbf
SBATCH --partition=mid_mdbf
SBATCH --cpus-per-task=16
SBATCH -o /cta/users/merve.kaftancioglu/alternative_scenario/trimmomatic_as/out-err/%x-%j-slurm.out
SBATCH -e /cta/users/merve.kaftancioglu/alternative_scenario/trimmomatic_as/out-err/%x-%j-slurm.err

 Load FastQC module
module load trimmomatic-0.38-gcc-9.2.0-3evdnxc

 List of sample IDs (folders)
sample_ids=("h3k27me3_normal_cd8_ENCSR720QRX" "h3k27me3_normal_nk_ENCSR639NIG" "h3k36me3_ms_b_ENCSR089VQL")

 Process each sample folder
for sample_id in "${sample_ids[@]}"; do
    input_folder="/cta/users/merve.kaftancioglu/alternative_scenario/seqs/$sample_id"
    output_folder="/cta/users/merve.kaftancioglu/alternative_scenario/trimmomatic_as/output/$sample_id"

     Create output folder if it doesn't exist
    mkdir -p "$output_folder"

     Process all .fastq.gz files in the specified folder
    for input_file in "$input_folder"/*.fastq.gz; do
        output_file="$output_folder/output_paired_$(basename "$input_file")"
        trimmomatic SE -threads 16 -phred33 -trimlog "$output_folder/trimmomatic_log.txt" "$input_file" "$output_file" HEADCROP:14
    done
done

""" !/bin/bash
SBATCH --job-name=trimmomatic_single
SBATCH --nodes=1
SBATCH --account=users
SBATCH --qos=mid_mdbf
SBATCH --partition=mid_mdbf
SBATCH --cpus-per-task=16
SBATCH -o /cta/users/merve.kaftancioglu/alternative_scenario/trimmomatic_as/out-err/%x-%j-slurm.out
SBATCH -e /cta/users/merve.kaftancioglu/alternative_scenario/trimmomatic_as/out-err/%x-%j-slurm.err

 Load FastQC module
module load trimmomatic-0.38-gcc-9.2.0-3evdnxc

 List of sample IDs (folders)
sample_ids=("h3k36me3_ms_b_ENCSR119PSR" "h3k36me3_ms_b_ENCSR238WFK" "h3k36me3_ms_b_ENCSR987OPY")

 Process each sample folder
for sample_id in "${sample_ids[@]}"; do
    input_folder="/cta/users/merve.kaftancioglu/alternative_scenario/seqs/$sample_id"
    output_folder="/cta/users/merve.kaftancioglu/alternative_scenario/trimmomatic_as/output/$sample_id"

     Create output folder if it doesn't exist
    mkdir -p "$output_folder"

     Process all .fastq.gz files in the specified folder
    for input_file in "$input_folder"/*.fastq.gz; do
        output_file="$output_folder/output_paired_$(basename "$input_file")"
        trimmomatic SE -threads 16 -phred33 -trimlog "$output_folder/trimmomatic_log.txt" "$input_file" "$output_file" HEADCROP:14
    done
done

""" !/bin/bash
SBATCH --job-name=trimmomatic_single
SBATCH --nodes=1
SBATCH --account=users
SBATCH --qos=mid_mdbf
SBATCH --partition=mid_mdbf
SBATCH --cpus-per-task=16
SBATCH -o /cta/users/merve.kaftancioglu/alternative_scenario/trimmomatic_as/out-err/%x-%j-slurm.out
SBATCH -e /cta/users/merve.kaftancioglu/alternative_scenario/trimmomatic_as/out-err/%x-%j-slurm.err

 Load FastQC module
module load trimmomatic-0.38-gcc-9.2.0-3evdnxc

 List of sample IDs (folders)
sample_ids=( "h3k36me3_ms_cd4_ENCSR276NGH" "h3k36me3_ms_cd4_ENCSR330CQU" "h3k36me3_ms_cd4_ENCSR471WZE")

 Process each sample folder
for sample_id in "${sample_ids[@]}"; do
    input_folder="/cta/users/merve.kaftancioglu/alternative_scenario/seqs/$sample_id"
    output_folder="/cta/users/merve.kaftancioglu/alternative_scenario/trimmomatic_as/output/$sample_id"

     Create output folder if it doesn't exist
    mkdir -p "$output_folder"

     Process all .fastq.gz files in the specified folder
    for input_file in "$input_folder"/*.fastq.gz; do
        output_file="$output_folder/output_paired_$(basename "$input_file")"
        trimmomatic SE -threads 16 -phred33 -trimlog "$output_folder/trimmomatic_log.txt" "$input_file" "$output_file" HEADCROP:14
    done
done

""" !/bin/bash
SBATCH --job-name=trimmomatic_single
SBATCH --nodes=1
SBATCH --account=users
SBATCH --qos=mid_mdbf
SBATCH --partition=mid_mdbf
SBATCH --cpus-per-task=16
SBATCH -o /cta/users/merve.kaftancioglu/alternative_scenario/trimmomatic_as/out-err/%x-%j-slurm.out
SBATCH -e /cta/users/merve.kaftancioglu/alternative_scenario/trimmomatic_as/out-err/%x-%j-slurm.err

 Load FastQC module
module load trimmomatic-0.38-gcc-9.2.0-3evdnxc

 List of sample IDs (folders)
sample_ids=("h3k36me3_ms_cd4_ENCSR473TGK" "h3k36me3_ms_cd4_ENCSR482VIB" "h3k36me3_ms_cd4_ENCSR532PXR")

 Process each sample folder
for sample_id in "${sample_ids[@]}"; do
    input_folder="/cta/users/merve.kaftancioglu/alternative_scenario/seqs/$sample_id"
    output_folder="/cta/users/merve.kaftancioglu/alternative_scenario/trimmomatic_as/output/$sample_id"

     Create output folder if it doesn't exist
    mkdir -p "$output_folder"

     Process all .fastq.gz files in the specified folder
    for input_file in "$input_folder"/*.fastq.gz; do
        output_file="$output_folder/output_paired_$(basename "$input_file")"
        trimmomatic SE -threads 16 -phred33 -trimlog "$output_folder/trimmomatic_log.txt" "$input_file" "$output_file" HEADCROP:14
    done
done

""" !/bin/bash
SBATCH --job-name=trimmomatic_single
SBATCH --nodes=1
SBATCH --account=users
SBATCH --qos=mid_mdbf
SBATCH --partition=mid_mdbf
SBATCH --cpus-per-task=16
SBATCH -o /cta/users/merve.kaftancioglu/alternative_scenario/trimmomatic_as/out-err/%x-%j-slurm.out
SBATCH -e /cta/users/merve.kaftancioglu/alternative_scenario/trimmomatic_as/out-err/%x-%j-slurm.err

 Load FastQC module
module load trimmomatic-0.38-gcc-9.2.0-3evdnxc

 List of sample IDs (folders)
sample_ids=("h3k36me3_ms_cd4_ENCSR785PKM" "h3k36me3_ms_cd4_ENCSR865FTW" "h3k36me3_ms_cd4_ENCSR898VJE")

 Process each sample folder
for sample_id in "${sample_ids[@]}"; do
    input_folder="/cta/users/merve.kaftancioglu/alternative_scenario/seqs/$sample_id"
    output_folder="/cta/users/merve.kaftancioglu/alternative_scenario/trimmomatic_as/output/$sample_id"

     Create output folder if it doesn't exist
    mkdir -p "$output_folder"

     Process all .fastq.gz files in the specified folder
    for input_file in "$input_folder"/*.fastq.gz; do
        output_file="$output_folder/output_paired_$(basename "$input_file")"
        trimmomatic SE -threads 16 -phred33 -trimlog "$output_folder/trimmomatic_log.txt" "$input_file" "$output_file" HEADCROP:14
    done
done

""" !/bin/bash
SBATCH --job-name=trimmomatic_single
SBATCH --nodes=1
SBATCH --account=users
SBATCH --qos=mid_mdbf
SBATCH --partition=mid_mdbf
SBATCH --cpus-per-task=16
SBATCH -o /cta/users/merve.kaftancioglu/alternative_scenario/trimmomatic_as/out-err/%x-%j-slurm.out
SBATCH -e /cta/users/merve.kaftancioglu/alternative_scenario/trimmomatic_as/out-err/%x-%j-slurm.err

 Load FastQC module
module load trimmomatic-0.38-gcc-9.2.0-3evdnxc

 List of sample IDs (folders)
sample_ids=("h3k36me3_ms_cd8_ENCSR239OWL" "h3k36me3_ms_cd8_ENCSR435MSB" "h3k36me3_ms_cd8_ENCSR631VIW")

 Process each sample folder
for sample_id in "${sample_ids[@]}"; do
    input_folder="/cta/users/merve.kaftancioglu/alternative_scenario/seqs/$sample_id"
    output_folder="/cta/users/merve.kaftancioglu/alternative_scenario/trimmomatic_as/output/$sample_id"

     Create output folder if it doesn't exist
    mkdir -p "$output_folder"

     Process all .fastq.gz files in the specified folder
    for input_file in "$input_folder"/*.fastq.gz; do
        output_file="$output_folder/output_paired_$(basename "$input_file")"
        trimmomatic SE -threads 16 -phred33 -trimlog "$output_folder/trimmomatic_log.txt" "$input_file" "$output_file" HEADCROP:14
    done
done

""" !/bin/bash
SBATCH --job-name=fastQC
SBATCH --nodes=1
SBATCH --account=users
SBATCH --qos=mid_mdbf
SBATCH --partition=mid_mdbf
SBATCH --cpus-per-task=16
SBATCH -o /cta/users/merve.kaftancioglu/alternative_scenario/post-trim_fastQC_as/out-err/%x-%j-slurm.out
SBATCH -e /cta/users/merve.kaftancioglu/alternative_scenario/post-trim_fastQC_as/out-err/%x-%j-slurm.err

 Load FastQC module
module load fastqc-0.11.7-gcc-9.2.0-nsjzjof

sample_ids=("h3k27ac_ms_b_ENCSR295QZX" "h3k27ac_ms_b_ENCSR364OIK" "h3k27ac_ms_b_ENCSR538URI" "h3k27ac_ms_b_ENCSR541PLO" "h3k27ac_ms_b_ENCSR617GMR")

for sample_id in "${sample_ids[@]}"; do
  echo "Predicting ${sample_id}..."
  input_directory="/cta/users/merve.kaftancioglu/alternative_scenario/trimmomatic_as/output/${sample_id}"
  output_dir="/cta/users/merve.kaftancioglu/alternative_scenario/post-trim_fastQC_as/output/${sample_id}"
  mkdir -p "${output_dir}"
   Fix: Expand the variable before using it
  output_dir_expanded="${output_directory}"
  fastqc --outdir "$output_dir_expanded" "$input_directory"/*.fastq.gz -threads 16
  echo "Finished ${sample_id}"
done


echo "DONE!""" "

""" !/bin/bash
SBATCH --job-name=fastQC
SBATCH --nodes=1
SBATCH --account=users
SBATCH --qos=mid_mdbf
SBATCH --partition=mid_mdbf
SBATCH --cpus-per-task=16
SBATCH -o /cta/users/merve.kaftancioglu/alternative_scenario/post-trim_fastQC_as/out-err/%x-%j-slurm.out
SBATCH -e /cta/users/merve.kaftancioglu/alternative_scenario/post-trim_fastQC_as/out-err/%x-%j-slurm.err

 Load FastQC module
module load fastqc-0.11.7-gcc-9.2.0-nsjzjof

sample_ids=("h3k27ac_ms_cd4_ENCSR200SSJ" "h3k27ac_ms_cd4_ENCSR322MTA" "h3k27ac_ms_cd4_ENCSR331WMS" "h3k27ac_ms_cd4_ENCSR350UKV" "h3k27ac_ms_cd4_ENCSR474PYR")

for sample_id in "${sample_ids[@]}"; do
  echo "Predicting ${sample_id}..."
  input_directory="/cta/users/merve.kaftancioglu/alternative_scenario/trimmomatic_as/output/${sample_id}"
  output_dir="/cta/users/merve.kaftancioglu/alternative_scenario/post-trim_fastQC_as/output/${sample_id}"
  mkdir -p "${output_dir}"
  fastqc --outdir "$output_dir" "$input_directory"/*.fastq.gz -threads 16
  echo "Finished ${sample_id}"
done

echo "DONE!""" "


""" !/bin/bash
SBATCH --job-name=fastQC
SBATCH --nodes=1
SBATCH --account=users
SBATCH --qos=mid_mdbf
SBATCH --partition=mid_mdbf
SBATCH --cpus-per-task=16
SBATCH -o /cta/users/merve.kaftancioglu/alternative_scenario/fastQC_as/out-err/%x-%j-slurm.out
SBATCH -e /cta/users/merve.kaftancioglu/alternative_scenario/fastQC_as/out-err/%x-%j-slurm.err


 Load FastQC module
module load fastqc-0.11.7-gcc-9.2.0-nsjzjof

substrate_ids=("h3k36me3_ms_cd8_ENCSR652WFO" "h3k36me3_ms_cd8_ENCSR757FGN" "h3k36me3_ms_nk_ENCSR158VSE" "h3k36me3_normal_cd8_ENCSR303SQG" "h3k4me1_ms_b_ENCSR084FXT" "h3k4me1_ms_cd4_ENCSR315MYP" "h3k4me1_ms_cd8_ENCSR815YZL" "h3k4me1_ms_cd8_ENCSR861FEC" "h3k4me1_ms_nk_ENCSR718LCI" "h3k4me1_ms_nk_ENCSR806HUT" "h3k4me1_normal_b_ENCSR156BXM" "h3k4me3_ms_b_ENCSR923EGO" "h3k4me3_ms_cd4_ENCSR180NCM" "h3k4me3_ms_cd4_ENCSR341QLC" "h3k4me3_ms_cd4_ENCSR603LTN" "h3k4me3_ms_cd4_ENCSR802MXQ" "h3k4me3_ms_cd8_ENCSR278QHR" "h3k4me3_ms_cd8_ENCSR516CKJ" "h3k4me3_ms_nk_ENCSR234BAD" "h3k4me3_ms_nk_ENCSR703TWY" "h3k4me3_normal_cd8_ENCSR123ZAT" "h3k4me3_normal_nk_ENCSR394JFQ" "h3k9me3_ms_cd4_ENCSR057IZD" "h3k9me3_ms_cd4_ENCSR550DPT" "h3k9me3_ms_nk_ENCSR611YIA" "h3k9me3_ms_nk_ENCSR667HUI" "h3k9me3_normal_b_ENCSR445LTM" "h3k9me3_normal_nk_ENCSR025UNZ")

for substrate_id in "${substrate_ids[@]}"; do
    echo "Predicting ${substrate_id}..."
    input_directory="/cta/users/merve.kaftancioglu/alternative_scenario/seqs/${substrate_id}"
    output_dir="/cta/users/merve.kaftancioglu/alternative_scenario/seqs/${substrate_id}"    
    mkdir -p "${output_dir}"
    fastqc --outdir "$output_directory" "$input_directory"/*.fastq.gz -threads 16
    echo "Finished ${substrate_id}"
done


echo "DONE!""" "



""" !/bin/bash
SBATCH --job-name=bwa_single
SBATCH --nodes=1
SBATCH --account=users
SBATCH --qos=cuda
SBATCH --partition=cuda
SBATCH --cpus-per-task=16
SBATCH -o /cta/users/merve.kaftancioglu/alternative_scenario/bwa_as/out-err/%x-%j-slurm.out
SBATCH -e /cta/users/merve.kaftancioglu/alternative_scenario/bwa_as/out-err/%x-%j-slurm.err

 Load BWA module
module load bwa-0.7.17-gcc-9.2.0-xwkclqy

 Define directories (consider adjusting paths)
input_dir="/cta/users/merve.kaftancioglu/alternative_scenario/trimmomatic_as/output"
output_dir="/cta/users/merve.kaftancioglu/alternative_scenario/bwa_as/output"
reference_genome="/cta/users/merve.kaftancioglu/alternative_scenario/bwa_as/reference/human_g1k_v37.fasta"

 Check if BWA index exists
if [ ! -f "$reference_genome.bwt" ]; then
  echo "Error: BWA index not found. Please create the index using: bwa mem index $reference_genome"
  exit 1
fi

 Define sample IDs (modify as needed)
sample_ids=("h3k36me_ms_cd4_ENCSR532PXR")

 Loop through all sample subdirectories
for sample_dir in "$input_dir"/*; do
   Check if directory name matches any sample ID
  if [[ " ${sample_ids[@]} " =~ " $(basename "$sample_dir") " ]]; then
    echo "Processing sample: $(basename "$sample_dir")"

     Create output directory for the sample (if it doesn't exist)
    mkdir -p "$output_dir/$(basename "$sample_dir")"

     Find FASTQ files within the sample directory
     Adjust the pattern below if your file names differ
    fastq_files=( "$sample_dir"/*.fastq.gz )

     Check if any FASTQ files found
    if [[ ${fastq_files[@]} -eq 0 ]]; then
      echo "Warning: No FASTQ files found in $sample_dir"
      continue  Skip to next iteration if no FASTQ files
    fi

     BWA mem command with paired-end or single-end handling
    if [[ ${fastq_files[@]} -eq 2 ]]; then
       Paired-end reads (assuming *_1.fastq.gz and *_2.fastq.gz naming)
      bwa mem $reference_genome ${fastq_files[0]} ${fastq_files[1]} > "$output_dir/$(basename "$sample_dir").sam"
    else
       Single-end read (assuming only one FASTQ file)
       Extract only the sample name from the directory path
      sample_name=$(basename "$sample_dir" /)
      bwa mem $reference_genome ${fastq_files[0]} > "$output_dir/$sample_name.sam"
    fi
  fi
done


""" !/bin/bash
SBATCH --job-name=sam_to_bam
SBATCH --nodes=1
SBATCH --account=users
SBATCH --qos=cuda
SBATCH --partition=cuda
SBATCH --cpus-per-task=16
SBATCH -o /cta/users/merve.kaftancioglu/alternative_scenario/samtools_as/out-err/%x-%j-slurm.out
SBATCH -e /cta/users/merve.kaftancioglu/alternative_scenario/samtools_as/out-err/%x-%j-slurm.err

 Load samtools module
module load samtools-1.9-gcc-9.2.0-w7pulwi

 Define directories (consider adjusting paths)
input_dir="/cta/users/merve.kaftancioglu/alternative_scenario/bwa_as/output"
output_dir="/cta/users/merve.kaftancioglu/alternative_scenario/samtools_as/output"

 Loop through all subdirectories in the input directory
for subdir in $input_dir/*; do
   Check if it's a directory (avoid hidden folders)
  if [ -d "$subdir" ]; then
     Get the sample name from the subdirectory name
    sample_name=$(basename "$subdir")

     Look for a .sam file within the subdirectory
    sam_file="$subdir/$sample_name.sam"

     Check if the SAM file exists
    if [ -f "$sam_file" ]; then
       Convert SAM to BAM using samtools view
      samtools view -bS "$sam_file" > "$output_dir/$sample_name.bam"

       Check for conversion errors (optional)
      if [ $? -ne 0 ]; then
        echo "Error converting $sam_file to BAM. See slurm error logs for details."
      fi
    fi
  fi
done

echo "DONE!""" "

""" !/bin/bash
SBATCH --job-name=star_alignment
SBATCH --nodes=1
SBATCH --account=users
SBATCH --qos=mid_mdbf
SBATCH --partition=mid_mdbf
SBATCH --cpus-per-task=16
SBATCH -o /cta/users/merve.kaftancioglu/alternative_scenario/bowtie2_as/out-err/%x-%j-slurm.out
SBATCH -e /cta/users/merve.kaftancioglu/alternative_scenario/bowtie2_as/out-err/%x-%j-slurm.err

 Load STAR module
module load star-2.7.0e-gcc-9.2.0-vynasg3

 Define directories
genome_dir="/cta/users/merve.kaftancioglu/alternative_scenario/star_as/index"   Update with your genome directory
input_dir="/cta/users/merve.kaftancioglu/alternative_scenario/trimmomatic_as/output"   Update with your input directory
output_dir="/cta/users/merve.kaftancioglu/alternative_scenario/star_as/output"   Update with your output directory

 Run STAR for each folder within the input directory
for sample_dir in "$input_dir"/*; do
   Check if it's a directory (avoid hidden files or non-directories)
  if [[ -d "$sample_dir" ]]; then
     Extract sample ID from directory name (modify if needed)
    sample_id=$(basename "$sample_dir")

     Define input and output paths with sample ID
    input_fastq="$sample_dir"/*.fastq.gz
    output_bam="$output_dir/$sample_id/$sample_id.bam"

     Run STAR for the current sample
    STAR --genomeDir "$genome_dir" \
         --runThreadN 16 \   Adjust threads based on your resource allocation
         --readFilesIn "$input_fastq" \   No colon after flag name
         --outFileNamePrefix "$output_dir/$sample_id/" \   Include trailing slash for proper path
         --outSAMtype BAM SortedByCoordinate \   No colon after flag name
         --outSAMunmapped Within \
         --outSAMattributes Standard \
         --alignEndsType EndToEnd
  fi
done

""" !/bin/bash
SBATCH --job-name=move
SBATCH --nodes=1
SBATCH --account=users
SBATCH --qos=cuda
SBATCH --partition=cuda
SBATCH --cpus-per-task=16
SBATCH -o /cta/users/merve.kaftancioglu/alternative_scenario/fastQC_as/out-err/%x-%j-slurm.out
SBATCH -e /cta/users/merve.kaftancioglu/alternative_scenario/fastQC_as/out-err/%x-%j-slurm.err

 Dry run (set to 'cp' for actual move)
DRY_RUN="cp"  

 Replace with the actual source directory path
 Use curly brackets only if $sample_id contains special characters or spaces
SOURCE_DIR="/cta/users/merve.kaftancioglu/alternative_scenario/trimmomatic_as/output/{$sample_id}"

 Destination directory 
DEST_DIR="/cta/users/merve.kaftancioglu/alternative_scenario/post-trim_fastQC_as/output/{$sample_id}"

 Create destination directory if it doesn't exist
mkdir -p "$DEST_DIR"

 Find zip files (ensure there are some in the source directory)
zip_files=$(find "$SOURCE_DIR" -type f -name "*.zip")

 Check if any zip files were found
if [ -z "$zip_files" ]; then
  echo "No zip files found in source directory: $SOURCE_DIR"
  exit 0
fi

 Loop through zip files and copy them
for file in $zip_files
do
  $DRY_RUN cp "$file" "$DEST_DIR"
done

 ~/.bashrc: executed by bash(1) for non-login shells.
 see /usr/share/doc/bash/examples/startup-files (in the package bash-doc)
 for examples

 If not running interactively, don't do anything
case $- in
    *i*) ;;
      *) return;;
esac

 don't put duplicate lines or lines starting with space in the history.
 See bash(1) for more options
HISTCONTROL=ignoreboth

 append to the history file, don't overwrite it
shopt -s histappend

 for setting history length see HISTSIZE and HISTFILESIZE in bash(1)
HISTSIZE=1000
HISTFILESIZE=2000

 check the window size after each command and, if necessary,
 update the values of LINES and COLUMNS.
shopt -s checkwinsize

 If set, the pattern "**" used in a pathname expansion context will
 match all files and zero or more directories and subdirectories.
shopt -s globstar

 make less more friendly for non-text input files, see lesspipe(1)
[ -x /usr/bin/lesspipe ] && eval "$(SHELL=/bin/sh lesspipe)"

 set variable identifying the chroot you work in (used in the prompt below)
if [ -z "${debian_chroot:-}" ] && [ -r /etc/debian_chroot ]; then
    debian_chroot=$(cat /etc/debian_chroot)
fi

 set a fancy prompt (non-color, unless we know we "want" color)
case "$TERM" in
    xterm-color|*-256color) color_prompt=yes;;
esac

 uncomment for a colored prompt, if the terminal has the capability; turned
 off by default to not distract the user: the focus in a terminal window
 should be on the output of commands, not on the prompt
force_color_prompt=yes

if [ -n "$force_color_prompt" ]; then
    if [ -x /usr/bin/tput ] && tput setaf 1 >&/dev/null; then
	 We have color support; assume it's compliant with Ecma-48
	 (ISO/IEC-6429). (Lack of such support is extremely rare, and such
	 a case would tend to support setf rather than setaf.)
	color_prompt=yes
    else
	color_prompt=
    fi
fi

if [ "$color_prompt" = yes ]; then
    PS1='${debian_chroot:+($debian_chroot)}\[\033[01;32m\]\u@\h\[\033[00m\]:\[\033[01;34m\]\w\[\033[00m\]\$ '
else
    PS1='${debian_chroot:+($debian_chroot)}\u@\h:\w\$ '
fi
unset color_prompt force_color_prompt

 If this is an xterm set the title to user@host:dir
case "$TERM" in
xterm*|rxvt*)
    PS1="\[\e]0;${debian_chroot:+($debian_chroot)}\u@\h: \w\a\]$PS1"
    ;;
*)
    ;;
esac

 enable color support of ls and also add handy aliases
if [ -x /usr/bin/dircolors ]; then
    test -r ~/.dircolors && eval "$(dircolors -b ~/.dircolors)" || eval "$(dircolors -b)"
    alias ls='ls --color=auto'
    alias dir='dir --color=auto'
    alias vdir='vdir --color=auto'

    alias grep='grep --color=auto'
    alias fgrep='fgrep --color=auto'
    alias egrep='egrep --color=auto'
fi

 colored GCC warnings and errors
export GCC_COLORS='error=01;31:warning=01;35:note=01;36:caret=01;32:locus=01:quote=01'

 some more ls aliases
alias ll='ls -alF'
alias la='ls -A'
alias l='ls -CF'

 Add an "alert" alias for long running commands.  Use like so:
   sleep 10; alert
alias alert='notify-send --urgency=low -i "$([ $? = 0 ] && echo terminal || echo error)" "$(history|tail -n1|sed -e '\''s/^\s*[0-9]\+\s*//;s/[;&|]\s*alert$//'\'')"'

 Alias definitions.
 You may want to put all your additions into a separate file like
 ~/.bash_aliases, instead of adding them here directly.
 See /usr/share/doc/bash-doc/examples in the bash-doc package.

if [ -f ~/.bash_aliases ]; then
    . ~/.bash_aliases
fi

 enable programmable completion features (you don't need to enable
 this, if it's already enabled in /etc/bash.bashrc and /etc/profile
 sources /etc/bash.bashrc).
if ! shopt -oq posix; then
  if [ -f /usr/share/bash-completion/bash_completion ]; then
    . /usr/share/bash-completion/bash_completion
  elif [ -f /etc/bash_completion ]; then
    . /etc/bash_completion
  fi
fi

 >>> conda initialize >>>
 !! Contents within this block are managed by 'conda init' !!
__conda_setup="$('/home/merve/anaconda3/bin/conda' 'shell.bash' 'hook' 2> /dev/null)"
if [ $? -eq 0 ]; then
    eval "$__conda_setup"
else
    if [ -f "/home/merve/anaconda3/etc/profile.d/conda.sh" ]; then
        . "/home/merve/anaconda3/etc/profile.d/conda.sh"
    else
        export PATH="/home/merve/anaconda3/bin:$PATH"
    fi
fi
unset __conda_setup
 <<< conda initialize <<<


""" !/bin/bash
SBATCH --job-name=star_single
SBATCH --nodes=1
SBATCH --account=users
SBATCH --qos=mid_mdbf
SBATCH --partition=mid_mdbf
SBATCH --cpus-per-task=16
SBATCH -o /cta/users/merve.kaftancioglu/alternative_scenario/star_as/out-err/%x-%j-slurm.out
SBATCH -e /cta/users/merve.kaftancioglu/alternative_scenario/star_as/out-err/%x-%j-slurm.err

 Define paths 
genome_fasta="/cta/users/merve.kaftancioglu/alternative_scenario/star_as/reference/human_g1k_v37.fasta"
star_index_dir="/cta/users/merve.kaftancioglu/alternative_scenario/star_as/index"

 Check if genome fasta file exists
if [ ! -f "$genome_fasta" ]; then
  echo "Error: Genome fasta file not found at: $genome_fasta"
  exit 1
fi

 Check if star_index_dir exists (and create if not)
if [ ! -d "$star_index_dir" ]; then
  echo "Creating star_genome_index directory: $star_index_dir"
  mkdir -p "$star_index_dir"
fi

 Define number of threads (adjust based on your system)
threads=8

 STAR command
star_cmd="STAR \
  --runThreadN $threads \
  --runMode genomeGenerate \
  --genomeDir $star_index_dir \
  --genomeFastaFiles $genome_fasta"

 Run STAR and handle potential errors
echo "Building STAR genome index..."
$star_cmd

if [ $? -ne 0 ]; then
  echo "Error: STAR failed to generate genome index."
  exit 1
fi

echo "STAR genome index created successfully in: $star_index_dir"

""" !/bin/bash
SBATCH --job-name=bowtie2_alignment
SBATCH --nodes=1
SBATCH --account=users
SBATCH --qos=mid_mdbf
SBATCH --partition=mid_mdbf
SBATCH --cpus-per-task=16
SBATCH -o /cta/users/merve.kaftancioglu/alternative_scenario/bowtie2_as/out-err/%x-%j-slurm.out
SBATCH -e /cta/users/merve.kaftancioglu/alternative_scenario/bowtie2_as/out-err/%x-%j-slurm.err

 Load Bowtie2 module
module load bowtie2-2.3.5.1-gcc-9.2.0-ae4dozk

 Sample IDs 
sample_ids=("h3k27ac_ms_b_ENCSR295QZX" "h3k27ac_ms_b_ENCSR364OIK" "h3k27ac_ms_b_ENCSR538URI")

 Directories
input_dir="/cta/users/merve.kaftancioglu/alternative_scenario/trimmomatic_as/output"
output_dir="/cta/users/merve.kaftancioglu/alternative_scenario/bowtie2_as/output"
index_dir="/cta/users/merve.kaftancioglu/alternative_scenario/bowtie2_as/index"  

 Create Bowtie2 index 
bowtie2-build "$reference_genome" "$index_dir"

 Loop through samples
for sample_id in "${sample_ids[@]}"; do
  echo "Processing ${sample_id}..."

   Create output directory 
  mkdir -p "${output_dir}/${sample_id}"

   Combine input and output paths with sample ID
  input_fastq="${input_dir}/${sample_id}"/*.fastq.gz
  output_bam="${output_dir}/${sample_id}/${sample_id}.bam"

   Run Bowtie2 with trimming
  bowtie2 --threads 16 -p 16 --trim3 30 -x "$index_dir" -U "$input_fastq" | samtools view -bS - > "$output_bam"
  
   Convert BAM to SAM
  input_bam="$output_bam"
  output_sam="${output_dir}/${sample_id}/${sample_id}.sam"
  samtools view -h -o "$output_sam" "$input_b

  echo "Finished processing ${sample_id}"
done

echo "All samples processed!"

""" !/bin/bash
SBATCH --job-name=star_single
SBATCH --nodes=1
SBATCH --account=users
SBATCH --qos=mid_mdbf
SBATCH --partition=mid_mdbf
SBATCH --cpus-per-task=16
SBATCH -o /cta/users/merve.kaftancioglu/alternative_scenario/bwa_as/out-err/%x-%j-slurm.out
SBATCH -e /cta/users/merve.kaftancioglu/alternative_scenario/bwa_as/out-err/%x-%j-slurm.err

 Load BWA module
module load bwa-0.7.17-gcc-9.2.0-xwkclqy

 Define directories
input_dir="/cta/users/merve.kaftancioglu/alternative_scenario/trimmomatic_as/output"   Folder containing sample subdirectories
output_dir="/cta/users/merve.kaftancioglu/alternative_scenario/bwa_as/output"   Output directory

 Path to reference genome FASTA file
reference_genome="/cta/users/merve.kaftancioglu/alternative_scenario/bwa_as/reference/human_g1k_v37.fasta"

 Create BWA index 
bwa index $reference_genome

 Loop through all sample subdirectories
for sample_dir in "$input_dir"/*; do
   Extract sample ID from directory name (assuming it matches)
  sample_id=$(basename "$sample_dir")

   Check if it's a directory (skip non-directories)
  if [[ -d "$sample_dir" ]]; then
    echo "Processing sample: $sample_id"

     Find FASTQ files within the sample directory
     Adjust the pattern below if your file names differ
    fastq_files=( "$sample_dir"/*.fastq.gz )

     Check if any FASTQ files found
    if [[ ${fastq_files[@]} -eq 0 ]]; then
      echo "Warning: No FASTQ files found in $sample_dir"
      continue   Skip to next iteration if no FASTQ files
    fi

     BWA command with paired-end or single-end handling
    if [[ ${fastq_files[@]} -eq 2 ]]; then
       Paired-end reads (assuming *_1.fastq.gz and *_2.fastq.gz naming)
      bwa mem $reference_genome ${fastq_files[0]} ${fastq_files[1]} > "$output_dir/$sample_id.sam"
    else
       Single-end read (assuming only one FASTQ file)
      bwa mem $reference_genome ${fastq_files[0]} > "$output_dir/$sample_id.sam"
    fi
  fi
done

 BWA Bwasw 
bwa bwasw index_prefix input_reads.fasta -t $SLURM_NTASKS_PER_NODE > bwa_bwasw_alignments.sam

""" !/bin/bash
SBATCH --job-name=bwa_single
SBATCH --nodes=1
SBATCH --account=users
SBATCH --qos=cuda
SBATCH --partition=cuda
SBATCH --cpus-per-task=16
SBATCH -o /cta/users/merve.kaftancioglu/alternative_scenario/bwa_as/out-err/%x-%j-slurm.out
SBATCH -e /cta/users/merve.kaftancioglu/alternative_scenario/bwa_as/out-err/%x-%j-slurm.err

 Load BWA module
module load bwa-0.7.17-gcc-9.2.0-xwkclqy

 Define directories
input_dir="/cta/users/merve.kaftancioglu/alternative_scenario/trimmomatic_as/output"   Folder containing sample subdirectories
output_dir="/cta/users/merve.kaftancioglu/alternative_scenario/bwa_as/output"   Output directory

 Path to reference genome FASTA file
reference_genome="/cta/users/merve.kaftancioglu/alternative_scenario/bwa_as/reference/human_g1k_v37.fasta"

 Assuming you already have the BWA index files created (remove the commented line below if you need to create them)
 bwa mem index $reference_genome   Only uncomment this line if you need to recreate the index

 Loop through all sample subdirectories
for sample_dir in "$input_dir"/*; do
   Extract sample ID from directory name (assuming it matches)
  sample_id=$(basename "$sample_dir")

   Check if it's a directory (skip non-directories)
  if [[ -d "$sample_dir" ]]; then
    echo "Processing sample: $sample_id"

     Create output directory for the sample (if it doesn't exist)
    mkdir -p "$output_dir/$sample_id"   Create directory with sample ID in the output directory

     Find FASTQ files within the sample directory
     Adjust the pattern below if your file names differ
    fastq_files=( "$sample_dir"/*.fastq.gz )

     Check if any FASTQ files found
    if [[ ${fastq_files[@]} -eq 0 ]]; then
      echo "Warning: No FASTQ files found in $sample_dir"
      continue   Skip to next iteration if no FASTQ files
    fi

     BWA mem command with paired-end or single-end handling
    if [[ ${fastq_files[@]} -eq 2 ]]; then
       Paired-end reads (assuming *_1.fastq.gz and *_2.fastq.gz naming)
      bwa mem $reference_genome ${fastq_files[0]} ${fastq_files[1]} > "$output_dir/$sample_id.sam"
    else
       Single-end read (assuming only one FASTQ file)
      bwa mem $reference_genome ${fastq_files[0]} > "$output_dir/$sample_dir.sam"   Use sample directory name for single-end output
    fi
  fi
done


""" !/bin/bash
SBATCH --job-name=bwa_single
SBATCH --nodes=1
SBATCH --account=users
SBATCH --qos=cuda
SBATCH --partition=cuda
SBATCH --cpus-per-task=16
SBATCH -o /cta/users/merve.kaftancioglu/alternative_scenario/bwa_as/out-err/%x-%j-slurm.out
SBATCH -e /cta/users/merve.kaftancioglu/alternative_scenario/bwa_as/out-err/%x-%j-slurm.err

 Load BWA module
module load bwa-0.7.17-gcc-9.2.0-xwkclqy

 Define directories (consider adjusting paths)
input_dir="/cta/users/merve.kaftancioglu/alternative_scenario/trimmomatic_as/output"
output_dir="/cta/users/merve.kaftancioglu/alternative_scenario/bwa_as/output"
reference_genome="/cta/users/merve.kaftancioglu/alternative_scenario/bwa_as/reference/human_g1k_v37.fasta"

 Check if BWA index exists
if [ ! -f "$reference_genome.bwt" ]; then
  echo "Error: BWA index not found. Please create the index using: bwa mem index $reference_genome"
  exit 1
fi

 Define sample IDs (modify as needed)
sample_ids=("h3k27ac_ms_b_ENCSR295QZX" "h3k27ac_ms_b_ENCSR538URI" "h3k27ac_ms_b_ENCSR617GMR" "h3k4me1_ms_b_ENCSR084FXT" "h3k4me1_ms_b_ENCSR480ITK" "h3k4me1_ms_b_ENCSR759DHN" "h3k4me1_ms_b_ENCSR900IAM" "h3k4me1_ms_b_ENCSR931JYC" "h3k4me3_ms_b_ENCSR260CRI" "h3k4me3_ms_b_ENCSR461QMZ" "h3k4me3_ms_b_ENCSR530AVK" "h3k4me3_ms_b_ENCSR713ZYF" "h3k4me3_ms_b_ENCSR923EGO" "h3k9me3_ms_b_ENCSR445DNH" "h3k9me3_ms_b_ENCSR486SZN" "h3k9me3_ms_b_ENCSR983MKX" "h3k27ac_ms_b_ENCSR364OIK" "h3k27ac_ms_b_ENCSR541PLO" "h3k36me3_ms_cd4_ENCSR471WZE" "h3k36me3_ms_cd4_ENCSR473TGK" "h3k4me1_ms_nk_ENCSR482UGV" "h3k4me1_ms_nk_ENCSR718LCI" "h3k4me1_ms_nk_ENCSR806HUT" "h3k4me3_ms_nk_ENCSR234BAD" "h3k4me3_ms_nk_ENCSR703TWY" "h3k4me3_ms_nk_ENCSR75")

 Loop through all sample subdirectories
for sample_dir in "$input_dir"/*; do
   Check if directory name matches any sample ID
  if [[ " ${sample_ids[@]} " =~ " $(basename "$sample_dir") " ]]; then
    echo "Processing sample: $(basename "$sample_dir")"

     Create output directory for the sample (if it doesn't exist)
    mkdir -p "<span class="math-inline">output\_dir/</span>(basename "$sample_dir")"

     Find FASTQ files within the sample directory
     Adjust the pattern below if your file names differ
    fastq_files=( "$sample_dir"/*.fastq.gz )

     Check if any FASTQ files found
    if [[ ${fastq_files[@]} -eq 0 ]]; then
      echo "Warning: No FASTQ files found in $sample_dir"
      continue  Skip to next iteration if no FASTQ files
    fi

     BWA mem command with paired-end or single-end handling
    if [[ ${fastq_files[@]} -eq 2 ]]; then
       Paired-end reads (assuming *_1.fastq.gz and *_2.fastq.gz naming)
      bwa mem $reference_genome ${fastq_files[0]} ${fastq_files[1]} > "$output_dir/$sample_id.sam"
    else
       Single-end read (assuming only one FASTQ file)
       Extract only the sample name from the directory path
      sample_name=$(basename "$sample_dir" /)
      bwa mem $reference_genome ${fastq_files[0]} > "$output_dir/$sample_name.sam"
    fi

     Convert SAM to BAM format
    samtools view -S -b "$output_dir/$sample_id.sam" > "$output_dir/$sample_id.bam"
    
     Remove the intermediate SAM file (optional, saves space)
    rm "$output_dir/$sample_id.sam"   Uncomment this line if you want to remove SAM files
  fi
done

echo "BWA alignment completed! Both SAM and BAM files generated."

""" !/bin/bash
SBATCH --job-name=bwa_single
SBATCH --nodes=1
SBATCH --account=users
SBATCH --qos=cuda
SBATCH --partition=cuda
SBATCH --cpus-per-task=16
SBATCH -o /cta/users/merve.kaftancioglu/alternative_scenario/bwa_as/out-err/%x-%j-slurm.out
SBATCH -e /cta/users/merve.kaftancioglu/alternative_scenario/bwa_as/out-err/%x-%j-slurm.err

 Load BWA module
module load bwa-0.7.17-gcc-9.2.0-xwkclqy

 Define directories (consider adjusting paths)
input_dir="/cta/users/merve.kaftancioglu/alternative_scenario/trimmomatic_as/output"
output_dir="/cta/users/merve.kaftancioglu/alternative_scenario/bwa_as/output"
reference_genome="/cta/users/merve.kaftancioglu/alternative_scenario/bwa_as/reference/human_g1k_v37.fasta"

 Check if BWA index exists
if [ ! -f "$reference_genome.bwt" ]; then
  echo "Error: BWA index not found. Please create the index using: bwa mem index $reference_genome"
  exit 1
fi

 Define sample IDs (modify as needed)
sample_ids=("h3k27ac_ms_b_ENCSR295QZX" "h3k27ac_ms_b_ENCSR538URI" "h3k27ac_ms_b_ENCSR617GMR" "h3k4me1_ms_b_ENCSR084FXT" "h3k4me1_ms_b_ENCSR480ITK" "h3k4me1_ms_b_ENCSR759DHN" "h3k4me1_ms_b_ENCSR900IAM" "h3k4me1_ms_b_ENCSR931JYC" "h3k4me3_ms_b_ENCSR260CRI" "h3k4me3_ms_b_ENCSR461QMZ" "h3k4me3_ms_b_ENCSR530AVK" "h3k4me3_ms_b_ENCSR713ZYF" "h3k4me3_ms_b_ENCSR923EGO" "h3k9me3_ms_b_ENCSR445DNH" "h3k9me3_ms_b_ENCSR486SZN" "h3k9me3_ms_b_ENCSR983MKX" "h3k27ac_ms_b_ENCSR364OIK" "h3k27ac_ms_b_ENCSR541PLO" "h3k36me3_ms_cd4_ENCSR471WZE" "h3k36me3_ms_cd4_ENCSR473TGK" "h3k4me1_ms_nk_ENCSR482UGV" "h3k4me1_ms_nk_ENCSR718LCI" "h3k4me1_ms_nk_ENCSR806HUT" "h3k4me3_ms_nk_ENCSR234BAD" "h3k4me3_ms_nk_ENCSR703TWY" "h3k4me3_ms_nk_ENCSR75")

 Loop through all sample subdirectories
for sample_dir in "$input_dir"/*; do
   Check if directory name matches any sample ID
  if [[ " ${sample_ids[@]} " =~ " $(basename "$sample_dir") " ]]; then
    echo "Processing sample: $(basename "$sample_dir")"

     Create output directory for the sample (if it doesn't exist)
    mkdir -p "$output_dir/$(basename "$sample_dir")"

     Find FASTQ files within the sample directory
     Adjust the pattern below if your file names differ
    fastq_files=( "$sample_dir"/*.fastq.gz )

     Check if any FASTQ files found
    if [[ ${fastq_files[@]} -eq 0 ]]; then
      echo "Warning: No FASTQ files found in $sample_dir"
      continue  Skip to next iteration if no FASTQ files
    fi

     BWA mem command with paired-end or single-end handling
    if [[ ${fastq_files[@]} -eq 2 ]]; then
       Paired-end reads (assuming *_1.fastq.gz and *_2.fastq.gz naming)
      sample_id=$(basename "$sample_dir")
      bwa mem "$reference_genome" "${fastq_files[0]}" "${fastq_files[1]}" > "$output_dir/$sample_id.sam"
    else
       Single-end read (assuming only one FASTQ file)
       Extract only the sample name from the directory path
      sample_name=$(basename "$sample_dir")
      bwa mem "$reference_genome" "${fastq_files[0]}" > "$output_dir/$sample_name.sam"
    fi

     Convert SAM to BAM format
    samtools view -S -b "$output_dir/$sample_id.sam" > "$output_dir/$sample_id.bam"
    
     Remove the intermediate SAM file (optional, saves space)
    rm "$output_dir/$sample_id.sam"   Uncomment this line if you want to remove SAM files
  fi
done

echo "BWA alignment completed! Both SAM and BAM files generated."

""" !/bin/bash
SBATCH --job-name=fastQC
SBATCH --nodes=1
SBATCH --account=users
SBATCH --qos=short_mdbf
SBATCH --partition=short_mdbf
SBATCH --cpus-per-task=8
SBATCH -o /cta/users/merve.kaftancioglu/fastQC_M2/job_outputs/%x-%j-slurm.out
SBATCH -e /cta/users/merve.kaftancioglu/fastQC_M2/job_outputs/%x-%j-slurm.err


input_directory="/cta/users/cta/users/merve.kaftancioglu/alternative_scenario/seqs/h3k27me3_ms"
output_directory="/cta/users/merve.kaftancioglu/fastQC_M2/job_outputs"

 Load FastQC module
module load fastqc-0.11.7-gcc-9.2.0-nsjzjof

 Run FastQC on .fastq files in input directory and write output to output directory
fastqc --outdir "$output_directory" "$input_directory"/*.fastq -threads 8

""" !/bin/bash
SBATCH --job-name=chipseq_pipeline
SBATCH --nodes=1
SBATCH --account=users
SBATCH --qos=cuda
SBATCH --partition=cuda
SBATCH --cpus-per-task=16
SBATCH -o /cta/users/merve.kaftancioglu/alternative_scenario/processed_seqs/out-err/%x-%j-slurm.out
SBATCH -e /cta/users/merve.kaftancioglu/alternative_scenario/processed_seqs/out-err/%x-%j-slurm.err

cd /cta/users/merve.kaftancioglu/alternative_scenario/nf-core/chipseq-master

nextflow run main.nf -profile singularity --outdir /output

""" !/bin/bash
SBATCH --job-name=trimmomatic
SBATCH --nodes=1
SBATCH --account=users
SBATCH --qos=short_mdbf
SBATCH --partition=short_mdbf
SBATCH --cpus-per-task=16
SBATCH -o /cta/users/merve.kaftancioglu/alternative_scenario/fastQC_as/out-err/%x-%j-slurm.out
SBATCH -e /cta/users/merve.kaftancioglu/alternative_scenario/fastQC_as/out-err/%x-%j-slurm.err

 Load FastQC module
module load fastqc-0.11.7-gcc-9.2.0-nsjzjof

substrate_ids=("h3k27ac_ms_b_ENCSR295QZX" "h3k27ac_ms_b_ENCSR364OIK" "h3k27ac_ms_b_ENCSR538URI" "h3k27ac_ms_b_ENCSR541PLO" "h3k27ac_ms_b_ENCSR617GMR" "h3k27ac_ms_cd4_ENCSR200SSJ" "h3k27ac_ms_cd4_ENCSR322MTA" "h3k27ac_ms_cd4_ENCSR331WMS" "h3k27ac_ms_cd4_ENCSR350UKV" "h3k27ac_ms_cd4_ENCSR474PYR" "h3k27ac_ms_cd4_ENCSR520QDR" "h3k27ac_ms_cd4_ENCSR540XNK" "h3k27ac_ms_cd4_ENCSR705VSO" "h3k27ac_ms_cd4_ENCSR832UMM" "h3k27ac_ms_cd8_ENCSR078ATS" "h3k27ac_ms_cd8_ENCSR348YRH" "h3k27ac_ms_cd8_ENCSR458TOW" "h3k27ac_ms_cd8_ENCSR476IPR" "h3k27ac_ms_cd8_ENCSR787HDF" "h3k27ac_ms_cd8_ENCSR923JIB" "h3k27ac_ms_nk_ENCSR246TTM" "h3k27ac_ms_nk_ENCSR469CFL" "h3k27ac_ms_nk_ENCSR746AIX" "h3k27ac_normal_b_ENCSR685KZA" "h3k27ac_normal_cd4_ENCSR819NCZ" "h3k27ac_normal_cd8_ENCSR976RWL" "h3k27ac_normal_nk_ENCSR977FMZ" "h3k27me3_ms_b_ENCSR009ZRH" "h3k27me3_ms_b_ENCSR182NLA" "h3k27me3_ms_b_ENCSR272YVX" "h3k27me3_ms_b_ENCSR649FUX" "h3k27me3_ms_b_ENCSR842WWX" "h3k27me3_ms_cd4_ENCSR277XYX" "h3k27me3_ms_cd4_ENCSR526TNC" "h3k27me3_ms_cd4_ENCSR592EKF" "h3k27me3_ms_cd4_ENCSR613UFD" "h3k27me3_ms_cd4_ENCSR740SDR" "h3k27me3_ms_cd4_ENCSR779JLY" "h3k27me3_ms_cd4_ENCSR993CTA" "h3k27me3_ms_cd8_ENCSR116FVG" "h3k27me3_ms_cd8_ENCSR122JCM" "h3k27me3_ms_cd8_ENCSR216ZVA" "h3k27me3_ms_cd8_ENCSR284IKS" "h3k27me3_ms_cd8_ENCSR385BOZ" "h3k27me3_ms_cd8_ENCSR521SFR" "h3k27me3_ms_nk_ENCSR469QVG" "h3k27me3_ms_nk_ENCSR565WDW" "h3k27me3_normal_b_ENCSR589LHR" "h3k27me3_normal_cd4_ENCSR068YVZ" "h3k27me3_normal_cd8_ENCSR720QRX" "h3k27me3_normal_nk_ENCSR639NIG" "h3k36me3_ms_b_ENCSR089VQL" "h3k36me3_ms_b_ENCSR119PSR" "h3k36me3_ms_b_ENCSR238WFK" "h3k36me3_ms_b_ENCSR987OPY" "h3k36me3_ms_cd4_ENCSR276NGH" "h3k36me3_ms_cd4_ENCSR330CQU" "h3k36me3_ms_cd4_ENCSR471WZE" "h3k36me3_ms_cd4_ENCSR473TGK" "h3k36me3_ms_cd4_ENCSR482VIB" "h3k36me3_ms_cd4_ENCSR532PXR" "h3k36me3_ms_cd4_ENCSR785PKM" "h3k36me3_ms_cd4_ENCSR865FTW" "h3k36me3_ms_cd4_ENCSR898VJE" "h3k36me3_ms_cd8_ENCSR239OWL" "h3k36me3_ms_cd8_ENCSR435MSB" "h3k36me3_ms_cd8_ENCSR631VIW" "h3k36me3_ms_cd8_ENCSR652WFO" "h3k36me3_ms_cd8_ENCSR757FGN" "h3k36me3_ms_nk_ENCSR158VSE" "h3k36me3_ms_nk_ENCSR245KON" "h3k36me3_ms_nk_ENCSR530YDY" "h3k36me3_normal_b_ENCSR831AXK" "h3k36me3_normal_cd4_ENCSR640OKB" "h3k36me3_normal_cd8_ENCSR303SQG" "h3k36me3_normal_nk_ENCSR056GJY" "h3k4me1_ms_b_ENCSR084FXT" "h3k4me1_ms_b_ENCSR480ITK" "h3k4me1_ms_b_ENCSR759DHN" "h3k4me1_ms_b_ENCSR900IAM" "h3k4me1_ms_b_ENCSR931JYC" "h3k4me1_ms_cd4_ENCSR036YVP" "h3k4me1_ms_cd4_ENCSR043JIA" "h3k4me1_ms_cd4_ENCSR093VUP" "h3k4me1_ms_cd4_ENCSR124IGC" "h3k4me1_ms_cd4_ENCSR315MYP" "h3k4me1_ms_cd4_ENCSR458QEY" "h3k4me1_ms_cd4_ENCSR485FBT" "h3k4me1_ms_cd4_ENCSR641ZFV" "h3k4me1_ms_cd4_ENCSR687CQX" "h3k4me1_ms_cd8_ENCSR231XAP" "h3k4me1_ms_cd8_ENCSR572XTB" "h3k4me1_ms_cd8_ENCSR788TEF" "h3k4me1_ms_cd8_ENCSR815YZL" "h3k4me1_ms_cd8_ENCSR861FEC" "h3k4me1_ms_cd8_ENCSR940PHE" "h3k4me1_ms_nk_ENCSR482UGV" "h3k4me1_ms_nk_ENCSR718LCI" "h3k4me1_ms_nk_ENCSR806HUT" "h3k4me1_normal_b_ENCSR156BXM" "h3k4me1_normal_cd4_ENCSR102SOR" "h3k4me1_normal_cd8_ENCSR217SHH" "h3k4me1_normal_nk_ENCSR277YKG" "h3k4me3_ms_b_ENCSR260CRI" "h3k4me3_ms_b_ENCSR461QMZ" "h3k4me3_ms_b_ENCSR530AVK" "h3k4me3_ms_b_ENCSR713ZYF" "h3k4me3_ms_b_ENCSR923EGO" "h3k4me3_ms_cd4_ENCSR180NCM" "h3k4me3_ms_cd4_ENCSR341QLC" "h3k4me3_ms_cd4_ENCSR482TGI" "h3k4me3_ms_cd4_ENCSR486XJK" "h3k4me3_ms_cd4_ENCSR496LKR" "h3k4me3_ms_cd4_ENCSR603LTN" "h3k4me3_ms_cd4_ENCSR802MXQ" "h3k4me3_ms_cd4_ENCSR878YHM" "h3k4me3_ms_cd4_ENCSR954ZLD" "h3k4me3_ms_cd8_ENCSR231ZZH" "h3k4me3_ms_cd8_ENCSR278QHR" "h3k4me3_ms_cd8_ENCSR516CKJ" "h3k4me3_ms_cd8_ENCSR535YYH" "h3k4me3_ms_cd8_ENCSR741XAE" "h3k4me3_ms_cd8_ENCSR848XJL" "h3k4me3_ms_nk_ENCSR234BAD" "h3k4me3_ms_nk_ENCSR703TWY" "h3k4me3_ms_nk_ENCSR755ZGS" "h3k4me3_normal_b_ENCSR791CAF" "h3k4me3_normal_cd4_ENCSR537KJA" "h3k4me3_normal_cd8_ENCSR123ZAT" "h3k4me3_normal_nk_ENCSR394JFQ" "h3k9me3_ms_b_ENCSR445DNH" "h3k9me3_ms_b_ENCSR486SZN" "h3k9me3_ms_b_ENCSR983MKX" "h3k9me3_ms_cd4_ENCSR057IZD" "h3k9me3_ms_cd4_ENCSR550DPT" "h3k9me3_ms_cd4_ENCSR677OEF" "h3k9me3_ms_cd4_ENCSR729ZXH" "h3k9me3_ms_cd4_ENCSR851RJV" "h3k9me3_ms_cd4_ENCSR919DFZ" "h3k9me3_ms_cd4_ENCSR953STZ" "h3k9me3_ms_cd4_ENCSR959VZU" "h3k9me3_ms_cd8_ENCSR101USF" "h3k9me3_ms_cd8_ENCSR354GNT" "h3k9me3_ms_cd8_ENCSR377QYB" "h3k9me3_ms_cd8_ENCSR733LCG" "h3k9me3_ms_cd8_ENCSR980KTW" "h3k9me3_ms_nk_ENCSR061ATV" "h3k9me3_ms_nk_ENCSR611YIA" "h3k9me3_ms_nk_ENCSR667HUI" "h3k9me3_normal_b_ENCSR445LTM" "h3k9me3_normal_cd4_ENCSR433EWI" "h3k9me3_normal_cd8_ENCSR294HTM" "h3k9me3_normal_nk_ENCSR025UNZ")

for substrate_id in "${substrate_ids[@]}"; do
    echo "Predicting ${substrate_id}..."
    input_directory="/cta/users/merve.kaftancioglu/alternative_scenario/seqs/${substrate_id}"
    output_dir="/cta/users/merve.kaftancioglu/alternative_scenario/fastQC_as/output/${substrate_id}"    
    mkdir -p "${output_dir}"
    trimmomatic PE -threads 16 -phred33 -trimlog "${output_dir}/trimmomatic_log.txt" "$input_directory"/R1.fastq.gz "$input_directory"/R2.fastq.gz "${output_dir}/output_forward_paired.fastq.gz" "${output_dir}/output_forward_unpaired.fastq.gz" "${output_dir}/output_reverse_paired.fastq.gz" "${output_dir}/output_reverse_unpaired.fastq.gz" ILLUMINACLIP:/path/to/adapters.fa:2:30:10 LEADING:3 TRAILING:3 SLIDINGWINDOW:4:15 MINLEN:36
    echo "Finished ${substrate_id}"
done


""" !/bin/bash
SBATCH --job-name=multiQC
SBATCH --nodes=1
SBATCH --account=users
SBATCH --qos=cuda
SBATCH --partition=cuda
SBATCH --cpus-per-task=16
SBATCH -o /cta/users/merve.kaftancioglu/alternative_scenario/post-trim_fastQC_as/output/out-err/%x-%j-slurm.out
SBATCH -e /cta/users/merve.kaftancioglu/alternative_scenario/post-trim_fastQC_as/output/out-err/%x-%j-slurm.err

 Specify the main directory containing the folders to merge
main_directory="/cta/users/merve.kaftancioglu/alternative_scenario/post-trim_fastQC_as/output_merged"

 Specify the directory where you want to merge the files
merge_directory="/cta/users/merve.kaftancioglu/alternative_scenario/post-trim_fastQC_as/output_merged"

 Copy all files from subdirectories to the merge directory
find "$main_directory" -type f -exec cp {} "$merge_directory" \;

echo "All files merged into: $merge_directory"

""" !/bin/bash
SBATCH --job-name=star_single
SBATCH --nodes=1
SBATCH --account=users
SBATCH --qos=mid_mdbf
SBATCH --partition=mid_mdbf
SBATCH --cpus-per-task=16
SBATCH -o /cta/users/merve.kaftancioglu/alternative_scenario/star_as/out-err/%x-%j-slurm.out
SBATCH -e /cta/users/merve.kaftancioglu/alternative_scenario/star_as/out-err/%x-%j-slurm.err

 Load STAR module
module load star-2.7.0e-gcc-9.2.0-vynasg3

 Define paths (replace with your actual paths)
genome_fasta="/cta/users/merve.kaftancioglu/alternative_scenario/star_as/reference/human_g1k_v37.fasta"
star_index_dir="/cta/users/merve.kaftancioglu/star_as/index"

 Check if genome fasta file exists
if [ ! -f "$genome_fasta" ]; then
  echo "Error: Genome fasta file not found at: $genome_fasta"
  exit 1
fi

 Check if star_index_dir exists (and create if not)
if [ ! -d "$star_index_dir" ]; then
  echo "Creating star_genome_index directory: $star_index_dir"
  mkdir -p "$star_index_dir"
fi

 Define number of threads (adjust based on your system)
threads=8

 STAR command
star_cmd="STAR \
  --runThreadN $threads \
  --runMode genomeGenerate \
  --genomeDir $star_index_dir \
  --genomeFastaFiles $genome_fasta"

 Run STAR and handle potential errors
echo "Building STAR genome index..."
$star_cmd

if [ $? -ne 0 ]; then
  echo "Error: STAR failed to generate genome index."
  exit 1
fi

echo "STAR genome index created successfully in: $star_index_dir"

""" !/bin/bash
SBATCH --job-name=nf-core
SBATCH --nodes=1
SBATCH --account=users
SBATCH --qos=cuda
SBATCH --partition=cuda
SBATCH --cpus-per-task=16
SBATCH -o /cta/users/merve.kaftancioglu/alternative_scenario/processed_seqs/out-err/%x-%j-slurm.out
SBATCH -e /cta/users/merve.kaftancioglu/alternative_scenario/processed_seqs/out-err/%x-%j-slurm.err

 Load necessary modules (if applicable for your cluster)
module load docker

 Specify the Docker image (replace with nf-core/chipseq version you downloaded)
docker_image=nf-core/chipseq:latest

 Navigate to the directory containing your Slurm script and pipeline files
cd /cta/users/merve.kaftancioglu/alternative_scenario/processed_seqs/all_data

 Create symbolic links for all fastq files in the current directory
for file in *.fastq; do ln -s "$file" all_fastq_files; done

 Run the nextflow command within the container
docker run --rm -v /cta/users/merve.kaftancioglu/alternative_scenario/trimmomatic_as/output

 Create a directory to hold symbolic links (soft links) to your fastq files (optional)
 This step is optional, you can comment it out if you don't want to create a separate directory
 mkdir all_fastq_files

 Loop through all subdirectories under processed_seqs
for sample_dir in */ ; do
   Navigate to the current subdirectory
  cd "$sample_dir"

   Loop through all fastq.gz files in the current subdirectory
  for file in *.fastq.gz; do
     Create a symbolic link (optional, uncomment if using the directory)
     ln -s "$file" ../all_fastq_files   Uncomment if using all_fastq_files directory
    
     Construct the full path to the fastq file (alternative to symbolic links)
    full_fastq_path=$(pwd)/"$file"
    
     You can echo the full path for verification (optional)
     echo "Found fastq: $full_fastq_path"
  done
  
   Move back to the processed_seqs directory after processing each subdirectory
  cd ..
done

""" !/bin/bash
SBATCH --job-name=picard
SBATCH --nodes=1
SBATCH --account=users
SBATCH --qos=cuda
SBATCH --partition=cuda
SBATCH --cpus-per-task=16
SBATCH -o /cta/users/merve.kaftancioglu/alternative_scenario/picard_as/out-err/%x-%j-slurm.out
SBATCH -e /cta/users/merve.kaftancioglu/alternative_scenario/picard_as/out-err/%x-%j-slurm.err

 Load BWA module

 Load required modules
module load picard samtools

 Define input and output file variables
MAIN_DIR="/cta/users/merve.kaftancioglu/alternative_scenario/samtools_as/output"   Replace with path to main directory
REALIGNED_BAMS_DIR="$MAIN_DIR/realigned_bams"   Output directory for realigned BAMs
TEMP_DIR="$MAIN_DIR/temp"   Optional temporary directory (create if needed)

 **Step 1: Mark duplicates on each individual BAM**

 Loop through subdirectories
for subdir in $MAIN_DIR/*; do
   Check if subdirectory is a directory (avoid hidden files)
  if [[ -d "$subdir" ]]; then
     Extract subdirectory name
    subdir_name=$(basename "$subdir")

     Find BAM files with wildcard
    bam_file="$subdir/$subdir_name*.bam"

     Check if BAM file exists
    if [[ -f "$bam_file" ]]; then
       Mark duplicates and place output in subdirectory
      mark_dups_cmd="picard MarkDuplicates \
        I=$bam_file \
        O=$subdir/$subdir_name.marked_dups.bam \
        M=$subdir/$subdir_name.marked_dups_metrics.txt \
        CREATE_INDEX=true"

      echo "Marking duplicates for $subdir_name..."
      $mark_dups_cmd
    else
      echo "No BAM file found in $subdir."
    fi
  fi
done

 **Step 2: Merge marked duplicate BAMs (optional)**

 List of temporary marked duplicate BAMs
marked_dups_bams="$TEMP_DIR/$SAMPLE_NAME*.marked_dups.bam"

 Merge marked duplicate BAMs
merge_bam_cmd="picard MergeSamFiles \
  O=$REALIGNED_BAMS_DIR/$SAMPLE_NAME.merged.marked_dups.bam \
  $marked_dups_bams"

echo "Merging marked duplicate BAMs..."
$merge_bam_cmd

 **Step 3: Re-mark duplicates on the merged BAM**

 Re-mark duplicates on the merged BAM
re_mark_dups_cmd="picard MarkDuplicates \
  I=$REALIGNED_BAMS_DIR/$SAMPLE_NAME.merged.marked_dups.bam \
  O=$REALIGNED_BAMS_DIR/$SAMPLE_NAME.realigned.bam \
  M=$REALIGNED_BAMS_DIR/$SAMPLE_NAME.realigned_metrics.txt \
  CREATE_INDEX=true"

echo "Re-marking duplicates on merged BAM..."
$re_mark_dups_cmd

 **Optional Step: Clean up temporary files (if used)**
 rm -rf $TEMP_DIR   Uncomment to remove temporary files

echo "DONE!""" "

 **Step 3: Re-mark duplicates on merged BAMs (optional)**

 If you performed merging, adjust commands to work with merged BAMs.

 **Optional Cleanup (if using temporary directory)**
 rm -rf $TEMP_DIR   Uncomment to remove temporary files

echo "DONE!""" "


""" !/bin/bash
SBATCH --job-name=bowtie2_alignment
SBATCH --nodes=1
SBATCH --account=users
SBATCH --qos=mid_mdbf
SBATCH --partition=mid_mdbf
SBATCH --cpus-per-task=16
SBATCH -o /cta/users/merve.kaftancioglu/alternative_scenario/bowtie2_as/out-err/%x-%j-slurm.out
SBATCH -e /cta/users/merve.kaftancioglu/alternative_scenario/bowtie2_as/out-err/%x-%j-slurm.err

 Load Bowtie2 module
module load bowtie2-2.3.5.1-gcc-9.2.0-ae4dozk

 Define directories
input_dir="/cta/users/merve.kaftancioglu/alternative_scenario/trimmomatic_as/output"
output_dir="/cta/users/merve.kaftancioglu/alternative_scenario/bowtie2_as/output"

 Path to reference genome FASTA file
reference_genome="/cta/users/merve.kaftancioglu/alternative_scenario/bowtie2_as/reference/human_g1k_v37.fasta"

 Create Bowtie2 index (if not already done)
bowtie2-build "$reference_genome" /cta/users/merve.kaftancioglu/alternative_scenario/bowtie2_as/index

 Loop through all FASTQ files in the input directory using find
for input_fastq in $(find "${input_dir}" -type f -name "*.fastq.gz" ! -name "*.*"); do

  sample_id="${input_fastq*/}"  
  sample_id="${sample_id//.fastq.gz}" 

  echo "Processing ${sample_id}..."

   Create output directory (if it doesn't exist)
  mkdir -p "${output_dir}/${sample_id}"

   Output SAM file path
  output_sam="${output_dir}/${sample_id}/${sample_id}.sam"

   Run Bowtie2 with the current FASTQ file, specifying output format (-U) as SAM
  bowtie2 --threads 16 -p 16 --trim3 30 -x /cta/users/merve.kaftancioglu/alternative_scenario/bowtie2_as/index -U $input_fastq -U $output_sam
   Note the change: -U $output_sam specifies the output SAM file
done

""" !/bin/bash
SBATCH --job-name=bowtie2_alignment
SBATCH --nodes=1
SBATCH --account=users
SBATCH --qos=mid_mdbf
SBATCH --partition=mid_mdbf
SBATCH --cpus-per-task=16
SBATCH -o /cta/users/merve.kaftancioglu/alternative_scenario/bowtie2_as/out-err/%x-%j-slurm.out
SBATCH -e /cta/users/merve.kaftancioglu/alternative_scenario/bowtie2_as/out-err/%x-%j-slurm.err

 Load Bowtie2 module
module load bowtie2-2.3.5.1-gcc-9.2.0-ae4dozk

 Define directories
input_dir="/cta/users/merve.kaftancioglu/alternative_scenario/trimmomatic_as/output"
output_dir="/cta/users/merve.kaftancioglu/alternative_scenario/bowtie2_as/output"

 Path to reference genome FASTA file
reference_genome="/cta/users/merve.kaftancioglu/alternative_scenario/bowtie2_as/reference/human_g1k_v37.fasta"

 Create Bowtie2 index (if not already done)
bowtie2-build "$reference_genome" /cta/users/merve.kaftancioglu/alternative_scenario/bowtie2_as/index

 Loop through all FASTQ files in the input directory
for input_fastq in "${input_dir}"/*.fastq.gz; do
   Extract sample ID from filename (assuming filenames start with sample ID)
   Modify this part to match your specific filename format if needed
  sample_id="${input_fastq*/}"   Double  removes everything before the last slash (/)
  sample_id="${sample_id//.fastq.gz}"   Remove .fastq.gz extension

  echo "Processing ${sample_id}..."

   Create output directory (if it doesn't exist)
  mkdir -p "${output_dir}/${sample_id}"

   Output BAM file path
  output_bam="${output_dir}/${sample_id}/${sample_id}.bam"

   Run Bowtie2 with the current FASTQ file
  bowtie2 -x /cta/users/merve.kaftancioglu/alternative_scenario/bowtie2_as/index -U "$input_fastq" -S "$output_bam"
done

""" !/bin/bash
SBATCH --job-name=bowtie2_alignment
SBATCH --nodes=1
SBATCH --account=users
SBATCH --qos=mid_mdbf
SBATCH --partition=mid_mdbf
SBATCH --cpus-per-task=16
SBATCH -o /cta/users/merve.kaftancioglu/alternative_scenario/bowtie2_as/out-err/%x-%j-slurm.out
SBATCH -e /cta/users/merve.kaftancioglu/alternative_scenario/bowtie2_as/out-err/%x-%j-slurm.err

 Load Bowtie2 module
module load bowtie2-2.3.5.1-gcc-9.2.0-ae4dozk

 Define directories
input_dir="/cta/users/merve.kaftancioglu/alternative_scenario/trimmomatic_as/output"
output_dir="/cta/users/merve.kaftancioglu/alternative_scenario/bowtie2_as/output"

 Path to reference genome FASTA file
reference_genome="/cta/users/merve.kaftancioglu/alternative_scenario/bowtie2_as/reference/human_g1k_v37.fasta"

 Create Bowtie2 index (if not already done)
bowtie2-build "$reference_genome" /cta/users/merve.kaftancioglu/alternative_scenario/bowtie2_as/index

 Loop through all FASTQ files in the input directory using find
for input_fastq in $(find "${input_dir}" -type f -name "*.fastq.gz" ! -name "*.*"); do

  sample_id="${input_fastq*/}"  
  sample_id="${sample_id//.fastq.gz}" 

  echo "Processing ${sample_id}..."

   Create output directory (if it doesn't exist)
  mkdir -p "${output_dir}/${sample_id}"

   Output BAM file path
  output_bam="${output_dir}/${sample_id}/${sample_id}.bam"

   Run Bowtie2 with the current FASTQ file
  bowtie2 --threads 16 -p 16 --trim3 30 -x reference_genome -U $input_fastq -S $output_bam
   bowtie2 -x /cta/users/merve.kaftancioglu/alternative_scenario/bowtie2_as/index -U "$input_fastq" -S "$output_bam"
done

""" !/bin/bash
SBATCH --job-name=bowtie2_alignment
SBATCH --nodes=1
SBATCH --account=users
SBATCH --qos=mid_mdbf
SBATCH --partition=mid_mdbf
SBATCH --cpus-per-task=16
SBATCH -o /cta/users/merve.kaftancioglu/alternative_scenario/bowtie2_as/out-err/%x-%j-slurm.out
SBATCH -e /cta/users/merve.kaftancioglu/alternative_scenario/bowtie2_as/out-err/%x-%j-slurm.err

 Load Bowtie2 module
module load bowtie2-2.3.5.1-gcc-9.2.0-ae4dozk

 Define directories
input_dir="/cta/users/merve.kaftancioglu/alternative_scenario/trimmomatic_as/output"
output_dir="/cta/users/merve.kaftancioglu/alternative_scenario/bowtie2_as/output"

 Path to reference genome FASTA file
reference_genome="/cta/users/merve.kaftancioglu/alternative_scenario/bowtie2_as/reference/human_g1k_v37.fasta"

 Create Bowtie2 index (if not already done)
bowtie2-build "$reference_genome" /cta/users/merve.kaftancioglu/alternative_scenario/bowtie2_as/index

 Loop through all FASTQ files in the input directory using find
for input_fastq in $(find "${input_dir}" -type f -name "*.fastq.gz" ! -name "*.*"); do

  sample_id="${input_fastq*/}"  
  sample_id="${sample_id//.fastq.gz}" 

  echo "Processing ${sample_id}..."

   Create output directory (if it doesn't exist)
  mkdir -p "${output_dir}/${sample_id}"

   Output SAM file path
  output_sam="${output_dir}/${sample_id}/${sample_id}.sam"

   Run Bowtie2 with the current FASTQ file, specifying output format (-S) as SAM
  bowtie2 --threads 16 -p 16 --trim3 30 -x /cta/users/merve.kaftancioglu/alternative_scenario/bowtie2_as/index -U $input_fastq -S $output_sam
done

""" !/bin/bash
SBATCH --job-name=bowtie2_longread_alignment
SBATCH --nodes=1
SBATCH --account=users
SBATCH --qos=mid_mdbf
SBATCH --partition=mid_mdbf
SBATCH --cpus-per-task=16
SBATCH -o /cta/users/merve.kaftancioglu/alternative_scenario/bowtie2_as/out-err/%x-%j-slurm.out
SBATCH -e /cta/users/merve.kaftancioglu/alternative_scenario/bowtie2_as/out-err/%x-%j-slurm.err

 Load Bowtie2 module
module load bowtie2-2.3.5.1-gcc-9.2.0-ae4dozk

 Define directories (consider adjusting paths)
input_dir="/cta/users/merve.kaftancioglu/alternative_scenario/trimmomatic_as/output"
output_dir="/cta/users/merve.kaftancioglu/alternative_scenario/bowtie2_as/longread_output"   Separate output directory for long reads

 Path to reference genome FASTA file
reference_genome="/cta/users/merve.kaftancioglu/alternative_scenario/bowtie2_as/reference/human_g1k_v37.fasta"

 Create Bowtie2 index (if not already done)
bowtie2-build "$reference_genome" /cta/users/merve.kaftancioglu/alternative_scenario/bowtie2_as/index

 Loop through all FASTQ files in the input directory
for input_fastq in $(find "${input_dir}" -type f -name "*.fastq.gz" ! -name "*.*"); do

  sample_id="${input_fastq*/}"  
  sample_id="${sample_id//.fastq.gz}" 

  echo "Processing ${sample_id}..."

   Create output directory (if it doesn't exist)
  mkdir -p "${output_dir}/${sample_id}"

   Output SAM file path
  output_sam="${output_dir}/${sample_id}/${sample_id}.sam"

   Run Bowtie2 with the current FASTQ file, specifying output format (-S) as SAM and "--very-sensitive" flag
  bowtie2 --threads 16 -p 16 --trim3 30 --very-sensitive -x /cta/users/merve.kaftancioglu/alternative_scenario/bowtie2_as/index -U $input_fastq -S $output_sam

   Optional: Check for successful completion (add this block if desired)
  if [ $? -eq 0 ]; then
    echo "Alignment for ${sample_id} completed successfully."
  else
    echo "Error occurred during alignment for ${sample_id}. Check slurm error logs for details."
  fi
done

@ECHO OFF

REM Get user inputs
SET /P pipeline_name=Enter the nf-core pipeline name (should be chipseq): 
IF NOT "%pipeline_name%"=="chipseq" (
  ECHO Error: Please enter "chipseq" as the pipeline name.
  GOTO EXIT
)

SET /P data_format=Enter data format (single-end or paired-end): 
SET /P input_path=/cta/users/merve.kaftancioglu/alternative_scenario/seqs: 
SET /P genome_path= /cta/users/merve.kaftancioglu/alternative_scenario/nf-core/reference/human_g1k_v37.fasta:
 SET /P blacklist_path=Enter the path to your blacklist BED file (optional): 
SET /P output_path=/cta/users/merve.kaftancioglu/alternative_scenario/nf-core/chipseq-2.0.0/output: 
SET /P profile=singularity:3.8.4  

REM Build the command with user inputs
SET command=nextflow run nf-core/chipseq -r 2.0.0 

REM Add options based on data format (only for paired-end)
IF "%data_format%"=="paired-end" (
  SET command=%command% --paired
)

SET command=%command% --input %input_path% --genome %genome_path% --outdir %output_path% -profile %profile%

REM Add optional blacklist path (uncomment if needed)
IF NOT "%blacklist_path%"=="" (
  SET command=%command% --blacklist %blacklist_path%
)

REM Echo the final command for verification
ECHO The command to be executed:
ECHO %command%

REM Pause before execution (optional)
PAUSE

REM Run the pipeline command
%command%

ECHO Pipeline execution finished!

PAUSE

:EXIT



""" !/bin/bash
SBATCH --job-name=bwa_single
SBATCH --nodes=1
SBATCH --account=users
SBATCH --qos=cuda
SBATCH --partition=cuda
SBATCH --cpus-per-task=16
SBATCH -o /cta/users/merve.kaftancioglu/alternative_scenario/processed_seqs/out-err/%x-%j-slurm.out
SBATCH -e /cta/users/merve.kaftancioglu/alternative_scenario/processed_seqs/out-err/%x-%j-slurm.err

 Load BWA module
module load rsync-3.1.3-gcc-9.2.0-mrwcim2

 Define directories (consider adjusting paths)
input_dir="/cta/users/merve.kaftancioglu/alternative_scenario/processed_seqs/all_data"
output_dir="/cta/users/merve.kaftancioglu/alternative_scenario/bwa_as/output"
reference_genome="/cta/users/merve.kaftancioglu/alternative_scenario/bwa_as/reference/human_g1k_v37.fasta"

 Check if BWA index exists
if [ ! -f "$reference_genome.bwt" ]; then
  echo "Error: BWA index not found. Please create the index using: bwa mem index $reference_genome"
  exit 1
fi

 Define sample IDs (modify as needed)
sample_ids=("h3k27ac_ms_b_ENCSR295QZX" "h3k27ac_ms_b_ENCSR538URI" "h3k27ac_ms_b_ENCSR617GMR" "h3k4me1_ms_b_ENCSR084FXT" "h3k4me1_ms_b_ENCSR480ITK" "h3k4me1_ms_b_ENCSR759DHN" "h3k4me1_ms_b_ENCSR900IAM" "h3k4me1_ms_b_ENCSR931JYC" "h3k4me3_ms_b_ENCSR260CRI" "h3k4me3_ms_b_ENCSR461QMZ" "h3k4me3_ms_b_ENCSR530AVK" "h3k4me3_ms_b_ENCSR713ZYF" "h3k4me3_ms_b_ENCSR923EGO" "h3k9me3_ms_b_ENCSR445DNH" "h3k9me3_ms_b_ENCSR486SZN" "h3k9me3_ms_b_ENCSR983MKX" "h3k27ac_ms_b_ENCSR364OIK" "h3k27ac_ms_b_ENCSR541PLO" "h3k36me3_ms_cd4_ENCSR471WZE" "h3k36me3_ms_cd4_ENCSR473TGK" "h3k4me1_ms_nk_ENCSR482UGV" "h3k4me1_ms_nk_ENCSR718LCI" "h3k4me1_ms_nk_ENCSR806HUT" "h3k4me3_ms_nk_ENCSR234BAD" "h3k4me3_ms_nk_ENCSR703TWY" "h3k4me3_ms_nk_ENCSR755ZGS" "h3k9me3_ms_nk_ENCSR061ATV" "h3k9me3_ms_nk_ENCSR611YIA" "h3k9me3_ms_nk_ENCSR667HUI" "h3k27ac_ms_nk_ENCSR246TTM")

 Loop through all sample subdirectories
for sample_dir in "$input_dir"/*; do
   Check if directory name matches any sample ID
  if [[ " ${sample_ids[@]} " =~ " $(basename "$sample_dir") " ]]; then
    echo "Processing sample: $(basename "$sample_dir")"

     Create output directory for the sample (if it doesn't exist)
    mkdir -p "$output_dir/$(basename "$sample_dir")"

     Find FASTQ files within the sample directory
     Adjust the pattern below if your file names differ
    fastq_files=( "$sample_dir"/*.fastq.gz )

     Check if any FASTQ files found
    if [[ ${fastq_files[@]} -eq 0 ]]; then
      echo "Warning: No FASTQ files found in $sample_dir"
      continue  Skip to next iteration if no FASTQ files
    fi

     BWA mem command with paired-end or single-end handling
    if [[ ${fastq_files[@]} -eq 2 ]]; then
       Paired-end reads (assuming *_1.fastq.gz and *_2.fastq.gz naming)
      bwa mem $reference_genome ${fastq_files[0]} ${fastq_files[1]} > "$output_dir/$(basename "$sample_dir").sam"
    else
       Single-end read (assuming only one FASTQ file)
       Extract only the sample name from the directory path
      sample_name=$(basename "$sample_dir" /)
      bwa mem $reference_genome ${fastq_files[0]} > "$output_dir/$sample_name.sam"
    fi
  fi
done

""" !/bin/bash
SBATCH --job-name=trimmomatic_single
SBATCH --nodes=1
SBATCH --account=users
SBATCH --qos=mid_mdbf
SBATCH --partition=mid_mdbf
SBATCH --cpus-per-task=16
SBATCH -o /cta/users/merve.kaftancioglu/alternative_scenario/trimmomatic_as/out-err/%x-%j-slurm.out
SBATCH -e /cta/users/merve.kaftancioglu/alternative_scenario/trimmomatic_as/out-err/%x-%j-slurm.err

 Load FastQC module
module load trimmomatic-0.38-gcc-9.2.0-3evdnxc

 List of sample IDs (folders)
sample_ids=("h3k27ac_ms_b_ENCSR541PLO" "h3k27ac_ms_b_ENCSR617GMR" "h3k27ac_ms_cd4_ENCSR200SSJ" "h3k27ac_ms_cd4_ENCSR322MTA" "h3k27ac_ms_cd4_ENCSR331WMS" "h3k27ac_ms_cd4_ENCSR350UKV" "h3k27ac_ms_cd4_ENCSR474PYR" "h3k27ac_ms_cd4_ENCSR520QDR" "h3k27ac_ms_cd4_ENCSR540XNK" "h3k27ac_ms_cd4_ENCSR705VSO" "h3k27ac_ms_cd4_ENCSR832UMM" "h3k27ac_ms_cd8_ENCSR078ATS" "h3k27ac_ms_cd8_ENCSR348YRH" "h3k27ac_ms_cd8_ENCSR458TOW" "h3k27ac_ms_cd8_ENCSR476IPR" "h3k27ac_ms_cd8_ENCSR787HDF" "h3k27ac_ms_cd8_ENCSR923JIB" "h3k27ac_ms_nk_ENCSR246TTM" "h3k27ac_ms_nk_ENCSR469CFL" "h3k27ac_ms_nk_ENCSR746AIX" "h3k27ac_normal_b_ENCSR685KZA" "h3k27ac_normal_cd4_ENCSR819NCZ" "h3k27ac_normal_cd8_ENCSR976RWL" "h3k27ac_normal_nk_ENCSR977FMZ" "h3k27me3_ms_b_ENCSR009ZRH" "h3k27me3_ms_b_ENCSR182NLA" "h3k27me3_ms_b_ENCSR272YVX" "h3k27me3_ms_b_ENCSR649FUX" "h3k27me3_ms_b_ENCSR842WWX" "h3k27me3_ms_cd4_ENCSR277XYX" "h3k27me3_ms_cd4_ENCSR526TNC" "h3k27me3_ms_cd4_ENCSR592EKF" "h3k27me3_ms_cd4_ENCSR613UFD" "h3k27me3_ms_cd4_ENCSR740SDR" "h3k27me3_ms_cd4_ENCSR779JLY" "h3k27me3_ms_cd4_ENCSR993CTA" "h3k27me3_ms_cd8_ENCSR116FVG" "h3k27me3_ms_cd8_ENCSR122JCM" "h3k27me3_ms_cd8_ENCSR216ZVA" "h3k27me3_ms_cd8_ENCSR284IKS" "h3k27me3_ms_cd8_ENCSR385BOZ" "h3k27me3_ms_cd8_ENCSR521SFR" "h3k27me3_ms_nk_ENCSR469QVG" "h3k27me3_ms_nk_ENCSR565WDW" "h3k27me3_normal_b_ENCSR589LHR" "h3k27me3_normal_cd4_ENCSR068YVZ" "h3k27me3_normal_cd8_ENCSR720QRX" "h3k27me3_normal_nk_ENCSR639NIG" "h3k36me3_ms_b_ENCSR089VQL" "h3k36me3_ms_b_ENCSR119PSR" "h3k36me3_ms_b_ENCSR238WFK" "h3k36me3_ms_b_ENCSR987OPY" "h3k36me3_ms_cd4_ENCSR276NGH" "h3k36me3_ms_cd4_ENCSR330CQU" "h3k36me3_ms_cd4_ENCSR471WZE" "h3k36me3_ms_cd4_ENCSR473TGK" "h3k36me3_ms_cd4_ENCSR482VIB" "h3k36me3_ms_cd4_ENCSR532PXR" "h3k36me3_ms_cd4_ENCSR785PKM" "h3k36me3_ms_cd4_ENCSR865FTW" "h3k36me3_ms_cd4_ENCSR898VJE" "h3k36me3_ms_cd8_ENCSR239OWL" "h3k36me3_ms_cd8_ENCSR435MSB" "h3k36me3_ms_cd8_ENCSR631VIW" "h3k36me3_ms_cd8_ENCSR652WFO" "h3k36me3_ms_cd8_ENCSR757FGN" "h3k36me3_ms_nk_ENCSR158VSE" "h3k36me3_ms_nk_ENCSR245KON" "h3k36me3_ms_nk_ENCSR530YDY" "h3k36me3_normal_b_ENCSR831AXK" "h3k36me3_normal_cd4_ENCSR640OKB" "h3k36me3_normal_cd8_ENCSR303SQG" "h3k36me3_normal_nk_ENCSR056GJY" "h3k4me1_ms_b_ENCSR084FXT" "h3k4me1_ms_b_ENCSR480ITK" "h3k4me1_ms_b_ENCSR759DHN" "h3k4me1_ms_b_ENCSR900IAM" "h3k4me1_ms_b_ENCSR931JYC" "h3k4me1_ms_cd4_ENCSR036YVP" "h3k4me1_ms_cd4_ENCSR043JIA" "h3k4me1_ms_cd4_ENCSR093VUP" "h3k4me1_ms_cd4_ENCSR124IGC" "h3k4me1_ms_cd4_ENCSR315MYP" "h3k4me1_ms_cd4_ENCSR458QEY" "h3k4me1_ms_cd4_ENCSR485FBT" "h3k4me1_ms_cd4_ENCSR641ZFV" "h3k4me1_ms_cd4_ENCSR687CQX" "h3k4me1_ms_cd8_ENCSR231XAP" "h3k4me1_ms_cd8_ENCSR572XTB" "h3k4me1_ms_cd8_ENCSR788TEF" "h3k4me1_ms_cd8_ENCSR815YZL" "h3k4me1_ms_cd8_ENCSR861FEC" "h3k4me1_ms_cd8_ENCSR940PHE" "h3k4me1_ms_nk_ENCSR482UGV" "h3k4me1_ms_nk_ENCSR718LCI" "h3k4me1_ms_nk_ENCSR806HUT" "h3k4me1_normal_b_ENCSR156BXM" "h3k4me1_normal_cd4_ENCSR102SOR" "h3k4me1_normal_cd8_ENCSR217SHH" "h3k4me1_normal_nk_ENCSR277YKG" "h3k4me3_ms_b_ENCSR260CRI" "h3k4me3_ms_b_ENCSR461QMZ" "h3k4me3_ms_b_ENCSR530AVK" "h3k4me3_ms_b_ENCSR713ZYF" "h3k4me3_ms_b_ENCSR923EGO" "h3k4me3_ms_cd4_ENCSR180NCM" "h3k4me3_ms_cd4_ENCSR341QLC" "h3k4me3_ms_cd4_ENCSR482TGI" "h3k4me3_ms_cd4_ENCSR486XJK" "h3k4me3_ms_cd4_ENCSR496LKR" "h3k4me3_ms_cd4_ENCSR603LTN" "h3k4me3_ms_cd4_ENCSR802MXQ" "h3k4me3_ms_cd4_ENCSR878YHM" "h3k4me3_ms_cd4_ENCSR954ZLD" "h3k4me3_ms_cd8_ENCSR231ZZH" "h3k4me3_ms_cd8_ENCSR278QHR" "h3k4me3_ms_cd8_ENCSR516CKJ" "h3k4me3_ms_cd8_ENCSR535YYH" "h3k4me3_ms_cd8_ENCSR741XAE" "h3k4me3_ms_cd8_ENCSR848XJL" "h3k4me3_ms_nk_ENCSR234BAD" "h3k4me3_ms_nk_ENCSR703TWY" "h3k4me3_ms_nk_ENCSR755ZGS" "h3k4me3_normal_b_ENCSR791CAF" "h3k4me3_normal_cd4_ENCSR537KJA" "h3k4me3_normal_cd8_ENCSR123ZAT" "h3k4me3_normal_nk_ENCSR394JFQ" "h3k9me3_ms_b_ENCSR445DNH" "h3k9me3_ms_b_ENCSR486SZN" "h3k9me3_ms_b_ENCSR983MKX" "h3k9me3_ms_cd4_ENCSR057IZD" "h3k9me3_ms_cd4_ENCSR550DPT" "h3k9me3_ms_cd4_ENCSR677OEF" "h3k9me3_ms_cd4_ENCSR729ZXH" "h3k9me3_ms_cd4_ENCSR851RJV" "h3k9me3_ms_cd4_ENCSR919DFZ" "h3k9me3_ms_cd4_ENCSR953STZ" "h3k9me3_ms_cd4_ENCSR959VZU" "h3k9me3_ms_cd8_ENCSR101USF" "h3k9me3_ms_cd8_ENCSR354GNT" "h3k9me3_ms_cd8_ENCSR377QYB" "h3k9me3_ms_cd8_ENCSR733LCG" "h3k9me3_ms_cd8_ENCSR980KTW" "h3k9me3_ms_nk_ENCSR061ATV" "h3k9me3_ms_nk_ENCSR611YIA" "h3k9me3_ms_nk_ENCSR667HUI" "h3k9me3_normal_b_ENCSR445LTM" "h3k9me3_normal_cd4_ENCSR433EWI" "h3k9me3_normal_cd8_ENCSR294HTM" "h3k9me3_normal_nk_ENCSR025UNZ")

 Process each sample folder
for sample_id in "${sample_ids[@]}"; do
    input_folder="/cta/users/merve.kaftancioglu/alternative_scenario/seqs/$sample_id"
    output_folder="/cta/users/merve.kaftancioglu/alternative_scenario/trimmomatic_as/output/$sample_id"

     Create output folder if it doesn't exist
    mkdir -p "$output_folder"

     Process all .fastq.gz files in the specified folder
    for input_file in "$input_folder"/*.fastq.gz; do
        output_file="$output_folder/output_paired_$(basename "$input_file")"
        trimmomatic SE -threads 16 -phred33 -trimlog "$output_folder/trimmomatic_log.txt" "$input_file" "$output_file" HEADCROP:14
    done
done

""" !/bin/bash
SBATCH --job-name=trimmomatic_single
SBATCH --nodes=1
SBATCH --account=users
SBATCH --qos=mid_mdbf
SBATCH --partition=mid_mdbf
SBATCH --cpus-per-task=16
SBATCH -o /cta/users/merve.kaftancioglu/alternative_scenario/trimmomatic_as/out-err/%x-%j-slurm.out
SBATCH -e /cta/users/merve.kaftancioglu/alternative_scenario/trimmomatic_as/out-err/%x-%j-slurm.err

 Load FastQC module
module load trimmomatic-0.38-gcc-9.2.0-3evdnxc

 List of sample IDs (folders)
sample_ids=("h3k27ac_ms_cd4_ENCSR331WMS" "h3k27ac_ms_cd4_ENCSR350UKV" "h3k27ac_ms_cd4_ENCSR474PYR")

 Process each sample folder
for sample_id in "${sample_ids[@]}"; do
    input_folder="/cta/users/merve.kaftancioglu/alternative_scenario/seqs/$sample_id"
    output_folder="/cta/users/merve.kaftancioglu/alternative_scenario/trimmomatic_as/output/$sample_id"

     Create output folder if it doesn't exist
    mkdir -p "$output_folder"

     Process all .fastq.gz files in the specified folder
    for input_file in "$input_folder"/*.fastq.gz; do
        output_file="$output_folder/output_paired_$(basename "$input_file")"
        trimmomatic SE -threads 16 -phred33 -trimlog "$output_folder/trimmomatic_log.txt" "$input_file" "$output_file" HEADCROP:14
    done
done

""" !/bin/bash
SBATCH --job-name=trimmomatic_single
SBATCH --nodes=1
SBATCH --account=users
SBATCH --qos=mid_mdbf
SBATCH --partition=mid_mdbf
SBATCH --cpus-per-task=16
SBATCH -o /cta/users/merve.kaftancioglu/alternative_scenario/trimmomatic_as/out-err/%x-%j-slurm.out
SBATCH -e /cta/users/merve.kaftancioglu/alternative_scenario/trimmomatic_as/out-err/%x-%j-slurm.err

 Load FastQC module
module load trimmomatic-0.38-gcc-9.2.0-3evdnxc

 List of sample IDs (folders)
sample_ids=("h3k27ac_ms_cd4_ENCSR520QDR" "h3k27ac_ms_cd4_ENCSR540XNK" "h3k27ac_ms_cd4_ENCSR705VSO")

 Process each sample folder
for sample_id in "${sample_ids[@]}"; do
    input_folder="/cta/users/merve.kaftancioglu/alternative_scenario/seqs/$sample_id"
    output_folder="/cta/users/merve.kaftancioglu/alternative_scenario/trimmomatic_as/output/$sample_id"

     Create output folder if it doesn't exist
    mkdir -p "$output_folder"

     Process all .fastq.gz files in the specified folder
    for input_file in "$input_folder"/*.fastq.gz; do
        output_file="$output_folder/output_paired_$(basename "$input_file")"
        trimmomatic SE -threads 16 -phred33 -trimlog "$output_folder/trimmomatic_log.txt" "$input_file" "$output_file" HEADCROP:14
    done
done

""" !/bin/bash
SBATCH --job-name=trimmomatic_single
SBATCH --nodes=1
SBATCH --account=users
SBATCH --qos=mid_mdbf
SBATCH --partition=mid_mdbf
SBATCH --cpus-per-task=16
SBATCH -o /cta/users/merve.kaftancioglu/alternative_scenario/trimmomatic_as/out-err/%x-%j-slurm.out
SBATCH -e /cta/users/merve.kaftancioglu/alternative_scenario/trimmomatic_as/out-err/%x-%j-slurm.err

 Load FastQC module
module load trimmomatic-0.38-gcc-9.2.0-3evdnxc

 List of sample IDs (folders)
sample_ids=("h3k27ac_ms_cd8_ENCSR458TOW" "h3k27ac_ms_cd8_ENCSR476IPR" "h3k27ac_ms_cd8_ENCSR787HDF")

 Process each sample folder
for sample_id in "${sample_ids[@]}"; do
    input_folder="/cta/users/merve.kaftancioglu/alternative_scenario/seqs/$sample_id"
    output_folder="/cta/users/merve.kaftancioglu/alternative_scenario/trimmomatic_as/output/$sample_id"

     Create output folder if it doesn't exist
    mkdir -p "$output_folder"

     Process all .fastq.gz files in the specified folder
    for input_file in "$input_folder"/*.fastq.gz; do
        output_file="$output_folder/output_paired_$(basename "$input_file")"
        trimmomatic SE -threads 16 -phred33 -trimlog "$output_folder/trimmomatic_log.txt" "$input_file" "$output_file" HEADCROP:14
    done
done

""" !/bin/bash
SBATCH --job-name=trimmomatic_single
SBATCH --nodes=1
SBATCH --account=users
SBATCH --qos=mid_mdbf
SBATCH --partition=mid_mdbf
SBATCH --cpus-per-task=16
SBATCH -o /cta/users/merve.kaftancioglu/alternative_scenario/trimmomatic_as/out-err/%x-%j-slurm.out
SBATCH -e /cta/users/merve.kaftancioglu/alternative_scenario/trimmomatic_as/out-err/%x-%j-slurm.err

 Load FastQC module
module load trimmomatic-0.38-gcc-9.2.0-3evdnxc

 List of sample IDs (folders)
sample_ids=("h3k27ac_ms_cd8_ENCSR923JIB" "h3k27ac_ms_nk_ENCSR246TTM" "h3k27ac_ms_nk_ENCSR469CFL")

 Process each sample folder
for sample_id in "${sample_ids[@]}"; do
    input_folder="/cta/users/merve.kaftancioglu/alternative_scenario/seqs/$sample_id"
    output_folder="/cta/users/merve.kaftancioglu/alternative_scenario/trimmomatic_as/output/$sample_id"

     Create output folder if it doesn't exist
    mkdir -p "$output_folder"

     Process all .fastq.gz files in the specified folder
    for input_file in "$input_folder"/*.fastq.gz; do
        output_file="$output_folder/output_paired_$(basename "$input_file")"
        trimmomatic SE -threads 16 -phred33 -trimlog "$output_folder/trimmomatic_log.txt" "$input_file" "$output_file" HEADCROP:14
    done
done

""" !/bin/bash
SBATCH --job-name=trimmomatic_single
SBATCH --nodes=1
SBATCH --account=users
SBATCH --qos=mid_mdbf
SBATCH --partition=mid_mdbf
SBATCH --cpus-per-task=16
SBATCH -o /cta/users/merve.kaftancioglu/alternative_scenario/trimmomatic_as/out-err/%x-%j-slurm.out
SBATCH -e /cta/users/merve.kaftancioglu/alternative_scenario/trimmomatic_as/out-err/%x-%j-slurm.err

 Load FastQC module
module load trimmomatic-0.38-gcc-9.2.0-3evdnxc

 List of sample IDs (folders)
sample_ids=("h3k27ac_ms_nk_ENCSR746AIX" "h3k27ac_normal_b_ENCSR685KZA" "h3k27ac_normal_cd4_ENCSR819NCZ")

 Process each sample folder
for sample_id in "${sample_ids[@]}"; do
    input_folder="/cta/users/merve.kaftancioglu/alternative_scenario/seqs/$sample_id"
    output_folder="/cta/users/merve.kaftancioglu/alternative_scenario/trimmomatic_as/output/$sample_id"

     Create output folder if it doesn't exist
    mkdir -p "$output_folder"

     Process all .fastq.gz files in the specified folder
    for input_file in "$input_folder"/*.fastq.gz; do
        output_file="$output_folder/output_paired_$(basename "$input_file")"
        trimmomatic SE -threads 16 -phred33 -trimlog "$output_folder/trimmomatic_log.txt" "$input_file" "$output_file" HEADCROP:14
    done
done

""" !/bin/bash
SBATCH --job-name=trimmomatic_single
SBATCH --nodes=1
SBATCH --account=users
SBATCH --qos=mid_mdbf
SBATCH --partition=mid_mdbf
SBATCH --cpus-per-task=16
SBATCH -o /cta/users/merve.kaftancioglu/alternative_scenario/trimmomatic_as/out-err/%x-%j-slurm.out
SBATCH -e /cta/users/merve.kaftancioglu/alternative_scenario/trimmomatic_as/out-err/%x-%j-slurm.err

 Load FastQC module
module load trimmomatic-0.38-gcc-9.2.0-3evdnxc

 List of sample IDs (folders)
sample_ids=("h3k27ac_normal_cd8_ENCSR976RWL" "h3k27ac_normal_nk_ENCSR977FMZ" "h3k27me3_ms_b_ENCSR009ZRH")

 Process each sample folder
for sample_id in "${sample_ids[@]}"; do
    input_folder="/cta/users/merve.kaftancioglu/alternative_scenario/seqs/$sample_id"
    output_folder="/cta/users/merve.kaftancioglu/alternative_scenario/trimmomatic_as/output/$sample_id"

     Create output folder if it doesn't exist
    mkdir -p "$output_folder"

     Process all .fastq.gz files in the specified folder
    for input_file in "$input_folder"/*.fastq.gz; do
        output_file="$output_folder/output_paired_$(basename "$input_file")"
        trimmomatic SE -threads 16 -phred33 -trimlog "$output_folder/trimmomatic_log.txt" "$input_file" "$output_file" HEADCROP:14
    done
done

""" !/bin/bash
SBATCH --job-name=trimmomatic_single
SBATCH --nodes=1
SBATCH --account=users
SBATCH --qos=mid_mdbf
SBATCH --partition=mid_mdbf
SBATCH --cpus-per-task=16
SBATCH -o /cta/users/merve.kaftancioglu/alternative_scenario/trimmomatic_as/out-err/%x-%j-slurm.out
SBATCH -e /cta/users/merve.kaftancioglu/alternative_scenario/trimmomatic_as/out-err/%x-%j-slurm.err

 Load FastQC module
module load trimmomatic-0.38-gcc-9.2.0-3evdnxc

 List of sample IDs (folders)
sample_ids=("h3k27me3_ms_b_ENCSR182NLA" "h3k27me3_ms_b_ENCSR272YVX" "h3k27me3_ms_b_ENCSR649FUX")

 Process each sample folder
for sample_id in "${sample_ids[@]}"; do
    input_folder="/cta/users/merve.kaftancioglu/alternative_scenario/seqs/$sample_id"
    output_folder="/cta/users/merve.kaftancioglu/alternative_scenario/trimmomatic_as/output/$sample_id"

     Create output folder if it doesn't exist
    mkdir -p "$output_folder"

     Process all .fastq.gz files in the specified folder
    for input_file in "$input_folder"/*.fastq.gz; do
        output_file="$output_folder/output_paired_$(basename "$input_file")"
        trimmomatic SE -threads 16 -phred33 -trimlog "$output_folder/trimmomatic_log.txt" "$input_file" "$output_file" HEADCROP:14
    done
done

""" !/bin/bash
SBATCH --job-name=fastQC
SBATCH --nodes=1
SBATCH --account=users
SBATCH --qos=cuda
SBATCH --partition=cuda
SBATCH --cpus-per-task=16
SBATCH -o /cta/users/merve.kaftancioglu/alternative_scenario/fastQC_as/out-err/%x-%j-slurm.out
SBATCH -e /cta/users/merve.kaftancioglu/alternative_scenario/fastQC_as/out-err/%x-%j-slurm.err


 Load FastQC module
module load fastqc-0.11.7-gcc-9.2.0-nsjzjof

substrate_ids=("h3k27ac_ms_b_ENCSR295QZX" "h3k27ac_ms_b_ENCSR364OIK" "h3k27ac_ms_b_ENCSR538URI" "h3k27ac_ms_b_ENCSR541PLO" "h3k27ac_ms_b_ENCSR617GMR" "h3k27ac_ms_cd4_ENCSR200SSJ" "h3k27ac_ms_cd4_ENCSR322MTA" "h3k27ac_ms_cd4_ENCSR331WMS" "h3k27ac_ms_cd4_ENCSR350UKV" "h3k27ac_ms_cd4_ENCSR474PYR" "h3k27ac_ms_cd4_ENCSR520QDR" "h3k27ac_ms_cd4_ENCSR540XNK" "h3k27ac_ms_cd4_ENCSR705VSO" "h3k27ac_ms_cd4_ENCSR832UMM" "h3k27ac_ms_cd8_ENCSR078ATS" "h3k27ac_ms_cd8_ENCSR348YRH" "h3k27ac_ms_cd8_ENCSR458TOW" "h3k27ac_ms_cd8_ENCSR476IPR" "h3k27ac_ms_cd8_ENCSR787HDF" "h3k27ac_ms_cd8_ENCSR923JIB" "h3k27ac_ms_nk_ENCSR246TTM" "h3k27ac_ms_nk_ENCSR469CFL" "h3k27ac_ms_nk_ENCSR746AIX" "h3k27ac_normal_b_ENCSR685KZA" "h3k27ac_normal_cd4_ENCSR819NCZ" "h3k27ac_normal_cd8_ENCSR976RWL" "h3k27ac_normal_nk_ENCSR977FMZ" "h3k27me3_ms_b_ENCSR009ZRH" "h3k27me3_ms_b_ENCSR182NLA" "h3k27me3_ms_b_ENCSR272YVX" "h3k27me3_ms_b_ENCSR649FUX" "h3k27me3_ms_b_ENCSR842WWX" "h3k27me3_ms_cd4_ENCSR277XYX" "h3k27me3_ms_cd4_ENCSR526TNC" "h3k27me3_ms_cd4_ENCSR592EKF" "h3k27me3_ms_cd4_ENCSR613UFD" "h3k27me3_ms_cd4_ENCSR740SDR" "h3k27me3_ms_cd4_ENCSR779JLY" "h3k27me3_ms_cd4_ENCSR993CTA" "h3k27me3_ms_cd8_ENCSR116FVG" "h3k27me3_ms_cd8_ENCSR122JCM" "h3k27me3_ms_cd8_ENCSR216ZVA" "h3k27me3_ms_cd8_ENCSR284IKS" "h3k27me3_ms_cd8_ENCSR385BOZ" "h3k27me3_ms_cd8_ENCSR521SFR" "h3k27me3_ms_nk_ENCSR469QVG" "h3k27me3_ms_nk_ENCSR565WDW" "h3k27me3_normal_b_ENCSR589LHR" "h3k27me3_normal_cd4_ENCSR068YVZ" "h3k27me3_normal_cd8_ENCSR720QRX" "h3k27me3_normal_nk_ENCSR639NIG" "h3k36me3_ms_b_ENCSR089VQL" "h3k36me3_ms_b_ENCSR119PSR" "h3k36me3_ms_b_ENCSR238WFK" "h3k36me3_ms_b_ENCSR987OPY" "h3k36me3_ms_cd4_ENCSR276NGH" "h3k36me3_ms_cd4_ENCSR330CQU" "h3k36me3_ms_cd4_ENCSR471WZE" "h3k36me3_ms_cd4_ENCSR473TGK" "h3k36me3_ms_cd4_ENCSR482VIB" "h3k36me3_ms_cd4_ENCSR532PXR" "h3k36me3_ms_cd4_ENCSR785PKM" "h3k36me3_ms_cd4_ENCSR865FTW" "h3k36me3_ms_cd4_ENCSR898VJE" "h3k36me3_ms_cd8_ENCSR239OWL" "h3k36me3_ms_cd8_ENCSR435MSB" "h3k36me3_ms_cd8_ENCSR631VIW" "h3k36me3_ms_cd8_ENCSR652WFO" "h3k36me3_ms_cd8_ENCSR757FGN" "h3k36me3_ms_nk_ENCSR158VSE" "h3k36me3_ms_nk_ENCSR245KON" "h3k36me3_ms_nk_ENCSR530YDY" "h3k36me3_normal_b_ENCSR831AXK" "h3k36me3_normal_cd4_ENCSR640OKB" "h3k36me3_normal_cd8_ENCSR303SQG" "h3k36me3_normal_nk_ENCSR056GJY" "h3k4me1_ms_b_ENCSR084FXT" "h3k4me1_ms_b_ENCSR480ITK" "h3k4me1_ms_b_ENCSR759DHN" "h3k4me1_ms_b_ENCSR900IAM" "h3k4me1_ms_b_ENCSR931JYC" "h3k4me1_ms_cd4_ENCSR036YVP" "h3k4me1_ms_cd4_ENCSR043JIA" "h3k4me1_ms_cd4_ENCSR093VUP" "h3k4me1_ms_cd4_ENCSR124IGC" "h3k4me1_ms_cd4_ENCSR315MYP" "h3k4me1_ms_cd4_ENCSR458QEY" "h3k4me1_ms_cd4_ENCSR485FBT" "h3k4me1_ms_cd4_ENCSR641ZFV" "h3k4me1_ms_cd4_ENCSR687CQX" "h3k4me1_ms_cd8_ENCSR231XAP" "h3k4me1_ms_cd8_ENCSR572XTB" "h3k4me1_ms_cd8_ENCSR788TEF" "h3k4me1_ms_cd8_ENCSR815YZL" "h3k4me1_ms_cd8_ENCSR861FEC" "h3k4me1_ms_cd8_ENCSR940PHE" "h3k4me1_ms_nk_ENCSR482UGV" "h3k4me1_ms_nk_ENCSR718LCI" "h3k4me1_ms_nk_ENCSR806HUT" "h3k4me1_normal_b_ENCSR156BXM" "h3k4me1_normal_cd4_ENCSR102SOR" "h3k4me1_normal_cd8_ENCSR217SHH" "h3k4me1_normal_nk_ENCSR277YKG" "h3k4me3_ms_b_ENCSR260CRI" "h3k4me3_ms_b_ENCSR461QMZ" "h3k4me3_ms_b_ENCSR530AVK" "h3k4me3_ms_b_ENCSR713ZYF" "h3k4me3_ms_b_ENCSR923EGO" "h3k4me3_ms_cd4_ENCSR180NCM" "h3k4me3_ms_cd4_ENCSR341QLC" "h3k4me3_ms_cd4_ENCSR482TGI" "h3k4me3_ms_cd4_ENCSR486XJK" "h3k4me3_ms_cd4_ENCSR496LKR" "h3k4me3_ms_cd4_ENCSR603LTN" "h3k4me3_ms_cd4_ENCSR802MXQ" "h3k4me3_ms_cd4_ENCSR878YHM" "h3k4me3_ms_cd4_ENCSR954ZLD" "h3k4me3_ms_cd8_ENCSR231ZZH" "h3k4me3_ms_cd8_ENCSR278QHR" "h3k4me3_ms_cd8_ENCSR516CKJ" "h3k4me3_ms_cd8_ENCSR535YYH" "h3k4me3_ms_cd8_ENCSR741XAE" "h3k4me3_ms_cd8_ENCSR848XJL" "h3k4me3_ms_nk_ENCSR234BAD" "h3k4me3_ms_nk_ENCSR703TWY" "h3k4me3_ms_nk_ENCSR755ZGS" "h3k4me3_normal_b_ENCSR791CAF" "h3k4me3_normal_cd4_ENCSR537KJA" "h3k4me3_normal_cd8_ENCSR123ZAT" "h3k4me3_normal_nk_ENCSR394JFQ" "h3k9me3_ms_b_ENCSR445DNH" "h3k9me3_ms_b_ENCSR486SZN" "h3k9me3_ms_b_ENCSR983MKX" "h3k9me3_ms_cd4_ENCSR057IZD" "h3k9me3_ms_cd4_ENCSR550DPT" "h3k9me3_ms_cd4_ENCSR677OEF" "h3k9me3_ms_cd4_ENCSR729ZXH" "h3k9me3_ms_cd4_ENCSR851RJV" "h3k9me3_ms_cd4_ENCSR919DFZ" "h3k9me3_ms_cd4_ENCSR953STZ" "h3k9me3_ms_cd4_ENCSR959VZU" "h3k9me3_ms_cd8_ENCSR101USF" "h3k9me3_ms_cd8_ENCSR354GNT" "h3k9me3_ms_cd8_ENCSR377QYB" "h3k9me3_ms_cd8_ENCSR733LCG" "h3k9me3_ms_cd8_ENCSR980KTW" "h3k9me3_ms_nk_ENCSR061ATV" "h3k9me3_ms_nk_ENCSR611YIA" "h3k9me3_ms_nk_ENCSR667HUI" "h3k9me3_normal_b_ENCSR445LTM" "h3k9me3_normal_cd4_ENCSR433EWI" "h3k9me3_normal_cd8_ENCSR294HTM" "h3k9me3_normal_nk_ENCSR025UNZ")

for substrate_id in "${substrate_ids[@]}"; do
    echo "Predicting ${substrate_id}..."
    input_directory="/cta/users/merve.kaftancioglu/alternative_scenario/seqs/${substrate_id}"
    output_dir="/cta/users/merve.kaftancioglu/alternative_scenario/fastQC_as/output/${substrate_id}"    
    mkdir -p "${output_dir}"
    fastqc --outdir "$output_directory" "$input_directory"/*.fastq.gz -threads 16
    echo "Finished ${substrate_id}"
done


echo "DONE!""" "



""" !/bin/bash
SBATCH --job-name=conversion
SBATCH --nodes=1
SBATCH --account=users
SBATCH --qos=cuda
SBATCH --partition=cuda
SBATCH --cpus-per-task=16
SBATCH -o /cta/users/merve.kaftancioglu/alternative_scenario/bowtie2_as/output/out-err/%x-%j-slurm.out
SBATCH -e /cta/users/merve.kaftancioglu/alternative_scenario/bowtie2_as/output/out-err/%x-%j-slurm.err

 Set directory containing fastq files
fastq_dir="/cta/users/merve.kaftancioglu/alternative_scenario/seqs/h3k27ac_normal_cd4_ENCSR819NCZ"

 Set output directory (optional, use same directory by default)
fasta_dir="/cta/users/merve.kaftancioglu/alternative_scenario/seqs/h3k27ac_normal_cd4_ENCSR819NCZ_fasta"   Uncomment if you want a separate directory

 Loop through all fastq files in the directory
for fastq_file in "$fastq_dir"/*.fastq.gz; do
   Extract filename without extension
  fasta_filename="${fastq_file%.fastq.gz}.fasta"

   Convert fastq.gz to fasta using seqtk
  seqtk seq -a "$fastq_file" > "$fasta_filename"

  echo "Converted $fastq_file to $fasta_filename"
done

echo "All fastq files converted to fasta!"

""" !/bin/bash
SBATCH --job-name=trimmomatic_single
SBATCH --nodes=1
SBATCH --account=users
SBATCH --qos=mid_mdbf
SBATCH --partition=mid_mdbf
SBATCH --cpus-per-task=16
SBATCH -o /cta/users/merve.kaftancioglu/alternative_scenario/trimmomatic_as/out-err/%x-%j-slurm.out
SBATCH -e /cta/users/merve.kaftancioglu/alternative_scenario/trimmomatic_as/out-err/%x-%j-slurm.err

 Load FastQC module
module load trimmomatic-0.38-gcc-9.2.0-3evdnxc

 List of sample IDs (folders)
sample_ids=("h3k27me3_ms_b_ENCSR842WWX" "h3k27me3_ms_cd4_ENCSR277XYX" "h3k27me3_ms_cd4_ENCSR526TNC")

 Process each sample folder
for sample_id in "${sample_ids[@]}"; do
    input_folder="/cta/users/merve.kaftancioglu/alternative_scenario/seqs/$sample_id"
    output_folder="/cta/users/merve.kaftancioglu/alternative_scenario/trimmomatic_as/output/$sample_id"

     Create output folder if it doesn't exist
    mkdir -p "$output_folder"

     Process all .fastq.gz files in the specified folder
    for input_file in "$input_folder"/*.fastq.gz; do
        output_file="$output_folder/output_paired_$(basename "$input_file")"
        trimmomatic SE -threads 16 -phred33 -trimlog "$output_folder/trimmomatic_log.txt" "$input_file" "$output_file" HEADCROP:14
    done
done

""" !/bin/bash
SBATCH --job-name=trimmomatic_single
SBATCH --nodes=1
SBATCH --account=users
SBATCH --qos=mid_mdbf
SBATCH --partition=mid_mdbf
SBATCH --cpus-per-task=16
SBATCH -o /cta/users/merve.kaftancioglu/alternative_scenario/trimmomatic_as/out-err/%x-%j-slurm.out
SBATCH -e /cta/users/merve.kaftancioglu/alternative_scenario/trimmomatic_as/out-err/%x-%j-slurm.err

 Load FastQC module
module load trimmomatic-0.38-gcc-9.2.0-3evdnxc

 List of sample IDs (folders)
sample_ids=("h3k27me3_ms_cd4_ENCSR592EKF" "h3k27me3_ms_cd4_ENCSR613UFD" "h3k27me3_ms_cd4_ENCSR740SDR")

 Process each sample folder
for sample_id in "${sample_ids[@]}"; do
    input_folder="/cta/users/merve.kaftancioglu/alternative_scenario/seqs/$sample_id"
    output_folder="/cta/users/merve.kaftancioglu/alternative_scenario/trimmomatic_as/output/$sample_id"

     Create output folder if it doesn't exist
    mkdir -p "$output_folder"

     Process all .fastq.gz files in the specified folder
    for input_file in "$input_folder"/*.fastq.gz; do
        output_file="$output_folder/output_paired_$(basename "$input_file")"
        trimmomatic SE -threads 16 -phred33 -trimlog "$output_folder/trimmomatic_log.txt" "$input_file" "$output_file" HEADCROP:14
    done
done

""" !/bin/bash
SBATCH --job-name=trimmomatic_single
SBATCH --nodes=1
SBATCH --account=users
SBATCH --qos=mid_mdbf
SBATCH --partition=mid_mdbf
SBATCH --cpus-per-task=16
SBATCH -o /cta/users/merve.kaftancioglu/alternative_scenario/trimmomatic_as/out-err/%x-%j-slurm.out
SBATCH -e /cta/users/merve.kaftancioglu/alternative_scenario/trimmomatic_as/out-err/%x-%j-slurm.err

 Load FastQC module
module load trimmomatic-0.38-gcc-9.2.0-3evdnxc

 List of sample IDs (folders)
sample_ids=("h3k27me3_ms_cd4_ENCSR779JLY" "h3k27me3_ms_cd4_ENCSR993CTA" "h3k27me3_ms_cd8_ENCSR116FVG")

 Process each sample folder
for sample_id in "${sample_ids[@]}"; do
    input_folder="/cta/users/merve.kaftancioglu/alternative_scenario/seqs/$sample_id"
    output_folder="/cta/users/merve.kaftancioglu/alternative_scenario/trimmomatic_as/output/$sample_id"

     Create output folder if it doesn't exist
    mkdir -p "$output_folder"

     Process all .fastq.gz files in the specified folder
    for input_file in "$input_folder"/*.fastq.gz; do
        output_file="$output_folder/output_paired_$(basename "$input_file")"
        trimmomatic SE -threads 16 -phred33 -trimlog "$output_folder/trimmomatic_log.txt" "$input_file" "$output_file" HEADCROP:14
    done
done

""" !/bin/bash
SBATCH --job-name=trimmomatic_single
SBATCH --nodes=1
SBATCH --account=users
SBATCH --qos=mid_mdbf
SBATCH --partition=mid_mdbf
SBATCH --cpus-per-task=16
SBATCH -o /cta/users/merve.kaftancioglu/alternative_scenario/trimmomatic_as/out-err/%x-%j-slurm.out
SBATCH -e /cta/users/merve.kaftancioglu/alternative_scenario/trimmomatic_as/out-err/%x-%j-slurm.err

 Load FastQC module
module load trimmomatic-0.38-gcc-9.2.0-3evdnxc

 List of sample IDs (folders)
sample_ids=("h3k27me3_ms_cd8_ENCSR122JCM" "h3k27me3_ms_cd8_ENCSR216ZVA" "h3k27me3_ms_cd8_ENCSR284IKS")

 Process each sample folder
for sample_id in "${sample_ids[@]}"; do
    input_folder="/cta/users/merve.kaftancioglu/alternative_scenario/seqs/$sample_id"
    output_folder="/cta/users/merve.kaftancioglu/alternative_scenario/trimmomatic_as/output/$sample_id"

     Create output folder if it doesn't exist
    mkdir -p "$output_folder"

     Process all .fastq.gz files in the specified folder
    for input_file in "$input_folder"/*.fastq.gz; do
        output_file="$output_folder/output_paired_$(basename "$input_file")"
        trimmomatic SE -threads 16 -phred33 -trimlog "$output_folder/trimmomatic_log.txt" "$input_file" "$output_file" HEADCROP:14
    done
done

""" !/bin/bash
SBATCH --job-name=trimmomatic_single
SBATCH --nodes=1
SBATCH --account=users
SBATCH --qos=mid_mdbf
SBATCH --partition=mid_mdbf
SBATCH --cpus-per-task=16
SBATCH -o /cta/users/merve.kaftancioglu/alternative_scenario/trimmomatic_as/out-err/%x-%j-slurm.out
SBATCH -e /cta/users/merve.kaftancioglu/alternative_scenario/trimmomatic_as/out-err/%x-%j-slurm.err

 Load FastQC module
module load trimmomatic-0.38-gcc-9.2.0-3evdnxc

 List of sample IDs (folders)
sample_ids=("h3k27me3_ms_cd8_ENCSR385BOZ" "h3k27me3_ms_cd8_ENCSR521SFR" "h3k27me3_ms_nk_ENCSR469QVG")

 Process each sample folder
for sample_id in "${sample_ids[@]}"; do
    input_folder="/cta/users/merve.kaftancioglu/alternative_scenario/seqs/$sample_id"
    output_folder="/cta/users/merve.kaftancioglu/alternative_scenario/trimmomatic_as/output/$sample_id"

     Create output folder if it doesn't exist
    mkdir -p "$output_folder"

     Process all .fastq.gz files in the specified folder
    for input_file in "$input_folder"/*.fastq.gz; do
        output_file="$output_folder/output_paired_$(basename "$input_file")"
        trimmomatic SE -threads 16 -phred33 -trimlog "$output_folder/trimmomatic_log.txt" "$input_file" "$output_file" HEADCROP:14
    done
done

""" !/bin/bash
SBATCH --job-name=trimmomatic_single
SBATCH --nodes=1
SBATCH --account=users
SBATCH --qos=mid_mdbf
SBATCH --partition=mid_mdbf
SBATCH --cpus-per-task=16
SBATCH -o /cta/users/merve.kaftancioglu/alternative_scenario/trimmomatic_as/out-err/%x-%j-slurm.out
SBATCH -e /cta/users/merve.kaftancioglu/alternative_scenario/trimmomatic_as/out-err/%x-%j-slurm.err

 Load FastQC module
module load trimmomatic-0.38-gcc-9.2.0-3evdnxc

 List of sample IDs (folders)
sample_ids=("h3k27me3_ms_nk_ENCSR565WDW" "h3k27me3_normal_b_ENCSR589LHR" "h3k27me3_normal_cd4_ENCSR068YVZ")

 Process each sample folder
for sample_id in "${sample_ids[@]}"; do
    input_folder="/cta/users/merve.kaftancioglu/alternative_scenario/seqs/$sample_id"
    output_folder="/cta/users/merve.kaftancioglu/alternative_scenario/trimmomatic_as/output/$sample_id"

     Create output folder if it doesn't exist
    mkdir -p "$output_folder"

     Process all .fastq.gz files in the specified folder
    for input_file in "$input_folder"/*.fastq.gz; do
        output_file="$output_folder/output_paired_$(basename "$input_file")"
        trimmomatic SE -threads 16 -phred33 -trimlog "$output_folder/trimmomatic_log.txt" "$input_file" "$output_file" HEADCROP:14
    done
done

""" !/bin/bash
SBATCH --job-name=trimmomatic_single
SBATCH --nodes=1
SBATCH --account=users
SBATCH --qos=mid_mdbf
SBATCH --partition=mid_mdbf
SBATCH --cpus-per-task=16
SBATCH -o /cta/users/merve.kaftancioglu/alternative_scenario/trimmomatic_as/out-err/%x-%j-slurm.out
SBATCH -e /cta/users/merve.kaftancioglu/alternative_scenario/trimmomatic_as/out-err/%x-%j-slurm.err

 Load FastQC module
module load trimmomatic-0.38-gcc-9.2.0-3evdnxc

 List of sample IDs (folders)
sample_ids=("h3k27me3_normal_cd8_ENCSR720QRX" "h3k27me3_normal_nk_ENCSR639NIG" "h3k36me3_ms_b_ENCSR089VQL")

 Process each sample folder
for sample_id in "${sample_ids[@]}"; do
    input_folder="/cta/users/merve.kaftancioglu/alternative_scenario/seqs/$sample_id"
    output_folder="/cta/users/merve.kaftancioglu/alternative_scenario/trimmomatic_as/output/$sample_id"

     Create output folder if it doesn't exist
    mkdir -p "$output_folder"

     Process all .fastq.gz files in the specified folder
    for input_file in "$input_folder"/*.fastq.gz; do
        output_file="$output_folder/output_paired_$(basename "$input_file")"
        trimmomatic SE -threads 16 -phred33 -trimlog "$output_folder/trimmomatic_log.txt" "$input_file" "$output_file" HEADCROP:14
    done
done

""" !/bin/bash
SBATCH --job-name=trimmomatic_single
SBATCH --nodes=1
SBATCH --account=users
SBATCH --qos=mid_mdbf
SBATCH --partition=mid_mdbf
SBATCH --cpus-per-task=16
SBATCH -o /cta/users/merve.kaftancioglu/alternative_scenario/trimmomatic_as/out-err/%x-%j-slurm.out
SBATCH -e /cta/users/merve.kaftancioglu/alternative_scenario/trimmomatic_as/out-err/%x-%j-slurm.err

 Load FastQC module
module load trimmomatic-0.38-gcc-9.2.0-3evdnxc

 List of sample IDs (folders)
sample_ids=("h3k36me3_ms_b_ENCSR119PSR" "h3k36me3_ms_b_ENCSR238WFK" "h3k36me3_ms_b_ENCSR987OPY")

 Process each sample folder
for sample_id in "${sample_ids[@]}"; do
    input_folder="/cta/users/merve.kaftancioglu/alternative_scenario/seqs/$sample_id"
    output_folder="/cta/users/merve.kaftancioglu/alternative_scenario/trimmomatic_as/output/$sample_id"

     Create output folder if it doesn't exist
    mkdir -p "$output_folder"

     Process all .fastq.gz files in the specified folder
    for input_file in "$input_folder"/*.fastq.gz; do
        output_file="$output_folder/output_paired_$(basename "$input_file")"
        trimmomatic SE -threads 16 -phred33 -trimlog "$output_folder/trimmomatic_log.txt" "$input_file" "$output_file" HEADCROP:14
    done
done

""" !/bin/bash
SBATCH --job-name=trimmomatic_single
SBATCH --nodes=1
SBATCH --account=users
SBATCH --qos=mid_mdbf
SBATCH --partition=mid_mdbf
SBATCH --cpus-per-task=16
SBATCH -o /cta/users/merve.kaftancioglu/alternative_scenario/trimmomatic_as/out-err/%x-%j-slurm.out
SBATCH -e /cta/users/merve.kaftancioglu/alternative_scenario/trimmomatic_as/out-err/%x-%j-slurm.err

 Load FastQC module
module load trimmomatic-0.38-gcc-9.2.0-3evdnxc

 List of sample IDs (folders)
sample_ids=( "h3k36me3_ms_cd4_ENCSR276NGH" "h3k36me3_ms_cd4_ENCSR330CQU" "h3k36me3_ms_cd4_ENCSR471WZE")

 Process each sample folder
for sample_id in "${sample_ids[@]}"; do
    input_folder="/cta/users/merve.kaftancioglu/alternative_scenario/seqs/$sample_id"
    output_folder="/cta/users/merve.kaftancioglu/alternative_scenario/trimmomatic_as/output/$sample_id"

     Create output folder if it doesn't exist
    mkdir -p "$output_folder"

     Process all .fastq.gz files in the specified folder
    for input_file in "$input_folder"/*.fastq.gz; do
        output_file="$output_folder/output_paired_$(basename "$input_file")"
        trimmomatic SE -threads 16 -phred33 -trimlog "$output_folder/trimmomatic_log.txt" "$input_file" "$output_file" HEADCROP:14
    done
done

""" !/bin/bash
SBATCH --job-name=trimmomatic_single
SBATCH --nodes=1
SBATCH --account=users
SBATCH --qos=mid_mdbf
SBATCH --partition=mid_mdbf
SBATCH --cpus-per-task=16
SBATCH -o /cta/users/merve.kaftancioglu/alternative_scenario/trimmomatic_as/out-err/%x-%j-slurm.out
SBATCH -e /cta/users/merve.kaftancioglu/alternative_scenario/trimmomatic_as/out-err/%x-%j-slurm.err

 Load FastQC module
module load trimmomatic-0.38-gcc-9.2.0-3evdnxc

 List of sample IDs (folders)
sample_ids=("h3k36me3_ms_cd4_ENCSR473TGK" "h3k36me3_ms_cd4_ENCSR482VIB" "h3k36me3_ms_cd4_ENCSR532PXR")

 Process each sample folder
for sample_id in "${sample_ids[@]}"; do
    input_folder="/cta/users/merve.kaftancioglu/alternative_scenario/seqs/$sample_id"
    output_folder="/cta/users/merve.kaftancioglu/alternative_scenario/trimmomatic_as/output/$sample_id"

     Create output folder if it doesn't exist
    mkdir -p "$output_folder"

     Process all .fastq.gz files in the specified folder
    for input_file in "$input_folder"/*.fastq.gz; do
        output_file="$output_folder/output_paired_$(basename "$input_file")"
        trimmomatic SE -threads 16 -phred33 -trimlog "$output_folder/trimmomatic_log.txt" "$input_file" "$output_file" HEADCROP:14
    done
done

""" !/bin/bash
SBATCH --job-name=trimmomatic_single
SBATCH --nodes=1
SBATCH --account=users
SBATCH --qos=mid_mdbf
SBATCH --partition=mid_mdbf
SBATCH --cpus-per-task=16
SBATCH -o /cta/users/merve.kaftancioglu/alternative_scenario/trimmomatic_as/out-err/%x-%j-slurm.out
SBATCH -e /cta/users/merve.kaftancioglu/alternative_scenario/trimmomatic_as/out-err/%x-%j-slurm.err

 Load FastQC module
module load trimmomatic-0.38-gcc-9.2.0-3evdnxc

 List of sample IDs (folders)
sample_ids=("h3k36me3_ms_cd4_ENCSR785PKM" "h3k36me3_ms_cd4_ENCSR865FTW" "h3k36me3_ms_cd4_ENCSR898VJE")

 Process each sample folder
for sample_id in "${sample_ids[@]}"; do
    input_folder="/cta/users/merve.kaftancioglu/alternative_scenario/seqs/$sample_id"
    output_folder="/cta/users/merve.kaftancioglu/alternative_scenario/trimmomatic_as/output/$sample_id"

     Create output folder if it doesn't exist
    mkdir -p "$output_folder"

     Process all .fastq.gz files in the specified folder
    for input_file in "$input_folder"/*.fastq.gz; do
        output_file="$output_folder/output_paired_$(basename "$input_file")"
        trimmomatic SE -threads 16 -phred33 -trimlog "$output_folder/trimmomatic_log.txt" "$input_file" "$output_file" HEADCROP:14
    done
done

""" !/bin/bash
SBATCH --job-name=trimmomatic_single
SBATCH --nodes=1
SBATCH --account=users
SBATCH --qos=mid_mdbf
SBATCH --partition=mid_mdbf
SBATCH --cpus-per-task=16
SBATCH -o /cta/users/merve.kaftancioglu/alternative_scenario/trimmomatic_as/out-err/%x-%j-slurm.out
SBATCH -e /cta/users/merve.kaftancioglu/alternative_scenario/trimmomatic_as/out-err/%x-%j-slurm.err

 Load FastQC module
module load trimmomatic-0.38-gcc-9.2.0-3evdnxc

 List of sample IDs (folders)
sample_ids=("h3k36me3_ms_cd8_ENCSR239OWL" "h3k36me3_ms_cd8_ENCSR435MSB" "h3k36me3_ms_cd8_ENCSR631VIW")

 Process each sample folder
for sample_id in "${sample_ids[@]}"; do
    input_folder="/cta/users/merve.kaftancioglu/alternative_scenario/seqs/$sample_id"
    output_folder="/cta/users/merve.kaftancioglu/alternative_scenario/trimmomatic_as/output/$sample_id"

     Create output folder if it doesn't exist
    mkdir -p "$output_folder"

     Process all .fastq.gz files in the specified folder
    for input_file in "$input_folder"/*.fastq.gz; do
        output_file="$output_folder/output_paired_$(basename "$input_file")"
        trimmomatic SE -threads 16 -phred33 -trimlog "$output_folder/trimmomatic_log.txt" "$input_file" "$output_file" HEADCROP:14
    done
done

""" !/bin/bash
SBATCH --job-name=download_processed
SBATCH --nodes=1
SBATCH --account=users
SBATCH --qos=mid_mdbf
SBATCH --partition=mid_mdbf
SBATCH --cpus-per-task=16
SBATCH -o /cta/users/merve.kaftancioglu/alternative_scenario/processed_seqs/out-err/%x-%j-slurm.out
SBATCH -e /cta/users/merve.kaftancioglu/alternative_scenario/processed_seqs/out-err/%x-%j-slurm.err

ENCSR348YRH="/cta/users/merve.kaftancioglu/alternative_scenario/processed_seqs/processed_ENCSR348YRH"

wget -P "processed_ENCSR348YRH" https://www.encodeproject.org/metadata/?type=Experiment&%40id=%2Fexperiments%2FENCSR348YRH%2F&files.analyses.%40id=%2Fanalyses%2FENCAN186KEO%2F
wget -P "processed_ENCSR348YRH" https://www.encodeproject.org/files/ENCFF415QNO/@@download/ENCFF415QNO.bigWig
wget -P "processed_ENCSR348YRH" https://www.encodeproject.org/files/ENCFF239RTD/@@download/ENCFF239RTD.bed.gz
wget -P "processed_ENCSR348YRH" https://www.encodeproject.org/files/ENCFF789UUS/@@download/ENCFF789UUS.bam
wget -P "processed_ENCSR348YRH" https://www.encodeproject.org/files/ENCFF101DXS/@@download/ENCFF101DXS.bigWig
wget -P "processed_ENCSR348YRH" https://www.encodeproject.org/files/ENCFF902JZF/@@download/ENCFF902JZF.bam
wget -P "processed_ENCSR348YRH" https://www.encodeproject.org/files/ENCFF158YHJ/@@download/ENCFF158YHJ.bigBed

ENCSR350UKV="/cta/users/merve.kaftancioglu/alternative_scenario/processed_seqs/processed_ENCSR350UKV"

wget -P "processed_ENCSR350UKV" https://www.encodeproject.org/metadata/?type=Experiment&%40id=%2Fexperiments%2FENCSR350UKV%2F&files.analyses.%40id=%2Fanalyses%2FENCAN847UZS%2F
wget -P "processed_ENCSR350UKV" https://www.encodeproject.org/files/ENCFF987SQJ/@@download/ENCFF987SQJ.bam
wget -P "processed_ENCSR350UKV" https://www.encodeproject.org/files/ENCFF274JBU/@@download/ENCFF274JBU.bigWig
wget -P "processed_ENCSR350UKV" https://www.encodeproject.org/files/ENCFF729RXH/@@download/ENCFF729RXH.bam
wget -P "processed_ENCSR350UKV" https://www.encodeproject.org/files/ENCFF718ZSS/@@download/ENCFF718ZSS.bigBed
wget -P "processed_ENCSR350UKV" https://www.encodeproject.org/files/ENCFF979YMP/@@download/ENCFF979YMP.bigWig
wget -P "processed_ENCSR350UKV" https://www.encodeproject.org/files/ENCFF329MIW/@@download/ENCFF329MIW.bed.gz

""" !/bin/bash
SBATCH --job-name=trimmomatic_single
SBATCH --nodes=1
SBATCH --account=users
SBATCH --qos=mid_mdbf
SBATCH --partition=mid_mdbf
SBATCH --cpus-per-task=16
SBATCH -o /cta/users/merve.kaftancioglu/alternative_scenario/trimmomatic_as/out-err/%x-%j-slurm.out
SBATCH -e /cta/users/merve.kaftancioglu/alternative_scenario/trimmomatic_as/out-err/%x-%j-slurm.err

 Load FastQC module
module load trimmomatic-0.38-gcc-9.2.0-3evdnxc

 List of sample IDs (folders)
sample_ids=("h3k36me3_ms_cd8_ENCSR239OWL" "h3k27me3_ms_b_ENCSR842WWX" "h3k27me3_ms_cd4_ENCSR277XYX" "h3k27me3_ms_cd4_ENCSR526TNC")

 Process each sample folder
for sample_id in "${sample_ids[@]}"; do
    input_folder="/cta/users/merve.kaftancioglu/alternative_scenario/seqs/$sample_id"
    output_folder="/cta/users/merve.kaftancioglu/alternative_scenario/trimmomatic_as/output/$sample_id"

     Create output folder if it doesn't exist
    mkdir -p "$output_folder"

     Process all .fastq.gz files in the specified folder
    for input_file in "$input_folder"/*.fastq.gz; do
        output_file="$output_folder/output_paired_$(basename "$input_file")"
        trimmomatic SE -threads 16 -phred33 -trimlog "$output_folder/trimmomatic_log.txt" "$input_file" "$output_file" HEADCROP:14
    done
done

""" !/bin/bash
SBATCH --job-name=trimmomatic_single
SBATCH --nodes=1
SBATCH --account=users
SBATCH --qos=mid_mdbf
SBATCH --partition=mid_mdbf
SBATCH --cpus-per-task=16
SBATCH -o /cta/users/merve.kaftancioglu/alternative_scenario/trimmomatic_as/out-err/%x-%j-slurm.out
SBATCH -e /cta/users/merve.kaftancioglu/alternative_scenario/trimmomatic_as/out-err/%x-%j-slurm.err

 Load FastQC module
module load trimmomatic-0.38-gcc-9.2.0-3evdnxc

 List of sample IDs (folders)
sample_ids=("h3k27ac_ms_b_ENCSR617GMR" "h3k27ac_ms_cd8_ENCSR348YRH" "h3k27me3_ms_b_ENCSR649FUX")

 Process each sample folder
for sample_id in "${sample_ids[@]}"; do
    input_folder="/cta/users/merve.kaftancioglu/alternative_scenario/seqs/$sample_id"
    output_folder="/cta/users/merve.kaftancioglu/alternative_scenario/trimmomatic_as/output/$sample_id"

     Create output folder if it doesn't exist
    mkdir -p "$output_folder"

     Process all .fastq.gz files in the specified folder
    for input_file in "$input_folder"/*.fastq.gz; do
        output_file="$output_folder/output_paired_$(basename "$input_file")"
        trimmomatic SE -threads 16 -phred33 -trimlog "$output_folder/trimmomatic_log.txt" "$input_file" "$output_file" HEADCROP:14
    done
done

""" !/bin/bash
SBATCH --job-name=download_processed
SBATCH --nodes=1
SBATCH --account=users
SBATCH --qos=short_mdbf
SBATCH --partition=short_mdbf
SBATCH --cpus-per-task=16
SBATCH -o /cta/users/merve.kaftancioglu/alternative_scenario/processed_seqs/out-err/%x-%j-slurm.out
SBATCH -e /cta/users/merve.kaftancioglu/alternative_scenario/processed_seqs/out-err/%x-%j-slurm.err

 Directory paths

ENCSR068YVZ="/cta/users/merve.kaftancioglu/alternative_scenario/processed_seqs/processed_ENCSR068YVZ"
ENCSR078ATS="/cta/users/merve.kaftancioglu/alternative_scenario/processed_seqs/processed_ENCSR078ATS"
ENCSR084FXT="/cta/users/merve.kaftancioglu/alternative_scenario/processed_seqs/processed_ENCSR084FXT"

ENCSR089VQL="/cta/users/merve.kaftancioglu/alternative_scenario/processed_seqs/processed_ENCSR089VQL"
ENCSR093VUP="/cta/users/merve.kaftancioglu/alternative_scenario/processed_seqs/processed_ENCSR093VUP"
ENCSR101USF="/cta/users/merve.kaftancioglu/alternative_scenario/processed_seqs/processed_ENCSR101USF"
ENCSR102SOR="/cta/users/merve.kaftancioglu/alternative_scenario/processed_seqs/processed_ENCSR102SOR"
ENCSR116FVG="/cta/users/merve.kaftancioglu/alternative_scenario/processed_seqs/processed_ENCSR116FVG"

ENCSR119PSR="/cta/users/merve.kaftancioglu/alternative_scenario/processed_seqs/processed_ENCSR119PSR"
ENCSR122JCM="/cta/users/merve.kaftancioglu/alternative_scenario/processed_seqs/processed_ENCSR122JCM"
ENCSR123ZAT="/cta/users/merve.kaftancioglu/alternative_scenario/processed_seqs/processed_ENCSR123ZAT"
ENCSR124IGC="/cta/users/merve.kaftancioglu/alternative_scenario/processed_seqs/processed_ENCSR124IGC"
ENCSR156BXM="/cta/users/merve.kaftancioglu/alternative_scenario/processed_seqs/processed_ENCSR156BXM"

ENCSR158VSE="/cta/users/merve.kaftancioglu/alternative_scenario/processed_seqs/processed_ENCSR158VSE"
ENCSR180NCM="/cta/users/merve.kaftancioglu/alternative_scenario/processed_seqs/processed_ENCSR180NCM"
ENCSR182NLA="/cta/users/merve.kaftancioglu/alternative_scenario/processed_seqs/processed_ENCSR182NLA"
ENCSR200SSJ="/cta/users/merve.kaftancioglu/alternative_scenario/processed_seqs/processed_ENCSR200SSJ"
ENCSR216ZVA="/cta/users/merve.kaftancioglu/alternative_scenario/processed_seqs/processed_ENCSR216ZVA"

ENCSR217SHH="/cta/users/merve.kaftancioglu/alternative_scenario/processed_seqs/processed_ENCSR217SHH"
ENCSR231XAP="/cta/users/merve.kaftancioglu/alternative_scenario/processed_seqs/processed_ENCSR231XAP"
ENCSR231ZZH="/cta/users/merve.kaftancioglu/alternative_scenario/processed_seqs/processed_ENCSR231ZZH"
ENCSR234BAD="/cta/users/merve.kaftancioglu/alternative_scenario/processed_seqs/processed_ENCSR234BAD"
ENCSR238WFK="/cta/users/merve.kaftancioglu/alternative_scenario/processed_seqs/processed_ENCSR238WFK"

ENCSR239OWL="/cta/users/merve.kaftancioglu/alternative_scenario/processed_seqs/processed_ENCSR239OWL"
ENCSR245KON="/cta/users/merve.kaftancioglu/alternative_scenario/processed_seqs/processed_ENCSR245KON"
ENCSR246TTM="/cta/users/merve.kaftancioglu/alternative_scenario/processed_seqs/processed_ENCSR246TTM"
ENCSR260CRI="/cta/users/merve.kaftancioglu/alternative_scenario/processed_seqs/processed_ENCSR260CRI"
ENCSR272YVX="/cta/users/merve.kaftancioglu/alternative_scenario/processed_seqs/processed_ENCSR272YVX"

ENCSR276NGH="/cta/users/merve.kaftancioglu/alternative_scenario/processed_seqs/processed_ENCSR276NGH"
ENCSR277XYX="/cta/users/merve.kaftancioglu/alternative_scenario/processed_seqs/processed_ENCSR277XYX"
ENCSR277YKG="/cta/users/merve.kaftancioglu/alternative_scenario/processed_seqs/processed_ENCSR277YKG"
ENCSR278QHR="/cta/users/merve.kaftancioglu/alternative_scenario/processed_seqs/processed_ENCSR278QHR"
ENCSR284IKS="/cta/users/merve.kaftancioglu/alternative_scenario/processed_seqs/processed_ENCSR284IKS"

ENCSR294HTM="/cta/users/merve.kaftancioglu/alternative_scenario/processed_seqs/processed_ENCSR294HTM"
ENCSR295QZX="/cta/users/merve.kaftancioglu/alternative_scenario/processed_seqs/processed_ENCSR295QZX"
ENCSR303SQG="/cta/users/merve.kaftancioglu/alternative_scenario/processed_seqs/processed_ENCSR303SQG"
ENCSR315MYP="/cta/users/merve.kaftancioglu/alternative_scenario/processed_seqs/processed_ENCSR315MYP"
ENCSR322MTA="/cta/users/merve.kaftancioglu/alternative_scenario/processed_seqs/processed_ENCSR322MTA"

ENCSR330CQU="/cta/users/merve.kaftancioglu/alternative_scenario/processed_seqs/processed_ENCSR330CQU"
ENCSR331WMS="/cta/users/merve.kaftancioglu/alternative_scenario/processed_seqs/processed_ENCSR331WMS"
ENCSR341QLC="/cta/users/merve.kaftancioglu/alternative_scenario/processed_seqs/processed_ENCSR341QLC"
ENCSR348YRH="/cta/users/merve.kaftancioglu/alternative_scenario/processed_seqs/processed_ENCSR348YRH"
ENCSR350UKV="/cta/users/merve.kaftancioglu/alternative_scenario/processed_seqs/processed_ENCSR350UKV"

ENCSR354GNT="/cta/users/merve.kaftancioglu/alternative_scenario/processed_seqs/processed_ENCSR354GNT"
ENCSR364OIK="/cta/users/merve.kaftancioglu/alternative_scenario/processed_seqs/processed_ENCSR364OIK"
ENCSR377QYB="/cta/users/merve.kaftancioglu/alternative_scenario/processed_seqs/processed_ENCSR377QYB"
ENCSR385BOZ="/cta/users/merve.kaftancioglu/alternative_scenario/processed_seqs/processed_ENCSR385BOZ"
ENCSR394JFQ="/cta/users/merve.kaftancioglu/alternative_scenario/processed_seqs/processed_ENCSR394JFQ"

ENCSR433EWI="/cta/users/merve.kaftancioglu/alternative_scenario/processed_seqs/processed_ENCSR433EWI"
ENCSR435MSB="/cta/users/merve.kaftancioglu/alternative_scenario/processed_seqs/processed_ENCSR435MSB"
ENCSR445DNH="/cta/users/merve.kaftancioglu/alternative_scenario/processed_seqs/processed_ENCSR445DNH"
ENCSR445LTM="/cta/users/merve.kaftancioglu/alternative_scenario/processed_seqs/processed_ENCSR445LTM"
ENCSR458QEY="/cta/users/merve.kaftancioglu/alternative_scenario/processed_seqs/processed_ENCSR458QEY"

ENCSR458TOW="/cta/users/merve.kaftancioglu/alternative_scenario/processed_seqs/processed_ENCSR458TOW"
ENCSR461QMZ="/cta/users/merve.kaftancioglu/alternative_scenario/processed_seqs/processed_ENCSR461QMZ"
ENCSR469CFL="/cta/users/merve.kaftancioglu/alternative_scenario/processed_seqs/processed_ENCSR469CFL"
ENCSR469QVG="/cta/users/merve.kaftancioglu/alternative_scenario/processed_seqs/processed_ENCSR469QVG"
ENCSR471WZE="/cta/users/merve.kaftancioglu/alternative_scenario/processed_seqs/processed_ENCSR471WZE"

ENCSR473TGK="/cta/users/merve.kaftancioglu/alternative_scenario/processed_seqs/processed_ENCSR473TGK"
ENCSR474PYR="/cta/users/merve.kaftancioglu/alternative_scenario/processed_seqs/processed_ENCSR474PYR"
ENCSR476IPR="/cta/users/merve.kaftancioglu/alternative_scenario/processed_seqs/processed_ENCSR476IPR"
ENCSR480ITK="/cta/users/merve.kaftancioglu/alternative_scenario/processed_seqs/processed_ENCSR480ITK"
ENCSR482TGI="/cta/users/merve.kaftancioglu/alternative_scenario/processed_seqs/processed_ENCSR482TGI"

ENCSR482UGV="/cta/users/merve.kaftancioglu/alternative_scenario/processed_seqs/processed_ENCSR482UGV"
ENCSR482VIB="/cta/users/merve.kaftancioglu/alternative_scenario/processed_seqs/processed_ENCSR482VIB"
ENCSR485FBT="/cta/users/merve.kaftancioglu/alternative_scenario/processed_seqs/processed_ENCSR485FBT"
ENCSR486SZN="/cta/users/merve.kaftancioglu/alternative_scenario/processed_seqs/processed_ENCSR486SZN"
ENCSR486XJK="/cta/users/merve.kaftancioglu/alternative_scenario/processed_seqs/processed_ENCSR486XJK"

ENCSR496LKR="/cta/users/merve.kaftancioglu/alternative_scenario/processed_seqs/processed_ENCSR496LKR"
ENCSR516CKJ="/cta/users/merve.kaftancioglu/alternative_scenario/processed_seqs/processed_ENCSR516CKJ"
ENCSR520QDR="/cta/users/merve.kaftancioglu/alternative_scenario/processed_seqs/processed_ENCSR520QDR"
ENCSR521SFR="/cta/users/merve.kaftancioglu/alternative_scenario/processed_seqs/processed_ENCSR521SFR"
ENCSR526TNC="/cta/users/merve.kaftancioglu/alternative_scenario/processed_seqs/processed_ENCSR526TNC"

ENCSR530AVK="/cta/users/merve.kaftancioglu/alternative_scenario/processed_seqs/processed_ENCSR530AVK"
ENCSR530YDY="/cta/users/merve.kaftancioglu/alternative_scenario/processed_seqs/processed_ENCSR530YDY"
ENCSR532PXR="/cta/users/merve.kaftancioglu/alternative_scenario/processed_seqs/processed_ENCSR532PXR"
ENCSR535YYH="/cta/users/merve.kaftancioglu/alternative_scenario/processed_seqs/processed_ENCSR535YYH"
ENCSR537KJA="/cta/users/merve.kaftancioglu/alternative_scenario/processed_seqs/processed_ENCSR537KJA"

ENCSR538URI="/cta/users/merve.kaftancioglu/alternative_scenario/processed_seqs/processed_ENCSR538URI"
ENCSR540XNK="/cta/users/merve.kaftancioglu/alternative_scenario/processed_seqs/processed_ENCSR540XNK"
ENCSR541PLO="/cta/users/merve.kaftancioglu/alternative_scenario/processed_seqs/processed_ENCSR541PLO"
ENCSR550DPT="/cta/users/merve.kaftancioglu/alternative_scenario/processed_seqs/processed_ENCSR550DPT"
ENCSR565WDW="/cta/users/merve.kaftancioglu/alternative_scenario/processed_seqs/processed_ENCSR565WDW"

ENCSR572XTB="/cta/users/merve.kaftancioglu/alternative_scenario/processed_seqs/processed_ENCSR572XTB"
ENCSR589LHR="/cta/users/merve.kaftancioglu/alternative_scenario/processed_seqs/processed_ENCSR589LHR"
ENCSR592EKF="/cta/users/merve.kaftancioglu/alternative_scenario/processed_seqs/processed_ENCSR592EKF"
ENCSR603LTN="/cta/users/merve.kaftancioglu/alternative_scenario/processed_seqs/processed_ENCSR603LTN"
ENCSR611YIA="/cta/users/merve.kaftancioglu/alternative_scenario/processed_seqs/processed_ENCSR611YIA"

ENCSR613UFD="/cta/users/merve.kaftancioglu/alternative_scenario/processed_seqs/processed_ENCSR613UFD"
ENCSR617GMR="/cta/users/merve.kaftancioglu/alternative_scenario/processed_seqs/processed_ENCSR617GMR"
ENCSR631VIW="/cta/users/merve.kaftancioglu/alternative_scenario/processed_seqs/processed_ENCSR631VIW"
ENCSR639NIG="/cta/users/merve.kaftancioglu/alternative_scenario/processed_seqs/processed_ENCSR639NIG"
ENCSR641ZFV="/cta/users/merve.kaftancioglu/alternative_scenario/processed_seqs/processed_ENCSR641ZFV"

ENCSR649FUX="/cta/users/merve.kaftancioglu/alternative_scenario/processed_seqs/processed_ENCSR649FUX"
ENCSR652WFO="/cta/users/merve.kaftancioglu/alternative_scenario/processed_seqs/processed_ENCSR652WFO"
ENCSR667HUI="/cta/users/merve.kaftancioglu/alternative_scenario/processed_seqs/processed_ENCSR667HUI"
ENCSR677OEF="/cta/users/merve.kaftancioglu/alternative_scenario/processed_seqs/processed_ENCSR677OEF"
ENCSR685KZA="/cta/users/merve.kaftancioglu/alternative_scenario/processed_seqs/processed_ENCSR685KZA"

ENCSR687CQX="/cta/users/merve.kaftancioglu/alternative_scenario/processed_seqs/processed_ENCSR687CQX"
ENCSR703TWY="/cta/users/merve.kaftancioglu/alternative_scenario/processed_seqs/processed_ENCSR703TWY"
ENCSR705VSO="/cta/users/merve.kaftancioglu/alternative_scenario/processed_seqs/processed_ENCSR705VSO"
ENCSR713ZYF="/cta/users/merve.kaftancioglu/alternative_scenario/processed_seqs/processed_ENCSR713ZYF"
ENCSR718LCI="/cta/users/merve.kaftancioglu/alternative_scenario/processed_seqs/processed_ENCSR718LCI"

ENCSR720QRX="/cta/users/merve.kaftancioglu/alternative_scenario/processed_seqs/processed_ENCSR720QRX"
ENCSR729ZXH="/cta/users/merve.kaftancioglu/alternative_scenario/processed_seqs/processed_ENCSR729ZXH"
ENCSR733LCG="/cta/users/merve.kaftancioglu/alternative_scenario/processed_seqs/processed_ENCSR733LCG"
ENCSR740SDR="/cta/users/merve.kaftancioglu/alternative_scenario/processed_seqs/processed_ENCSR740SDR"
ENCSR741XAE="/cta/users/merve.kaftancioglu/alternative_scenario/processed_seqs/processed_ENCSR741XAE"

ENCSR746AIX="/cta/users/merve.kaftancioglu/alternative_scenario/processed_seqs/processed_ENCSR746AIX"
ENCSR755ZGS="/cta/users/merve.kaftancioglu/alternative_scenario/processed_seqs/processed_ENCSR755ZGS"
ENCSR757FGN="/cta/users/merve.kaftancioglu/alternative_scenario/processed_seqs/processed_ENCSR757FGN"
ENCSR759DHN="/cta/users/merve.kaftancioglu/alternative_scenario/processed_seqs/processed_ENCSR759DHN"
ENCSR779JLY="/cta/users/merve.kaftancioglu/alternative_scenario/processed_seqs/processed_ENCSR779JLY"

ENCSR785PKM="/cta/users/merve.kaftancioglu/alternative_scenario/processed_seqs/processed_ENCSR785PKM"
ENCSR787HDF="/cta/users/merve.kaftancioglu/alternative_scenario/processed_seqs/processed_ENCSR787HDF"
ENCSR788TEF="/cta/users/merve.kaftancioglu/alternative_scenario/processed_seqs/processed_ENCSR788TEF"
ENCSR791CAF="/cta/users/merve.kaftancioglu/alternative_scenario/processed_seqs/processed_ENCSR791CAF"
ENCSR802MXQ="/cta/users/merve.kaftancioglu/alternative_scenario/processed_seqs/processed_ENCSR802MXQ"

ENCSR806HUT="/cta/users/merve.kaftancioglu/alternative_scenario/processed_seqs/processed_ENCSR806HUT"
ENCSR815YZL="/cta/users/merve.kaftancioglu/alternative_scenario/processed_seqs/processed_ENCSR815YZL"
ENCSR819NCZ="/cta/users/merve.kaftancioglu/alternative_scenario/processed_seqs/processed_ENCSR819NCZ"
ENCSR831AXK="/cta/users/merve.kaftancioglu/alternative_scenario/processed_seqs/processed_ENCSR831AXK"
ENCSR832UMM="/cta/users/merve.kaftancioglu/alternative_scenario/processed_seqs/processed_ENCSR832UMM"

ENCSR842WWX="/cta/users/merve.kaftancioglu/alternative_scenario/processed_seqs/processed_ENCSR842WWX"
ENCSR848XJL="/cta/users/merve.kaftancioglu/alternative_scenario/processed_seqs/processed_ENCSR848XJL"
ENCSR851RJV="/cta/users/merve.kaftancioglu/alternative_scenario/processed_seqs/processed_ENCSR851RJV"
ENCSR861FEC="/cta/users/merve.kaftancioglu/alternative_scenario/processed_seqs/processed_ENCSR861FEC"
ENCSR864OKB="/cta/users/merve.kaftancioglu/alternative_scenario/processed_seqs/processed_ENCSR864OKB"

ENCSR865FTW="/cta/users/merve.kaftancioglu/alternative_scenario/processed_seqs/processed_ENCSR865FTW"
ENCSR878YHM="/cta/users/merve.kaftancioglu/alternative_scenario/processed_seqs/processed_ENCSR878YHM"
ENCSR898VJE="/cta/users/merve.kaftancioglu/alternative_scenario/processed_seqs/processed_ENCSR898VJE"
ENCSR900IAM="/cta/users/merve.kaftancioglu/alternative_scenario/processed_seqs/processed_ENCSR900IAM"
ENCSR919DFZ="/cta/users/merve.kaftancioglu/alternative_scenario/processed_seqs/processed_ENCSR919DFZ"

ENCSR923EGO="/cta/users/merve.kaftancioglu/alternative_scenario/processed_seqs/processed_ENCSR923EGO"
ENCSR923JIB="/cta/users/merve.kaftancioglu/alternative_scenario/processed_seqs/processed_ENCSR923JIB"
ENCSR931JYC="/cta/users/merve.kaftancioglu/alternative_scenario/processed_seqs/processed_ENCSR931JYC"
ENCSR940PHE="/cta/users/merve.kaftancioglu/alternative_scenario/processed_seqs/processed_ENCSR940PHE"
ENCSR953STZ="/cta/users/merve.kaftancioglu/alternative_scenario/processed_seqs/processed_ENCSR953STZ"

ENCSR954ZLD="/cta/users/merve.kaftancioglu/alternative_scenario/processed_seqs/processed_ENCSR954ZLD"
ENCSR959VZU="/cta/users/merve.kaftancioglu/alternative_scenario/processed_seqs/processed_ENCSR959VZU"
ENCSR976RWL="/cta/users/merve.kaftancioglu/alternative_scenario/processed_seqs/processed_ENCSR976RWL"
ENCSR977FMZ="/cta/users/merve.kaftancioglu/alternative_scenario/processed_seqs/processed_ENCSR977FMZ"
ENCSR980KTW="/cta/users/merve.kaftancioglu/alternative_scenario/processed_seqs/processed_ENCSR980KTW"

ENCSR983MKX="/cta/users/merve.kaftancioglu/alternative_scenario/processed_seqs/processed_ENCSR983MKX"
ENCSR987OPY="/cta/users/merve.kaftancioglu/alternative_scenario/processed_seqs/processed_ENCSR987OPY"
ENCSR993CTA="/cta/users/merve.kaftancioglu/alternative_scenario/processed_seqs/processed_ENCSR993CTA"

""" !/bin/bash
SBATCH --job-name=trimmomatic_single
SBATCH --nodes=1
SBATCH --account=users
SBATCH --qos=mid_mdbf
SBATCH --partition=mid_mdbf
SBATCH --cpus-per-task=16
SBATCH -o /cta/users/merve.kaftancioglu/alternative_scenario/trimmomatic_as/out-err/%x-%j-slurm.out
SBATCH -e /cta/users/merve.kaftancioglu/alternative_scenario/trimmomatic_as/out-err/%x-%j-slurm.err

 Load FastQC module
module load trimmomatic-0.38-gcc-9.2.0-3evdnxc

 List of sample IDs (folders)
sample_ids=("h3k27ac_ms_b_ENCSR617GMR" "h3k27ac_ms_cd4_ENCSR832UMM" "h3k27ac_ms_cd8_ENCSR078ATS")

 Process each sample folder
for sample_id in "${sample_ids[@]}"; do
    input_folder="/cta/users/merve.kaftancioglu/alternative_scenario/seqs/$sample_id"
    output_folder="/cta/users/merve.kaftancioglu/alternative_scenario/trimmomatic_as/output/$sample_id"

     Create output folder if it doesn't exist
    mkdir -p "$output_folder"

     Process all .fastq.gz files in the specified folder
    for input_file in "$input_folder"/*.fastq.gz; do
        output_file="$output_folder/output_paired_$(basename "$input_file")"
        trimmomatic SE -threads 16 -phred33 -trimlog "$output_folder/trimmomatic_log.txt" "$input_file" "$output_file" HEADCROP:14
    done
done

""" !/bin/bash
SBATCH --job-name=trimmomatic_single
SBATCH --nodes=1
SBATCH --account=users
SBATCH --qos=mid_mdbf
SBATCH --partition=mid_mdbf
SBATCH --cpus-per-task=16
SBATCH -o /cta/users/merve.kaftancioglu/alternative_scenario/trimmomatic_as/out-err/%x-%j-slurm.out
SBATCH -e /cta/users/merve.kaftancioglu/alternative_scenario/trimmomatic_as/out-err/%x-%j-slurm.err

 Load FastQC module
module load trimmomatic-0.38-gcc-9.2.0-3evdnxc

 List of sample IDs (folders)
sample_ids=("h3k9me3_ms_cd4_ENCSR953STZ")

 Process each sample folder
for sample_id in "${sample_ids[@]}"; do
    input_folder="/cta/users/merve.kaftancioglu/alternative_scenario/seqs/$sample_id"
    output_folder="/cta/users/merve.kaftancioglu/alternative_scenario/trimmomatic_as/output/$sample_id"

     Create output folder if it doesn't exist
    mkdir -p "$output_folder"

     Process all .fastq.gz files in the specified folder
    for input_file in "$input_folder"/*.fastq.gz; do
        output_file="$output_folder/output_paired_$(basename "$input_file")"
        trimmomatic SE -threads 16 -phred33 -trimlog "$output_folder/trimmomatic_log.txt" "$input_file" "$output_file" HEADCROP:14
    done
done

!/usr/bin/env bash

 This script is used to build and deploy the GCR docker image for Picard
 The dockerhub image is built using a dockerhub automated build: https://hub.docker.com/r/broadinstitute/picard/builds

if [[ "$1" == "" ]]
then
    echo "Usage: build_push_docker.sh <git-tag>"
    exit 1
fi

declare -r TAG=${1}

declare -r PICARD_CLOUD_TAG=us.gcr.io/broad-gotc-prod/picard-cloud:${TAG}

echo "Will build and push the following docker images:"
echo ${PICARD_CLOUD_TAG}

read -p "Is this really what you want to do? " -n 1 -r
echo     (optional) move to a new line
if [[ $REPLY =~ ^[Yy]$ ]]
then
    docker build -t ${PICARD_CLOUD_TAG} --build-arg release=true .
    docker push ${PICARD_CLOUD_TAG}
fi

""" !/bin/bash
SBATCH --job-name=fastQC
SBATCH --nodes=1
SBATCH --account=users
SBATCH --qos=short_mdbf
SBATCH --partition=short_mdbf
SBATCH --cpus-per-task=8
SBATCH -o /cta/users/merve.kaftancioglu/_GWSTA-2024-Workshop-main/4_fastQC/%x-%j-slurm.out
SBATCH -e /cta/users/merve.kaftancioglu/_GWSTA-2024-Workshop-main/4_fastQC/%x-%j-slurm.err


input_directory="/cta/users/merve.kaftancioglu/_GWSTA-2024-Workshop-main/4_fastQC"
output_directory="/cta/users/merve.kaftancioglu/_GWSTA-2024-Workshop-main/4_fastQC"

 Load FastQC module
module load fastqc-0.11.7-gcc-9.2.0-nsjzjof

 Run FastQC on .fastq files in input directory and write output to output directory
fastqc --outdir "$output_directory" "$input_directory"/*.fastq -threads 8



""" !/bin/bash
SBATCH --job-name=fastQC
SBATCH --nodes=1
SBATCH --account=users
SBATCH --partition=short_mdbf
SBATCH --qos=short_mdbf 
SBATCH -o /cta/users/merve.kaftancioglu/_GWSTA-2024-Workshop-main/4_fastQC/%x-%j-slurm.out
SBATCH -e /cta/users/merve.kaftancioglu/_GWSTA-2024-Workshop-main/4_fastQC/%x-%j-slurm.err



 Define the directory to save the downloaded files
output_directory=/cta/users/merve.kaftancioglu/_GWSTA-2024-Workshop-main/4_fastQC/

 Create the output directory if it doesn't exist
mkdir -p "$output_directory"

 Change to the output directory
cd "$output_directory" || exit

 Download the files
wget -nc ftp://ftp.sra.ebi.ac.uk/vol1/fastq/SRR001/SRR001666/SRR001666_2.fastq.gz
wget -nc ftp://ftp.sra.ebi.ac.uk/vol1/fastq/SRR001/SRR001666/SRR001666_1.fastq.gz

""" !/bin/bash
SBATCH --job-name=bwa_align
SBATCH --nodes=1
SBATCH --account=users
SBATCH --qos=short_mdbf
SBATCH --partition=short_mdbf
SBATCH --cpus-per-task=8
SBATCH -o %x-%j-slurm.out
SBATCH -e %x-%j-slurm.err

 Load BWA module

module load bwa-0.7.17-gcc-9.2.0-xwkclqy


 Run bwa mem for alignment
bwa mem -t 8 index/Agy99.fasta P7741_R1.fastq.gz P7741_R2.fastq.gz > output.sam

 -t: threads
""" !/bin/bash
SBATCH --job-name=bwa
SBATCH --nodes=1
SBATCH --account=users
SBATCH --qos=short_mdbf
SBATCH --partition=short_mdbf
SBATCH --cpus-per-task=8
SBATCH -o %x-%j-slurm.out
SBATCH -e %x-%j-slurm.err

 Load BWA module

module load bwa-0.7.17-gcc-9.2.0-xwkclqy


 Run bwa on .fastq files and reference fasta file
bwa mem -t 8 index/sequence.fasta SRR4051933_1.fastq.gz SRR4051933_2.fastq.gz > output_co.sam

bwa mem -t 8 index/sequence.fasta SRR24878810_1.fastq.gz SRR24878810_2.fastq.gz > output_pre.sam

 -t: threads

!/usr/bin/env bash
 Copyright 2009 The Go Authors. All rights reserved.
 Use of this source code is governed by a BSD-style
 license that can be found in the LICENSE file.

set -e
if [ ! -f make.bash ]; then
	echo 'all.bash must be run from $GOROOT/src' 1>&2
	exit 1
fi
. ./make.bash "$@" --no-banner
bash run.bash --no-rebuild
$GOTOOLDIR/dist banner   print build info

""" !/bin/bash
SBATCH --job-name=trim
SBATCH --nodes=1
SBATCH --account=users
SBATCH --qos=short_mdbf
SBATCH --partition=short_mdbf
SBATCH --cpus-per-task=4
SBATCH -o %x-%j-slurm.out
SBATCH -e %x-%j-slurm.err

output_directory="reports"

 Load trimmomatic module
module load trimmomatic-0.38-gcc-9.2.0-3evdnxc

 Run trimmomatic
trimmomatic PE -threads 4 -phred33 -trimlog ./trimmomatic_log.txt SRR1553606_1.fastq.gz SRR1553606_2.fastq.gz  output_SRR1553606_forward_paired.fq.gz output_SRR1553606_forward_unpaired.fq.gz output_SRR1553606_reverse_paired.fq.gz output_SRR1553606_reverse_unpaired.fq.gz ILLUMINACLIP:TruSeq3-PE-2.fa:2:30:10 LEADING:3 TRAILING:3 SLIDINGWINDOW:4:20 MINLEN:30
