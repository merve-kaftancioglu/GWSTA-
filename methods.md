## Bash Scripts:
# FastQC (samples divided into 31 files): 

```
#!/bin/bash
#SBATCH --job-name=fastQC
#SBATCH --nodes=1
#SBATCH --account=users
#SBATCH --qos=mid_mdbf
#SBATCH --partition=mid_mdbf
#SBATCH --cpus-per-task=16
#SBATCH -o /cta/users/merve.kaftancioglu/alternative_scenario/post-trim_fastQC_as/out-err/%x-%j-slurm.out
#SBATCH -e /cta/users/merve.kaftancioglu/alternative_scenario/post-trim_fastQC_as/out-err/%x-%j-slurm.err

# Load FastQC module
module load fastqc-0.11.7-gcc-9.2.0-nsjzjof

sample_ids=("h3k27ac_ms_b_ENCSR295QZX" "h3k27ac_ms_b_ENCSR538URI" "h3k27ac_ms_b_ENCSR617GMR" "h3k4me1_ms_b_ENCSR084FXT" "h3k4me1_ms_b_ENCSR480ITK" "h3k4me1_ms_b_ENCSR759DHN" "h3k4me1_ms_b_ENCSR900IAM" "h3k4me1_ms_b_ENCSR931JYC" "h3k4me3_ms_b_ENCSR260CRI" "h3k4me3_ms_b_ENCSR461QMZ" "h3k4me3_ms_b_ENCSR530AVK" "h3k4me3_ms_b_ENCSR713ZYF" "h3k4me3_ms_b_ENCSR923EGO" "h3k9me3_ms_b_ENCSR445DNH" "h3k9me3_ms_b_ENCSR486SZN" "h3k9me3_ms_b_ENCSR983MKX" "h3k27ac_ms_b_ENCSR364OIK" "h3k27ac_ms_b_ENCSR541PLO" "h3k36me3_ms_cd4_ENCSR471WZE" "h3k36me3_ms_cd4_ENCSR473TGK" "h3k4me1_ms_nk_ENCSR482UGV" "h3k4me1_ms_nk_ENCSR718LCI" "h3k4me1_ms_nk_ENCSR806HUT" "h3k4me3_ms_nk_ENCSR234BAD" "h3k4me3_ms_nk_ENCSR703TWY" "h3k4me3_ms_nk_ENCSR755ZGS" "h3k9me3_ms_nk_ENCSR061ATV" "h3k9me3_ms_nk_ENCSR611YIA" "h3k9me3_ms_nk_ENCSR667HUI" "h3k27ac_ms_nk_ENCSR246TTM" "h3k27ac_ms_nk_ENCSR469CFL" "h3k27ac_ms_nk_ENCSR746AIX" "h3k27me3_ms_b_ENCSR009ZRH" "h3k27me3_ms_b_ENCSR182NLA" "h3k27me3_ms_b_ENCSR272YVX" "h3k27me3_ms_b_ENCSR649FUX" "h3k27me3_ms_b_ENCSR842WWX" "h3k36me3_ms_b_ENCSR089VQL" "h3k36me3_ms_b_ENCSR119PSR" "h3k36me3_ms_b_ENCSR238WFK" "h3k36me3_ms_b_ENCSR987OPY" "h3k4me1_ms_cd4_ENCSR036YVP" "h3k4me1_ms_cd4_ENCSR043JIA" "h3k4me1_ms_cd4_ENCSR093VUP" "h3k4me1_ms_cd4_ENCSR124IGC" "h3k4me1_ms_cd4_ENCSR315MYP" "h3k4me1_ms_cd4_ENCSR458QEY" "h3k4me1_ms_cd4_ENCSR485FBT" "h3k4me1_ms_cd4_ENCSR641ZFV" "h3k4me1_ms_cd4_ENCSR687CQX" "h3k4me1_ms_cd8_ENCSR231XAP" "h3k4me1_ms_cd8_ENCSR572XTB" "h3k4me1_ms_cd8_ENCSR788TEF" "h3k4me1_ms_cd8_ENCSR815YZL" "h3k4me1_ms_cd8_ENCSR861FEC" "h3k4me1_ms_cd8_ENCSR940PHE" "h3k4me3_ms_cd4_ENCSR180NCM" "h3k4me3_ms_cd4_ENCSR341QLC" "h3k4me3_ms_cd4_ENCSR482TGI" "h3k4me3_ms_cd4_ENCSR486XJK" "h3k4me3_ms_cd4_ENCSR496LKR" "h3k4me3_ms_cd4_ENCSR603LTN" "h3k4me3_ms_cd4_ENCSR802MXQ" "h3k4me3_ms_cd4_ENCSR878YHM" "h3k4me3_ms_cd4_ENCSR954ZLD" "h3k4me3_ms_cd8_ENCSR231ZZH" "h3k4me3_ms_cd8_ENCSR278QHR" "h3k4me3_ms_cd8_ENCSR516CKJ" "h3k4me3_ms_cd8_ENCSR535YYH" "h3k4me3_ms_cd8_ENCSR741XAE" "h3k4me3_ms_cd8_ENCSR848XJL" "h3k9me3_ms_cd4_ENCSR057IZD" "h3k9me3_ms_cd4_ENCSR550DPT" "h3k9me3_ms_cd4_ENCSR677OEF" "h3k9me3_ms_cd4_ENCSR729ZXH" "h3k9me3_ms_cd4_ENCSR851RJV" "h3k9me3_ms_cd4_ENCSR919DFZ" "h3k9me3_ms_cd4_ENCSR953STZ" "h3k9me3_ms_cd4_ENCSR959VZU" "h3k9me3_ms_cd8_ENCSR101USF" "h3k9me3_ms_cd8_ENCSR354GNT" "h3k9me3_ms_cd8_ENCSR377QYB" "h3k9me3_ms_cd8_ENCSR733LCG" "h3k9me3_ms_cd8_ENCSR980KTW" "h3k27ac_ms_cd4_ENCSR200SSJ" "h3k27ac_ms_cd4_ENCSR322MTA" "h3k27ac_ms_cd4_ENCSR331WMS" "h3k27ac_ms_cd4_ENCSR350UKV" "h3k27ac_ms_cd4_ENCSR474PYR" "h3k27ac_ms_cd4_ENCSR520QDR" "h3k27ac_ms_cd4_ENCSR540XNK" "h3k27ac_ms_cd4_ENCSR705VSO" "h3k27ac_ms_cd4_ENCSR832UMM" "h3k27ac_ms_cd8_ENCSR078ATS" "h3k27ac_ms_cd8_ENCSR348YRH" "h3k27ac_ms_cd8_ENCSR458TOW" "h3k27ac_ms_cd8_ENCSR476IPR" "h3k27ac_ms_cd8_ENCSR787HDF" "h3k27ac_ms_cd8_ENCSR923JIB" "h3k27me3_ms_nk_ENCSR469QVG" "h3k27me3_ms_nk_ENCSR565WDW" "h3k36me3_ms_nk_ENCSR158VSE" "h3k36me3_ms_nk_ENCSR245KON" "h3k36me3_ms_nk_ENCSR530YDY" "h3k36me_ms_cd4_ENCSR532PXR" "h3k27me3_ms_cd4_ENCSR277XYX" "h3k27me3_ms_cd4_ENCSR526TNC" "h3k27me3_ms_cd4_ENCSR592EKF" "h3k27me3_ms_cd4_ENCSR613UFD" "h3k27me3_ms_cd4_ENCSR740SDR" "h3k27me3_ms_cd4_ENCSR779JLY" "h3k27me3_ms_cd4_ENCSR993CTA" "h3k27me3_ms_cd8_ENCSR116FVG" "h3k27me3_ms_cd8_ENCSR122JCM" "h3k27me3_ms_cd8_ENCSR216ZVA" "h3k27me3_ms_cd8_ENCSR284IKS" "h3k27me3_ms_cd8_ENCSR385BOZ" "h3k27me3_ms_cd8_ENCSR521SFR" "h3k36me3_ms_cd4_ENCSR276NGH" "h3k36me3_ms_cd4_ENCSR330CQU" "h3k36me3_ms_cd4_ENCSR482VIB" "h3k36me3_ms_cd4_ENCSR785PKM" "h3k36me3_ms_cd4_ENCSR865FTW" "h3k36me3_ms_cd4_ENCSR898VJE" "h3k36me3_ms_cd8_ENCSR239OWL" "h3k36me3_ms_cd8_ENCSR435MSB" "h3k36me3_ms_cd8_ENCSR631VIW" "h3k36me3_ms_cd8_ENCSR652WFO" "h3k36me3_ms_cd8_ENCSR757FGN" "h3k4me1_normal_b_ENCSR156BXM" "h3k4me3_normal_b_ENCSR791CAF" "h3k9me3_normal_b_ENCSR445LTM" "h3k27ac_normal_b_ENCSR685KZA" "h3k4me1_normal_nk_ENCSR277YKG" "h3k4me3_normal_nk_ENCSR394JFQ" "h3k9me3_normal_nk_ENCSR025UNZ" "h3k27ac_normal_nk_ENCSR977FMZ" "h3k27me3_normal_b_ENCSR589LHR" "h3k36me3_normal_b_ENCSR831AXK" "h3k4me1_normal_cd4_ENCSR102SOR" "h3k4me1_normal_cd8_ENCSR217SHH" "h3k4me3_normal_cd4_ENCSR537KJA" "h3k4me3_normal_cd8_ENCSR123ZAT" "h3k9me3_normal_cd4_ENCSR433EWI" "h3k9me3_normal_cd8_ENCSR294HTM" "h3k27ac_normal_cd4_ENCSR819NCZ" "h3k27ac_normal_cd8_ENCSR976RWL" "h3k27me3_normal_nk_ENCSR639NIG" "h3k36me3_normal_nk_ENCSR056GJY" "h3k27me3_normal_cd4_ENCSR068YVZ" "h3k27me3_normal_cd8_ENCSR720QRX" "h3k36me3_normal_cd4_ENCSR864OKB" "h3k36me3_normal_cd8_ENCSR303SQG" "h3k27me3_normal_cd4_ENCSR002UMT" "h3k4me3_normal_cd4_ENCSR935ELX")

for sample_id in "${sample_ids[@]}"; do
  echo "Predicting ${sample_id}..."
  input_directory="/cta/users/merve.kaftancioglu/alternative_scenario/trimmomatic_as/output/${sample_id}"
  output_dir="/cta/users/merve.kaftancioglu/alternative_scenario/post-trim_fastQC_as/output/${sample_id}"
  mkdir -p "${output_dir}"
  fastqc --outdir "$output_dir" "$input_directory"/*.fastq.gz -threads 16
  echo "Finished ${sample_id}"
done

echo "DONE!"
```
<br>

# MultiQC (without bash):

```
multiqc -o /cta/users/merve.kaftancioglu/alternative_scenario/post-trim_fastQC_as/output_merged_all /cta/users/merve.kaftancioglu/alternative_scenario/post-trim_fastQC_as/output_merged_all/*_fastqc.zip
multiqc -o /cta/users/merve.kaftancioglu/alternative_scenario/multiqc_data/h3k27ac_ms_b_ENCSR295QZX /cta/users/merve.kaftancioglu/alternative_scenario/fastQC_as/output/h3k27ac_ms_b_ENCSR295QZX/*_fastqc.zip	 
multiqc -o /cta/users/merve.kaftancioglu/alternative_scenario/multiqc_data/h3k4me1_ms_b_ENCSR480ITK /cta/users/merve.kaftancioglu/alternative_scenario/fastQC_as/output/h3k4me1_ms_b_ENCSR480ITK/*_fastqc.zip
multiqc -o /cta/users/merve.kaftancioglu/alternative_scenario/multiqc_data/h3k27ac_ms_b_ENCSR364OIK /cta/users/merve.kaftancioglu/alternative_scenario/fastQC_as/output/h3k27ac_ms_b_ENCSR364OIK/*_fastqc.zip
multiqc -o /cta/users/merve.kaftancioglu/alternative_scenario/multiqc_data/h3k4me1_ms_b_ENCSR759DHN /cta/users/merve.kaftancioglu/alternative_scenario/fastQC_as/output/h3k4me1_ms_b_ENCSR759DHN/*_fastqc.zip
multiqc -o /cta/users/merve.kaftancioglu/alternative_scenario/multiqc_data/h3k27ac_ms_b_ENCSR538URI /cta/users/merve.kaftancioglu/alternative_scenario/fastQC_as/output/h3k27ac_ms_b_ENCSR538URI/*_fastqc.zip	 
multiqc -o /cta/users/merve.kaftancioglu/alternative_scenario/multiqc_data/h3k4me1_ms_b_ENCSR900IAM /cta/users/merve.kaftancioglu/alternative_scenario/fastQC_as/output/h3k4me1_ms_b_ENCSR900IAM/*_fastqc.zip
multiqc -o /cta/users/merve.kaftancioglu/alternative_scenario/multiqc_data/h3k27ac_ms_b_ENCSR541PLO /cta/users/merve.kaftancioglu/alternative_scenario/fastQC_as/output/h3k27ac_ms_b_ENCSR541PLO/*_fastqc.zip	 
multiqc -o /cta/users/merve.kaftancioglu/alternative_scenario/multiqc_data/h3k4me1_ms_b_ENCSR931JYC /cta/users/merve.kaftancioglu/alternative_scenario/fastQC_as/output/h3k4me1_ms_b_ENCSR931JYC/*_fastqc.zip
multiqc -o /cta/users/merve.kaftancioglu/alternative_scenario/multiqc_data/h3k27ac_ms_b_ENCSR617GMR /cta/users/merve.kaftancioglu/alternative_scenario/fastQC_as/output/h3k27ac_ms_b_ENCSR617GMR/*_fastqc.zip	 
multiqc -o /cta/users/merve.kaftancioglu/alternative_scenario/multiqc_data/h3k4me1_ms_cd4_ENCSR036YVP /cta/users/merve.kaftancioglu/alternative_scenario/fastQC_as/output/h3k4me1_ms_cd4_ENCSR036YVP/*_fastqc.zip
multiqc -o /cta/users/merve.kaftancioglu/alternative_scenario/multiqc_data/h3k27ac_ms_cd4_ENCSR200SSJ /cta/users/merve.kaftancioglu/alternative_scenario/fastQC_as/output/h3k27ac_ms_cd4_ENCSR200SSJ/*_fastqc.zip
multiqc -o /cta/users/merve.kaftancioglu/alternative_scenario/multiqc_data/h3k4me1_ms_cd4_ENCSR043JIA /cta/users/merve.kaftancioglu/alternative_scenario/fastQC_as/output/h3k4me1_ms_cd4_ENCSR043JIA/*_fastqc.zip
multiqc -o /cta/users/merve.kaftancioglu/alternative_scenario/multiqc_data/h3k27ac_ms_cd4_ENCSR322MTA /cta/users/merve.kaftancioglu/alternative_scenario/fastQC_as/output/h3k27ac_ms_cd4_ENCSR322MTA/*_fastqc.zip	 
multiqc -o /cta/users/merve.kaftancioglu/alternative_scenario/multiqc_data/h3k4me1_ms_cd4_ENCSR093VUP /cta/users/merve.kaftancioglu/alternative_scenario/fastQC_as/output/h3k4me1_ms_cd4_ENCSR093VUP/*_fastqc.zip
multiqc -o /cta/users/merve.kaftancioglu/alternative_scenario/multiqc_data/h3k27ac_ms_cd4_ENCSR331WMS /cta/users/merve.kaftancioglu/alternative_scenario/fastQC_as/output/h3k27ac_ms_cd4_ENCSR331WMS/*_fastqc.zip	 
multiqc -o /cta/users/merve.kaftancioglu/alternative_scenario/multiqc_data/h3k4me1_ms_cd4_ENCSR124IGC /cta/users/merve.kaftancioglu/alternative_scenario/fastQC_as/output/h3k4me1_ms_cd4_ENCSR124IGC/*_fastqc.zip
multiqc -o /cta/users/merve.kaftancioglu/alternative_scenario/multiqc_data/h3k27ac_ms_cd4_ENCSR350UKV /cta/users/merve.kaftancioglu/alternative_scenario/fastQC_as/output/h3k27ac_ms_cd4_ENCSR350UKV/*_fastqc.zip	 
multiqc -o /cta/users/merve.kaftancioglu/alternative_scenario/multiqc_data/h3k4me1_ms_cd4_ENCSR315MYP /cta/users/merve.kaftancioglu/alternative_scenario/fastQC_as/output/h3k4me1_ms_cd4_ENCSR315MYP/*_fastqc.zip
multiqc -o /cta/users/merve.kaftancioglu/alternative_scenario/multiqc_data/h3k27ac_ms_cd4_ENCSR474PYR /cta/users/merve.kaftancioglu/alternative_scenario/fastQC_as/output/h3k27ac_ms_cd4_ENCSR474PYR/*_fastqc.zip	 
multiqc -o /cta/users/merve.kaftancioglu/alternative_scenario/multiqc_data/h3k4me1_ms_cd4_ENCSR458QEY /cta/users/merve.kaftancioglu/alternative_scenario/fastQC_as/output/h3k4me1_ms_cd4_ENCSR458QEY/*_fastqc.zip
multiqc -o /cta/users/merve.kaftancioglu/alternative_scenario/multiqc_data/h3k27ac_ms_cd4_ENCSR520QDR /cta/users/merve.kaftancioglu/alternative_scenario/fastQC_as/output/h3k27ac_ms_cd4_ENCSR520QDR/*_fastqc.zip	 
multiqc -o /cta/users/merve.kaftancioglu/alternative_scenario/multiqc_data/h3k4me1_ms_cd4_ENCSR485FBT /cta/users/merve.kaftancioglu/alternative_scenario/fastQC_as/output/h3k4me1_ms_cd4_ENCSR485FBT/*_fastqc.zip
multiqc -o /cta/users/merve.kaftancioglu/alternative_scenario/multiqc_data/h3k27ac_ms_cd4_ENCSR540XNK /cta/users/merve.kaftancioglu/alternative_scenario/fastQC_as/output/h3k27ac_ms_cd4_ENCSR540XNK/*_fastqc.zip	 
multiqc -o /cta/users/merve.kaftancioglu/alternative_scenario/multiqc_data/h3k4me1_ms_cd4_ENCSR641ZFV/ /cta/users/merve.kaftancioglu/alternative_scenario/fastQC_as/output/h3k4me1_ms_cd4_ENCSR641ZFV/*_fastqc.zip
multiqc -o /cta/users/merve.kaftancioglu/alternative_scenario/multiqc_data/h3k27ac_ms_cd4_ENCSR705VSO /cta/users/merve.kaftancioglu/alternative_scenario/fastQC_as/output/h3k27ac_ms_cd4_ENCSR705VSO/*_fastqc.zip	 
multiqc -o /cta/users/merve.kaftancioglu/alternative_scenario/multiqc_data/h3k4me1_ms_cd4_ENCSR687CQX /cta/users/merve.kaftancioglu/alternative_scenario/fastQC_as/output/h3k4me1_ms_cd4_ENCSR687CQX/*_fastqc.zip
multiqc -o /cta/users/merve.kaftancioglu/alternative_scenario/multiqc_data/h3k27ac_ms_cd4_ENCSR832UMM /cta/users/merve.kaftancioglu/alternative_scenario/fastQC_as/output/h3k27ac_ms_cd4_ENCSR832UMM/*_fastqc.zip	 
multiqc -o /cta/users/merve.kaftancioglu/alternative_scenario/multiqc_data/h3k4me1_ms_cd8_ENCSR231XAP /cta/users/merve.kaftancioglu/alternative_scenario/fastQC_as/output/h3k4me1_ms_cd8_ENCSR231XAP/*_fastqc.zip
multiqc -o /cta/users/merve.kaftancioglu/alternative_scenario/multiqc_data/h3k27ac_ms_cd8_ENCSR078ATS /cta/users/merve.kaftancioglu/alternative_scenario/fastQC_as/output/h3k27ac_ms_cd8_ENCSR078ATS/*_fastqc.zip	 
multiqc -o /cta/users/merve.kaftancioglu/alternative_scenario/multiqc_data/h3k4me1_ms_cd8_ENCSR572XTB /cta/users/merve.kaftancioglu/alternative_scenario/fastQC_as/output/h3k4me1_ms_cd8_ENCSR572XTB/*_fastqc.zip
multiqc -o /cta/users/merve.kaftancioglu/alternative_scenario/multiqc_data/h3k27ac_ms_cd8_ENCSR348YRH /cta/users/merve.kaftancioglu/alternative_scenario/fastQC_as/output/h3k27ac_ms_cd8_ENCSR348YRH/*_fastqc.zip	 
multiqc -o /cta/users/merve.kaftancioglu/alternative_scenario/multiqc_data/h3k4me1_ms_cd8_ENCSR788TEF /cta/users/merve.kaftancioglu/alternative_scenario/fastQC_as/output/h3k4me1_ms_cd8_ENCSR788TEF/*_fastqc.zip
multiqc -o /cta/users/merve.kaftancioglu/alternative_scenario/multiqc_data/h3k27ac_ms_cd8_ENCSR458TOW /cta/users/merve.kaftancioglu/alternative_scenario/fastQC_as/output/h3k27ac_ms_cd8_ENCSR458TOW/*_fastqc.zip	 
multiqc -o /cta/users/merve.kaftancioglu/alternative_scenario/multiqc_data/h3k4me1_ms_cd8_ENCSR815YZL /cta/users/merve.kaftancioglu/alternative_scenario/fastQC_as/output/h3k4me1_ms_cd8_ENCSR815YZL/*_fastqc.zip
multiqc -o /cta/users/merve.kaftancioglu/alternative_scenario/multiqc_data/h3k27ac_ms_cd8_ENCSR476IPR /cta/users/merve.kaftancioglu/alternative_scenario/fastQC_as/output/h3k27ac_ms_cd8_ENCSR476IPR/*_fastqc.zip	 
multiqc -o /cta/users/merve.kaftancioglu/alternative_scenario/multiqc_data/h3k4me1_ms_cd8_ENCSR861FEC /cta/users/merve.kaftancioglu/alternative_scenario/fastQC_as/output/h3k4me1_ms_cd8_ENCSR861FEC/*_fastqc.zip
multiqc -o /cta/users/merve.kaftancioglu/alternative_scenario/multiqc_data/h3k27ac_ms_cd8_ENCSR787HDF /cta/users/merve.kaftancioglu/alternative_scenario/fastQC_as/output/h3k27ac_ms_cd8_ENCSR787HDF/*_fastqc.zip	 
multiqc -o /cta/users/merve.kaftancioglu/alternative_scenario/multiqc_data/h3k4me1_ms_cd8_ENCSR940PHE /cta/users/merve.kaftancioglu/alternative_scenario/fastQC_as/output/h3k4me1_ms_cd8_ENCSR940PHE/*_fastqc.zip
multiqc -o /cta/users/merve.kaftancioglu/alternative_scenario/multiqc_data/h3k27ac_ms_cd8_ENCSR923JIB /cta/users/merve.kaftancioglu/alternative_scenario/fastQC_as/output/h3k27ac_ms_cd8_ENCSR923JIB/*_fastqc.zip	 
multiqc -o /cta/users/merve.kaftancioglu/alternative_scenario/multiqc_data/h3k4me1_ms_nk_ENCSR482UGV /cta/users/merve.kaftancioglu/alternative_scenario/fastQC_as/output/h3k4me1_ms_nk_ENCSR482UGV/*_fastqc.zip
multiqc -o /cta/users/merve.kaftancioglu/alternative_scenario/multiqc_data/h3k27ac_ms_nk_ENCSR246TTM /cta/users/merve.kaftancioglu/alternative_scenario/fastQC_as/output/h3k27ac_ms_nk_ENCSR246TTM/*_fastqc.zip
multiqc -o /cta/users/merve.kaftancioglu/alternative_scenario/multiqc_data/h3k4me1_ms_nk_ENCSR718LCI /cta/users/merve.kaftancioglu/alternative_scenario/fastQC_as/output/h3k4me1_ms_nk_ENCSR718LCI/*_fastqc.zip
multiqc -o /cta/users/merve.kaftancioglu/alternative_scenario/multiqc_data/h3k27ac_ms_nk_ENCSR469CFL /cta/users/merve.kaftancioglu/alternative_scenario/fastQC_as/output/h3k27ac_ms_nk_ENCSR469CFL/*_fastqc.zip	 
multiqc -o /cta/users/merve.kaftancioglu/alternative_scenario/multiqc_data/h3k4me1_ms_nk_ENCSR806HUT /cta/users/merve.kaftancioglu/alternative_scenario/fastQC_as/output/h3k4me1_ms_nk_ENCSR806HUT/*_fastqc.zip
multiqc -o /cta/users/merve.kaftancioglu/alternative_scenario/multiqc_data/h3k27ac_ms_nk_ENCSR746AIX /cta/users/merve.kaftancioglu/alternative_scenario/fastQC_as/output/h3k27ac_ms_nk_ENCSR746AIX/*_fastqc.zip	 
multiqc -o /cta/users/merve.kaftancioglu/alternative_scenario/multiqc_data/h3k4me1_normal_b_ENCSR156BXM /cta/users/merve.kaftancioglu/alternative_scenario/fastQC_as/output/h3k4me1_normal_b_ENCSR156BXM/*_fastqc.zip
multiqc -o /cta/users/merve.kaftancioglu/alternative_scenario/multiqc_data/h3k27ac_normal_b_ENCSR685KZA /cta/users/merve.kaftancioglu/alternative_scenario/fastQC_as/output/h3k27ac_normal_b_ENCSR685KZA/*_fastqc.zip
multiqc -o /cta/users/merve.kaftancioglu/alternative_scenario/multiqc_data/h3k4me1_normal_cd4_ENCSR102SOR /cta/users/merve.kaftancioglu/alternative_scenario/fastQC_as/output/h3k4me1_normal_cd4_ENCSR102SOR/*_fastqc.zip
multiqc -o /cta/users/merve.kaftancioglu/alternative_scenario/multiqc_data/h3k27ac_normal_cd4_ENCSR819NCZ /cta/users/merve.kaftancioglu/alternative_scenario/fastQC_as/output/h3k27ac_normal_cd4_ENCSR819NCZ/*_fastqc.zip
multiqc -o /cta/users/merve.kaftancioglu/alternative_scenario/multiqc_data/h3k4me1_normal_cd8_ENCSR217SHH /cta/users/merve.kaftancioglu/alternative_scenario/fastQC_as/output/h3k4me1_normal_cd8_ENCSR217SHH/*_fastqc.zip
multiqc -o /cta/users/merve.kaftancioglu/alternative_scenario/multiqc_data/h3k27ac_normal_cd8_ENCSR976RWL /cta/users/merve.kaftancioglu/alternative_scenario/fastQC_as/output/h3k27ac_normal_cd8_ENCSR976RWL/*_fastqc.zip 
multiqc -o /cta/users/merve.kaftancioglu/alternative_scenario/multiqc_data/h3k4me1_normal_nk_ENCSR277YKG /cta/users/merve.kaftancioglu/alternative_scenario/fastQC_as/output/h3k4me1_normal_nk_ENCSR277YKG/*_fastqc.zip
multiqc -o /cta/users/merve.kaftancioglu/alternative_scenario/multiqc_data/h3k27ac_normal_nk_ENCSR977FMZ /cta/users/merve.kaftancioglu/alternative_scenario/fastQC_as/output/h3k27ac_normal_nk_ENCSR977FMZ/*_fastqc.zip	 
multiqc -o /cta/users/merve.kaftancioglu/alternative_scenario/multiqc_data/h3k4me3_ms_b_ENCSR260CRI /cta/users/merve.kaftancioglu/alternative_scenario/fastQC_as/output/h3k4me3_ms_b_ENCSR260CRI/*_fastqc.zip
multiqc -o /cta/users/merve.kaftancioglu/alternative_scenario/multiqc_data/output/h3k27me3_ms_b_ENCSR009ZRH /cta/users/merve.kaftancioglu/alternative_scenario/fastQC_as/output/h3k27me3_ms_b_ENCSR009ZRH/*_fastqc.zip
multiqc -o /cta/users/merve.kaftancioglu/alternative_scenario/multiqc_data/h3k27me3_ms_b_ENCSR182NLA /cta/users/merve.kaftancioglu/alternative_scenario/fastQC_as/output/h3k27me3_ms_b_ENCSR182NLA/*_fastqc.zip
multiqc -o /cta/users/merve.kaftancioglu/alternative_scenario/multiqc_data/h3k4me3_ms_b_ENCSR530AVK /cta/users/merve.kaftancioglu/alternative_scenario/fastQC_as/output/h3k4me3_ms_b_ENCSR530AVK/*_fastqc.zip
multiqc -o /cta/users/merve.kaftancioglu/alternative_scenario/multiqc_data/h3k27me3_ms_b_ENCSR272YVX /cta/users/merve.kaftancioglu/alternative_scenario/fastQC_as/output/h3k27me3_ms_b_ENCSR272YVX/*_fastqc.zip	 
multiqc -o /cta/users/merve.kaftancioglu/alternative_scenario/multiqc_data/h3k4me3_ms_b_ENCSR713ZYF /cta/users/merve.kaftancioglu/alternative_scenario/fastQC_as/output/h3k4me3_ms_b_ENCSR713ZYF/*_fastqc.zip
multiqc -o /cta/users/merve.kaftancioglu/alternative_scenario/multiqc_data/h3k27me3_ms_b_ENCSR649FUX /cta/users/merve.kaftancioglu/alternative_scenario/fastQC_as/output/h3k27me3_ms_b_ENCSR649FUX/*_fastqc.zip	 
multiqc -o /cta/users/merve.kaftancioglu/alternative_scenario/multiqc_data/h3k4me3_ms_b_ENCSR923EGO /cta/users/merve.kaftancioglu/alternative_scenario/fastQC_as/output/h3k4me3_ms_b_ENCSR923EGO/*_fastqc.zip
multiqc -o /cta/users/merve.kaftancioglu/alternative_scenario/multiqc_data/h3k27me3_ms_b_ENCSR842WWX /cta/users/merve.kaftancioglu/alternative_scenario/fastQC_as/output/h3k27me3_ms_b_ENCSR842WWX/*_fastqc.zip	 
multiqc -o /cta/users/merve.kaftancioglu/alternative_scenario/multiqc_data/h3k4me3_ms_cd4_ENCSR180NCM /cta/users/merve.kaftancioglu/alternative_scenario/fastQC_as/output/h3k4me3_ms_cd4_ENCSR180NCM/*_fastqc.zip
multiqc -o /cta/users/merve.kaftancioglu/alternative_scenario/multiqc_data/h3k27me3_ms_cd4_ENCSR277XYX /cta/users/merve.kaftancioglu/alternative_scenario/fastQC_as/output/h3k27me3_ms_cd4_ENCSR277XYX/*_fastqc.zip
multiqc -o /cta/users/merve.kaftancioglu/alternative_scenario/multiqc_data/h3k4me3_ms_cd4_ENCSR341QLC /cta/users/merve.kaftancioglu/alternative_scenario/fastQC_as/output/h3k4me3_ms_cd4_ENCSR341QLC/*_fastqc.zip
multiqc -o /cta/users/merve.kaftancioglu/alternative_scenario/multiqc_data/h3k27me3_ms_cd4_ENCSR526TNC /cta/users/merve.kaftancioglu/alternative_scenario/fastQC_as/output/h3k27me3_ms_cd4_ENCSR526TNC/*_fastqc.zip	 
multiqc -o /cta/users/merve.kaftancioglu/alternative_scenario/multiqc_data/h3k4me3_ms_cd4_ENCSR482TGI /cta/users/merve.kaftancioglu/alternative_scenario/fastQC_as/output/h3k4me3_ms_cd4_ENCSR482TGI/*_fastqc.zip
multiqc -o /cta/users/merve.kaftancioglu/alternative_scenario/multiqc_data/h3k27me3_ms_cd4_ENCSR592EKF /cta/users/merve.kaftancioglu/alternative_scenario/fastQC_as/output/h3k27me3_ms_cd4_ENCSR592EKF/*_fastqc.zip	 
multiqc -o /cta/users/merve.kaftancioglu/alternative_scenario/multiqc_data/h3k4me3_ms_cd4_ENCSR486XJK /cta/users/merve.kaftancioglu/alternative_scenario/fastQC_as/output/h3k4me3_ms_cd4_ENCSR486XJK/*_fastqc.zip
multiqc -o /cta/users/merve.kaftancioglu/alternative_scenario/multiqc_data/h3k27me3_ms_cd4_ENCSR613UFD /cta/users/merve.kaftancioglu/alternative_scenario/fastQC_as/output/h3k27me3_ms_cd4_ENCSR613UFD/*_fastqc.zip	 
multiqc -o /cta/users/merve.kaftancioglu/alternative_scenario/multiqc_data/h3k4me3_ms_cd4_ENCSR496LKR /cta/users/merve.kaftancioglu/alternative_scenario/fastQC_as/output/h3k4me3_ms_cd4_ENCSR496LKR/*_fastqc.zip
multiqc -o /cta/users/merve.kaftancioglu/alternative_scenario/multiqc_data/h3k27me3_ms_cd4_ENCSR740SDR /cta/users/merve.kaftancioglu/alternative_scenario/fastQC_as/output/h3k27me3_ms_cd4_ENCSR740SDR/*_fastqc.zip
multiqc -o /cta/users/merve.kaftancioglu/alternative_scenario/multiqc_data/h3k4me3_ms_cd4_ENCSR603LTN /cta/users/merve.kaftancioglu/alternative_scenario/fastQC_as/output/h3k4me3_ms_cd4_ENCSR603LTN/*_fastqc.zip
multiqc -o /cta/users/merve.kaftancioglu/alternative_scenario/multiqc_data/h3k27me3_ms_cd4_ENCSR779JLY /cta/users/merve.kaftancioglu/alternative_scenario/fastQC_as/output/h3k27me3_ms_cd4_ENCSR779JLY/*_fastqc.zip	 
multiqc -o /cta/users/merve.kaftancioglu/alternative_scenario/multiqc_data/h3k4me3_ms_cd4_ENCSR802MXQ /cta/users/merve.kaftancioglu/alternative_scenario/fastQC_as/output/h3k4me3_ms_cd4_ENCSR802MXQ/*_fastqc.zip
multiqc -o /cta/users/merve.kaftancioglu/alternative_scenario/multiqc_data/h3k27me3_ms_cd4_ENCSR993CTA /cta/users/merve.kaftancioglu/alternative_scenario/fastQC_as/output/h3k27me3_ms_cd4_ENCSR993CTA/*_fastqc.zip	 
multiqc -o /cta/users/merve.kaftancioglu/alternative_scenario/multiqc_data/h3k4me3_ms_cd4_ENCSR878YHM /cta/users/merve.kaftancioglu/alternative_scenario/fastQC_as/output/h3k4me3_ms_cd4_ENCSR878YHM/*_fastqc.zip
multiqc -o /cta/users/merve.kaftancioglu/alternative_scenario/multiqc_data/h3k27me3_ms_cd8_ENCSR116FVG /cta/users/merve.kaftancioglu/alternative_scenario/fastQC_as/output/h3k27me3_ms_cd8_ENCSR116FVG/*_fastqc.zip	 
multiqc -o /cta/users/merve.kaftancioglu/alternative_scenario/multiqc_data/h3k4me3_ms_cd4_ENCSR954ZLD /cta/users/merve.kaftancioglu/alternative_scenario/fastQC_as/output/h3k4me3_ms_cd4_ENCSR954ZLD/*_fastqc.zip
multiqc -o /cta/users/merve.kaftancioglu/alternative_scenario/multiqc_data/h3k27me3_ms_cd8_ENCSR122JCM /cta/users/merve.kaftancioglu/alternative_scenario/fastQC_as/output/h3k27me3_ms_cd8_ENCSR122JCM/*_fastqc.zip	 
multiqc -o /cta/users/merve.kaftancioglu/alternative_scenario/multiqc_data/h3k4me3_ms_cd8_ENCSR231ZZH /cta/users/merve.kaftancioglu/alternative_scenario/fastQC_as/output/h3k4me3_ms_cd8_ENCSR231ZZH/*_fastqc.zip
multiqc -o /cta/users/merve.kaftancioglu/alternative_scenario/multiqc_data/h3k27me3_ms_cd8_ENCSR216ZVA /cta/users/merve.kaftancioglu/alternative_scenario/fastQC_as/output/h3k27me3_ms_cd8_ENCSR216ZVA/*_fastqc.zip	 
multiqc -o /cta/users/merve.kaftancioglu/alternative_scenario/multiqc_data/h3k4me3_ms_cd8_ENCSR278QHR /cta/users/merve.kaftancioglu/alternative_scenario/fastQC_as/output/h3k4me3_ms_cd8_ENCSR278QHR/*_fastqc.zip
multiqc -o /cta/users/merve.kaftancioglu/alternative_scenario/multiqc_data/h3k27me3_ms_cd8_ENCSR284IKS /cta/users/merve.kaftancioglu/alternative_scenario/fastQC_as/output/h3k27me3_ms_cd8_ENCSR284IKS/*_fastqc.zip	 
multiqc -o /cta/users/merve.kaftancioglu/alternative_scenario/multiqc_data /h3k4me3_ms_cd8_ENCSR516CKJ/cta/users/merve.kaftancioglu/alternative_scenario/fastQC_as/output/h3k4me3_ms_cd8_ENCSR516CKJ/*_fastqc.zip
multiqc -o /cta/users/merve.kaftancioglu/alternative_scenario/multiqc_data/h3k27me3_ms_cd8_ENCSR385BOZ /cta/users/merve.kaftancioglu/alternative_scenario/fastQC_as/output/h3k27me3_ms_cd8_ENCSR385BOZ/*_fastqc.zip	 
multiqc -o /cta/users/merve.kaftancioglu/alternative_scenario/multiqc_data/h3k4me3_ms_cd8_ENCSR535YYH /cta/users/merve.kaftancioglu/alternative_scenario/fastQC_as/output/h3k4me3_ms_cd8_ENCSR535YYH/*_fastqc.zip
multiqc -o /cta/users/merve.kaftancioglu/alternative_scenario/multiqc_data/h3k27me3_ms_cd8_ENCSR521SFR /cta/users/merve.kaftancioglu/alternative_scenario/fastQC_as/output/h3k27me3_ms_cd8_ENCSR521SFR/*_fastqc.zip	 
multiqc -o /cta/users/merve.kaftancioglu/alternative_scenario/multiqc_data/h3k4me3_ms_cd8_ENCSR741XAE /cta/users/merve.kaftancioglu/alternative_scenario/fastQC_as/output/h3k4me3_ms_cd8_ENCSR741XAE/*_fastqc.zip
multiqc -o /cta/users/merve.kaftancioglu/alternative_scenario/multiqc_data/h3k27me3_ms_nk_ENCSR469QVG /cta/users/merve.kaftancioglu/alternative_scenario/fastQC_as/output/h3k27me3_ms_nk_ENCSR469QVG/*_fastqc.zip 
multiqc -o /cta/users/merve.kaftancioglu/alternative_scenario/multiqc_data/h3k4me3_ms_cd8_ENCSR848XJL /cta/users/merve.kaftancioglu/alternative_scenario/fastQC_as/output/h3k4me3_ms_cd8_ENCSR848XJL/*_fastqc.zip
multiqc -o /cta/users/merve.kaftancioglu/alternative_scenario/multiqc_data/h3k27me3_ms_nk_ENCSR565WDW /cta/users/merve.kaftancioglu/alternative_scenario/fastQC_as/output/h3k27me3_ms_nk_ENCSR565WDW/*_fastqc.zip	 
multiqc -o /cta/users/merve.kaftancioglu/alternative_scenario/multiqc_data/h3k4me3_ms_nk_ENCSR234BAD /cta/users/merve.kaftancioglu/alternative_scenario/fastQC_as/output/h3k4me3_ms_nk_ENCSR234BAD/*_fastqc.zip
multiqc -o /cta/users/merve.kaftancioglu/alternative_scenario/multiqc_data/h3k27me3_normal_b_ENCSR589LHR /cta/users/merve.kaftancioglu/alternative_scenario/fastQC_as/output/h3k27me3_normal_b_ENCSR589LHR/*_fastqc.zip	 
multiqc -o /cta/users/merve.kaftancioglu/alternative_scenario/multiqc_data/h3k4me3_ms_nk_ENCSR703TWY /cta/users/merve.kaftancioglu/alternative_scenario/fastQC_as/output/h3k4me3_ms_nk_ENCSR703TWY/*_fastqc.zip
multiqc -o /cta/users/merve.kaftancioglu/alternative_scenario/multiqc_data/h3k27me3_normal_cd4_ENCSR068YVZ /cta/users/merve.kaftancioglu/alternative_scenario/fastQC_as/output/h3k27me3_normal_cd4_ENCSR068YVZ/*_fastqc.zip
multiqc -o /cta/users/merve.kaftancioglu/alternative_scenario/multiqc_data/h3k4me3_ms_nk_ENCSR755ZGS /cta/users/merve.kaftancioglu/alternative_scenario/fastQC_as/output/h3k4me3_ms_nk_ENCSR755ZGS/*_fastqc.zip
multiqc -o /cta/users/merve.kaftancioglu/alternative_scenario/multiqc_data/h3k27me3_normal_cd8_ENCSR720QRX /cta/users/merve.kaftancioglu/alternative_scenario/fastQC_as/output/h3k27me3_normal_cd8_ENCSR720QRX/*_fastqc.zip  
multiqc -o /cta/users/merve.kaftancioglu/alternative_scenario/multiqc_data/h3k4me3_normal_b_ENCSR791CAF /cta/users/merve.kaftancioglu/alternative_scenario/fastQC_as/output/h3k4me3_normal_b_ENCSR791CAF/*_fastqc.zip
multiqc -o /cta/users/merve.kaftancioglu/alternative_scenario/multiqc_data/h3k27me3_normal_nk_ENCSR639NIG /cta/users/merve.kaftancioglu/alternative_scenario/fastQC_as/output/h3k27me3_normal_nk_ENCSR639NIG/*_fastqc.zip	 
multiqc -o /cta/users/merve.kaftancioglu/alternative_scenario/multiqc_data/h3k4me3_normal_cd4_ENCSR537KJA /cta/users/merve.kaftancioglu/alternative_scenario/fastQC_as/output/h3k4me3_normal_cd4_ENCSR537KJA/*_fastqc.zip
multiqc -o /cta/users/merve.kaftancioglu/alternative_scenario/multiqc_data/h3k36me3_ms_b_ENCSR089VQL /cta/users/merve.kaftancioglu/alternative_scenario/fastQC_as/output/h3k36me3_ms_b_ENCSR089VQL/*_fastqc.zip
multiqc -o /cta/users/merve.kaftancioglu/alternative_scenario/multiqc_data/h3k4me3_normal_cd8_ENCSR123ZAT /cta/users/merve.kaftancioglu/alternative_scenario/fastQC_as/output/h3k4me3_normal_cd8_ENCSR123ZAT/*_fastqc.zip
multiqc -o /cta/users/merve.kaftancioglu/alternative_scenario/multiqc_data/h3k36me3_ms_b_ENCSR119PSR /cta/users/merve.kaftancioglu/alternative_scenario/fastQC_as/output/h3k36me3_ms_b_ENCSR119PSR/*_fastqc.zip
multiqc -o /cta/users/merve.kaftancioglu/alternative_scenario/multiqc_data/h3k4me3_normal_nk_ENCSR394JFQ /cta/users/merve.kaftancioglu/alternative_scenario/fastQC_as/output/h3k4me3_normal_nk_ENCSR394JFQ/*_fastqc.zip
multiqc -o /cta/users/merve.kaftancioglu/alternative_scenario/multiqc_data/h3k36me3_ms_b_ENCSR238WFK /cta/users/merve.kaftancioglu/alternative_scenario/fastQC_as/output/h3k36me3_ms_b_ENCSR238WFK/*_fastqc.zip
multiqc -o /cta/users/merve.kaftancioglu/alternative_scenario/multiqc_data/h3k9me3_ms_b_ENCSR445DNH /cta/users/merve.kaftancioglu/alternative_scenario/fastQC_as/output/h3k9me3_ms_b_ENCSR445DNH/*_fastqc.zip
multiqc -o /cta/users/merve.kaftancioglu/alternative_scenario/multiqc_data/h3k36me3_ms_b_ENCSR987OPY /cta/users/merve.kaftancioglu/alternative_scenario/fastQC_as/output/h3k36me3_ms_b_ENCSR987OPY/*_fastqc.zip	 
multiqc -o /cta/users/merve.kaftancioglu/alternative_scenario/multiqc_data/h3k9me3_ms_b_ENCSR486SZN /cta/users/merve.kaftancioglu/alternative_scenario/fastQC_as/output/h3k9me3_ms_b_ENCSR486SZN/*_fastqc.zip
multiqc -o /cta/users/merve.kaftancioglu/alternative_scenario/multiqc_data/h3k36me3_ms_cd4_ENCSR276NGH/ /cta/users/merve.kaftancioglu/alternative_scenario/fastQC_as/output/h3k36me3_ms_cd4_ENCSR276NGH/*_fastqc.zip	 
multiqc -o /cta/users/merve.kaftancioglu/alternative_scenario/multiqc_data/h3k9me3_ms_b_ENCSR983MKX /cta/users/merve.kaftancioglu/alternative_scenario/fastQC_as/output/h3k9me3_ms_b_ENCSR983MKX/*_fastqc.zip
multiqc -o /cta/users/merve.kaftancioglu/alternative_scenario/multiqc_data/h3k36me3_ms_cd4_ENCSR330CQU /cta/users/merve.kaftancioglu/alternative_scenario/fastQC_as/output/h3k36me3_ms_cd4_ENCSR330CQU/*_fastqc.zip
multiqc -o /cta/users/merve.kaftancioglu/alternative_scenario/multiqc_data/h3k9me3_ms_cd4_ENCSR057IZD /cta/users/merve.kaftancioglu/alternative_scenario/fastQC_as/output/h3k9me3_ms_cd4_ENCSR057IZD/*_fastqc.zip
multiqc -o /cta/users/merve.kaftancioglu/alternative_scenario/multiqc_data/h3k36me3_ms_cd4_ENCSR471WZE /cta/users/merve.kaftancioglu/alternative_scenario/fastQC_as/output/h3k36me3_ms_cd4_ENCSR471WZE/*_fastqc.zip	 
multiqc -o /cta/users/merve.kaftancioglu/alternative_scenario/multiqc_data/h3k9me3_ms_cd4_ENCSR550DPT /cta/users/merve.kaftancioglu/alternative_scenario/fastQC_as/output/h3k9me3_ms_cd4_ENCSR550DPT/*_fastqc.zip
multiqc -o /cta/users/merve.kaftancioglu/alternative_scenario/multiqc_data/h3k36me3_ms_cd4_ENCSR473TGK /cta/users/merve.kaftancioglu/alternative_scenario/fastQC_as/output/h3k36me3_ms_cd4_ENCSR473TGK/*_fastqc.zip	 
multiqc -o /cta/users/merve.kaftancioglu/alternative_scenario/multiqc_data/h3k9me3_ms_cd4_ENCSR677OEF /cta/users/merve.kaftancioglu/alternative_scenario/fastQC_as/output/h3k9me3_ms_cd4_ENCSR677OEF/*_fastqc.zip
multiqc -o /cta/users/merve.kaftancioglu/alternative_scenario/multiqc_data/h3k36me3_ms_cd4_ENCSR482VIB /cta/users/merve.kaftancioglu/alternative_scenario/fastQC_as/output/h3k36me3_ms_cd4_ENCSR482VIB/*_fastqc.zip	 
multiqc -o /cta/users/merve.kaftancioglu/alternative_scenario/multiqc_data/h3k9me3_ms_cd4_ENCSR729ZXH /cta/users/merve.kaftancioglu/alternative_scenario/fastQC_as/output/h3k9me3_ms_cd4_ENCSR729ZXH/*_fastqc.zip
multiqc -o /cta/users/merve.kaftancioglu/alternative_scenario/multiqc_data/h3k36me3_ms_cd4_ENCSR532PXR /cta/users/merve.kaftancioglu/alternative_scenario/fastQC_as/output/h3k36me3_ms_cd4_ENCSR532PXR/*_fastqc.zip
multiqc -o /cta/users/merve.kaftancioglu/alternative_scenario/multiqc_data/h3k9me3_ms_cd4_ENCSR851RJV /cta/users/merve.kaftancioglu/alternative_scenario/fastQC_as/output/h3k9me3_ms_cd4_ENCSR851RJV/*_fastqc.zip
multiqc -o /cta/users/merve.kaftancioglu/alternative_scenario/multiqc_data/h3k36me3_ms_cd4_ENCSR785PKM /cta/users/merve.kaftancioglu/alternative_scenario/fastQC_as/output/h3k36me3_ms_cd4_ENCSR785PKM/*_fastqc.zip	 
multiqc -o /cta/users/merve.kaftancioglu/alternative_scenario/multiqc_data/h3k9me3_ms_cd4_ENCSR919DFZ /cta/users/merve.kaftancioglu/alternative_scenario/fastQC_as/output/h3k9me3_ms_cd4_ENCSR919DFZ/*_fastqc.zip
multiqc -o /cta/users/merve.kaftancioglu/alternative_scenario/multiqc_data/h3k36me3_ms_cd4_ENCSR865FTW /cta/users/merve.kaftancioglu/alternative_scenario/fastQC_as/output/h3k36me3_ms_cd4_ENCSR865FTW/*_fastqc.zip	 
multiqc -o /cta/users/merve.kaftancioglu/alternative_scenario/multiqc_data/h3k9me3_ms_cd4_ENCSR953STZ /cta/users/merve.kaftancioglu/alternative_scenario/fastQC_as/output/h3k9me3_ms_cd4_ENCSR953STZ/*_fastqc.zip
multiqc -o /cta/users/merve.kaftancioglu/alternative_scenario/multiqc_data/h3k36me3_ms_cd4_ENCSR898VJE /cta/users/merve.kaftancioglu/alternative_scenario/fastQC_as/output/h3k36me3_ms_cd4_ENCSR898VJE/*_fastqc.zip	 
multiqc -o /cta/users/merve.kaftancioglu/alternative_scenario/multiqc_data/h3k9me3_ms_cd4_ENCSR959VZU /cta/users/merve.kaftancioglu/alternative_scenario/fastQC_as/output/h3k9me3_ms_cd4_ENCSR959VZU/*_fastqc.zip
multiqc -o /cta/users/merve.kaftancioglu/alternative_scenario/multiqc_data/output/h3k36me3_ms_cd8_ENCSR239OWL /cta/users/merve.kaftancioglu/alternative_scenario/fastQC_as/output/h3k36me3_ms_cd8_ENCSR239OWL/*_fastqc.zip	 
multiqc -o /cta/users/merve.kaftancioglu/alternative_scenario/multiqc_data/h3k9me3_ms_cd8_ENCSR101USF /cta/users/merve.kaftancioglu/alternative_scenario/fastQC_as/output/h3k9me3_ms_cd8_ENCSR101USF/*_fastqc.zip
multiqc -o /cta/users/merve.kaftancioglu/alternative_scenario/multiqc_data/h3k36me3_ms_cd8_ENCSR435MSB /cta/users/merve.kaftancioglu/alternative_scenario/fastQC_as/output/h3k36me3_ms_cd8_ENCSR435MSB/*_fastqc.zip	 
multiqc -o /cta/users/merve.kaftancioglu/alternative_scenario/multiqc_data/h3k9me3_ms_cd8_ENCSR354GNT /cta/users/merve.kaftancioglu/alternative_scenario/fastQC_as/output/h3k9me3_ms_cd8_ENCSR354GNT/*_fastqc.zip
multiqc -o /cta/users/merve.kaftancioglu/alternative_scenario/multiqc_data/h3k36me3_ms_cd8_ENCSR631VIW /cta/users/merve.kaftancioglu/alternative_scenario/fastQC_as/output/h3k36me3_ms_cd8_ENCSR631VIW/*_fastqc.zip	 
multiqc -o /cta/users/merve.kaftancioglu/alternative_scenario/multiqc_data/h3k9me3_ms_cd8_ENCSR377QYB /cta/users/merve.kaftancioglu/alternative_scenario/fastQC_as/output/h3k9me3_ms_cd8_ENCSR377QYB/*_fastqc.zip
multiqc -o /cta/users/merve.kaftancioglu/alternative_scenario/multiqc_data/h3k36me3_ms_cd8_ENCSR652WFO /cta/users/merve.kaftancioglu/alternative_scenario/fastQC_as/output/h3k36me3_ms_cd8_ENCSR652WFO/*_fastqc.zip	 
multiqc -o /cta/users/merve.kaftancioglu/alternative_scenario/multiqc_data/h3k9me3_ms_cd8_ENCSR733LCG /cta/users/merve.kaftancioglu/alternative_scenario/fastQC_as/output/h3k9me3_ms_cd8_ENCSR733LCG/*_fastqc.zip
multiqc -o /cta/users/merve.kaftancioglu/alternative_scenario/multiqc_data/h3k36me3_ms_cd8_ENCSR757FGN /cta/users/merve.kaftancioglu/alternative_scenario/fastQC_as/output/h3k36me3_ms_cd8_ENCSR757FGN/*_fastqc.zip	 
multiqc -o /cta/users/merve.kaftancioglu/alternative_scenario/multiqc_data/h3k9me3_ms_cd8_ENCSR980KTW /cta/users/merve.kaftancioglu/alternative_scenario/fastQC_as/output/h3k9me3_ms_cd8_ENCSR980KTW/*_fastqc.zip
multiqc -o /cta/users/merve.kaftancioglu/alternative_scenario/multiqc_data/h3k36me3_ms_nk_ENCSR158VSE /cta/users/merve.kaftancioglu/alternative_scenario/fastQC_as/output/h3k36me3_ms_nk_ENCSR158VSE/*_fastqc.zip 
multiqc -o /cta/users/merve.kaftancioglu/alternative_scenario/multiqc_data/h3k9me3_ms_nk_ENCSR061ATV /cta/users/merve.kaftancioglu/alternative_scenario/fastQC_as/output/h3k9me3_ms_nk_ENCSR061ATV/*_fastqc.zip
multiqc -o /cta/users/merve.kaftancioglu/alternative_scenario/multiqc_data/h3k36me3_ms_nk_ENCSR245KON /cta/users/merve.kaftancioglu/alternative_scenario/fastQC_as/output/h3k36me3_ms_nk_ENCSR245KON/*_fastqc.zip	 
multiqc -o /cta/users/merve.kaftancioglu/alternative_scenario/multiqc_data/h3k9me3_ms_nk_ENCSR611YIA /cta/users/merve.kaftancioglu/alternative_scenario/fastQC_as/output/h3k9me3_ms_nk_ENCSR611YIA/*_fastqc.zip
multiqc -o /cta/users/merve.kaftancioglu/alternative_scenario/multiqc_data/h3k36me3_ms_nk_ENCSR530YDY /cta/users/merve.kaftancioglu/alternative_scenario/fastQC_as/output/h3k36me3_ms_nk_ENCSR530YDY/*_fastqc.zip	 
multiqc -o /cta/users/merve.kaftancioglu/alternative_scenario/multiqc_data/h3k9me3_ms_nk_ENCSR667HUI /cta/users/merve.kaftancioglu/alternative_scenario/fastQC_as/output/h3k9me3_ms_nk_ENCSR667HUI/*_fastqc.zip
multiqc -o /cta/users/merve.kaftancioglu/alternative_scenario/multiqc_data/h3k36me3_normal_b_ENCSR831AXK /cta/users/merve.kaftancioglu/alternative_scenario/fastQC_as/output/h3k36me3_normal_b_ENCSR831AXK/*_fastqc.zip	 
multiqc -o /cta/users/merve.kaftancioglu/alternative_scenario/multiqc_data/h3k9me3_normal_b_ENCSR445LTM /cta/users/merve.kaftancioglu/alternative_scenario/fastQC_as/output/h3k9me3_normal_b_ENCSR445LTM/*_fastqc.zip
multiqc -o /cta/users/merve.kaftancioglu/alternative_scenario/multiqc_data/h3k36me3_normal_cd4_ENCSR864OKB /cta/users/merve.kaftancioglu/alternative_scenario/fastQC_as/output/h3k36me3_normal_cd4_ENCSR864OKB/*_fastqc.zip
multiqc -o /cta/users/merve.kaftancioglu/alternative_scenario/multiqc_data/h3k9me3_normal_cd4_ENCSR433EWI /cta/users/merve.kaftancioglu/alternative_scenario/fastQC_as/output/h3k9me3_normal_cd4_ENCSR433EWI/*_fastqc.zip
multiqc -o /cta/users/merve.kaftancioglu/alternative_scenario/multiqc_data/h3k36me3_normal_cd8_ENCSR303SQG /cta/users/merve.kaftancioglu/alternative_scenario/fastQC_as/output/h3k36me3_normal_cd8_ENCSR303SQG/*_fastqc.zip
multiqc -o /cta/users/merve.kaftancioglu/alternative_scenario/multiqc_data/h3k9me3_normal_cd8_ENCSR294HTM /cta/users/merve.kaftancioglu/alternative_scenario/fastQC_as/output/h3k9me3_normal_cd8_ENCSR294HTM/*_fastqc.zip
multiqc -o /cta/users/merve.kaftancioglu/alternative_scenario/multiqc_data/h3k36me3_normal_nk_ENCSR056GJY /cta/users/merve.kaftancioglu/alternative_scenario/fastQC_as/output/h3k36me3_normal_nk_ENCSR056GJY/*_fastqc.zip	 
multiqc -o /cta/users/merve.kaftancioglu/alternative_scenario/multiqc_data/h3k9me3_normal_nk_ENCSR025UNZ /cta/users/merve.kaftancioglu/alternative_scenario/fastQC_as/output/h3k9me3_normal_nk_ENCSR025UNZ/*_fastqc.zip
multiqc -o /cta/users/merve.kaftancioglu/alternative_scenario/multiqc_data/h3k4me1_ms_b_ENCSR084FXT /cta/users/merve.kaftancioglu/alternative_scenario/fastQC_as/output/h3k4me1_ms_b_ENCSR084FXT/*_fastqc.zip

```

# Trimmomatic (divided into ~50 files):

```
#!/bin/bash
#SBATCH --job-name=trimmomatic_single
#SBATCH --nodes=1
#SBATCH --account=users
#SBATCH --qos=mid_mdbf
#SBATCH --partition=mid_mdbf
#SBATCH --cpus-per-task=16
#SBATCH -o /cta/users/merve.kaftancioglu/alternative_scenario/trimmomatic_as/out-err/%x-%j-slurm.out
#SBATCH -e /cta/users/merve.kaftancioglu/alternative_scenario/trimmomatic_as/out-err/%x-%j-slurm.err

# Load FastQC module
module load trimmomatic-0.38-gcc-9.2.0-3evdnxc

# List of sample IDs (folders)
sample_ids=("h3k27ac_ms_b_ENCSR295QZX" "h3k27ac_ms_b_ENCSR538URI" "h3k27ac_ms_b_ENCSR617GMR" "h3k4me1_ms_b_ENCSR084FXT" "h3k4me1_ms_b_ENCSR480ITK" "h3k4me1_ms_b_ENCSR759DHN" "h3k4me1_ms_b_ENCSR900IAM" "h3k4me1_ms_b_ENCSR931JYC" "h3k4me3_ms_b_ENCSR260CRI" "h3k4me3_ms_b_ENCSR461QMZ" "h3k4me3_ms_b_ENCSR530AVK" "h3k4me3_ms_b_ENCSR713ZYF" "h3k4me3_ms_b_ENCSR923EGO" "h3k9me3_ms_b_ENCSR445DNH" "h3k9me3_ms_b_ENCSR486SZN" "h3k9me3_ms_b_ENCSR983MKX" "h3k27ac_ms_b_ENCSR364OIK" "h3k27ac_ms_b_ENCSR541PLO" "h3k36me3_ms_cd4_ENCSR471WZE" "h3k36me3_ms_cd4_ENCSR473TGK" "h3k4me1_ms_nk_ENCSR482UGV" "h3k4me1_ms_nk_ENCSR718LCI" "h3k4me1_ms_nk_ENCSR806HUT" "h3k4me3_ms_nk_ENCSR234BAD" "h3k4me3_ms_nk_ENCSR703TWY" "h3k4me3_ms_nk_ENCSR755ZGS" "h3k9me3_ms_nk_ENCSR061ATV" "h3k9me3_ms_nk_ENCSR611YIA" "h3k9me3_ms_nk_ENCSR667HUI" "h3k27ac_ms_nk_ENCSR246TTM" "h3k27ac_ms_nk_ENCSR469CFL" "h3k27ac_ms_nk_ENCSR746AIX" "h3k27me3_ms_b_ENCSR009ZRH" "h3k27me3_ms_b_ENCSR182NLA" "h3k27me3_ms_b_ENCSR272YVX" "h3k27me3_ms_b_ENCSR649FUX" "h3k27me3_ms_b_ENCSR842WWX" "h3k36me3_ms_b_ENCSR089VQL" "h3k36me3_ms_b_ENCSR119PSR" "h3k36me3_ms_b_ENCSR238WFK" "h3k36me3_ms_b_ENCSR987OPY" "h3k4me1_ms_cd4_ENCSR036YVP" "h3k4me1_ms_cd4_ENCSR043JIA" "h3k4me1_ms_cd4_ENCSR093VUP" "h3k4me1_ms_cd4_ENCSR124IGC" "h3k4me1_ms_cd4_ENCSR315MYP" "h3k4me1_ms_cd4_ENCSR458QEY" "h3k4me1_ms_cd4_ENCSR485FBT" "h3k4me1_ms_cd4_ENCSR641ZFV" "h3k4me1_ms_cd4_ENCSR687CQX" "h3k4me1_ms_cd8_ENCSR231XAP" "h3k4me1_ms_cd8_ENCSR572XTB" "h3k4me1_ms_cd8_ENCSR788TEF" "h3k4me1_ms_cd8_ENCSR815YZL" "h3k4me1_ms_cd8_ENCSR861FEC" "h3k4me1_ms_cd8_ENCSR940PHE" "h3k4me3_ms_cd4_ENCSR180NCM" "h3k4me3_ms_cd4_ENCSR341QLC" "h3k4me3_ms_cd4_ENCSR482TGI" "h3k4me3_ms_cd4_ENCSR486XJK" "h3k4me3_ms_cd4_ENCSR496LKR" "h3k4me3_ms_cd4_ENCSR603LTN" "h3k4me3_ms_cd4_ENCSR802MXQ" "h3k4me3_ms_cd4_ENCSR878YHM" "h3k4me3_ms_cd4_ENCSR954ZLD" "h3k4me3_ms_cd8_ENCSR231ZZH" "h3k4me3_ms_cd8_ENCSR278QHR" "h3k4me3_ms_cd8_ENCSR516CKJ" "h3k4me3_ms_cd8_ENCSR535YYH" "h3k4me3_ms_cd8_ENCSR741XAE" "h3k4me3_ms_cd8_ENCSR848XJL" "h3k9me3_ms_cd4_ENCSR057IZD" "h3k9me3_ms_cd4_ENCSR550DPT" "h3k9me3_ms_cd4_ENCSR677OEF" "h3k9me3_ms_cd4_ENCSR729ZXH" "h3k9me3_ms_cd4_ENCSR851RJV" "h3k9me3_ms_cd4_ENCSR919DFZ" "h3k9me3_ms_cd4_ENCSR953STZ" "h3k9me3_ms_cd4_ENCSR959VZU" "h3k9me3_ms_cd8_ENCSR101USF" "h3k9me3_ms_cd8_ENCSR354GNT" "h3k9me3_ms_cd8_ENCSR377QYB" "h3k9me3_ms_cd8_ENCSR733LCG" "h3k9me3_ms_cd8_ENCSR980KTW" "h3k27ac_ms_cd4_ENCSR200SSJ" "h3k27ac_ms_cd4_ENCSR322MTA" "h3k27ac_ms_cd4_ENCSR331WMS" "h3k27ac_ms_cd4_ENCSR350UKV" "h3k27ac_ms_cd4_ENCSR474PYR" "h3k27ac_ms_cd4_ENCSR520QDR" "h3k27ac_ms_cd4_ENCSR540XNK" "h3k27ac_ms_cd4_ENCSR705VSO" "h3k27ac_ms_cd4_ENCSR832UMM" "h3k27ac_ms_cd8_ENCSR078ATS" "h3k27ac_ms_cd8_ENCSR348YRH" "h3k27ac_ms_cd8_ENCSR458TOW" "h3k27ac_ms_cd8_ENCSR476IPR" "h3k27ac_ms_cd8_ENCSR787HDF" "h3k27ac_ms_cd8_ENCSR923JIB" "h3k27me3_ms_nk_ENCSR469QVG" "h3k27me3_ms_nk_ENCSR565WDW" "h3k36me3_ms_nk_ENCSR158VSE" "h3k36me3_ms_nk_ENCSR245KON" "h3k36me3_ms_nk_ENCSR530YDY" "h3k36me_ms_cd4_ENCSR532PXR" "h3k27me3_ms_cd4_ENCSR277XYX" "h3k27me3_ms_cd4_ENCSR526TNC" "h3k27me3_ms_cd4_ENCSR592EKF" "h3k27me3_ms_cd4_ENCSR613UFD" "h3k27me3_ms_cd4_ENCSR740SDR" "h3k27me3_ms_cd4_ENCSR779JLY" "h3k27me3_ms_cd4_ENCSR993CTA" "h3k27me3_ms_cd8_ENCSR116FVG" "h3k27me3_ms_cd8_ENCSR122JCM" "h3k27me3_ms_cd8_ENCSR216ZVA" "h3k27me3_ms_cd8_ENCSR284IKS" "h3k27me3_ms_cd8_ENCSR385BOZ" "h3k27me3_ms_cd8_ENCSR521SFR" "h3k36me3_ms_cd4_ENCSR276NGH" "h3k36me3_ms_cd4_ENCSR330CQU" "h3k36me3_ms_cd4_ENCSR482VIB" "h3k36me3_ms_cd4_ENCSR785PKM" "h3k36me3_ms_cd4_ENCSR865FTW" "h3k36me3_ms_cd4_ENCSR898VJE" "h3k36me3_ms_cd8_ENCSR239OWL" "h3k36me3_ms_cd8_ENCSR435MSB" "h3k36me3_ms_cd8_ENCSR631VIW" "h3k36me3_ms_cd8_ENCSR652WFO" "h3k36me3_ms_cd8_ENCSR757FGN" "h3k4me1_normal_b_ENCSR156BXM" "h3k4me3_normal_b_ENCSR791CAF" "h3k9me3_normal_b_ENCSR445LTM" "h3k27ac_normal_b_ENCSR685KZA" "h3k4me1_normal_nk_ENCSR277YKG" "h3k4me3_normal_nk_ENCSR394JFQ" "h3k9me3_normal_nk_ENCSR025UNZ" "h3k27ac_normal_nk_ENCSR977FMZ" "h3k27me3_normal_b_ENCSR589LHR" "h3k36me3_normal_b_ENCSR831AXK" "h3k4me1_normal_cd4_ENCSR102SOR" "h3k4me1_normal_cd8_ENCSR217SHH" "h3k4me3_normal_cd4_ENCSR537KJA" "h3k4me3_normal_cd8_ENCSR123ZAT" "h3k9me3_normal_cd4_ENCSR433EWI" "h3k9me3_normal_cd8_ENCSR294HTM" "h3k27ac_normal_cd4_ENCSR819NCZ" "h3k27ac_normal_cd8_ENCSR976RWL" "h3k27me3_normal_nk_ENCSR639NIG" "h3k36me3_normal_nk_ENCSR056GJY" "h3k27me3_normal_cd4_ENCSR068YVZ" "h3k27me3_normal_cd8_ENCSR720QRX" "h3k36me3_normal_cd4_ENCSR864OKB" "h3k36me3_normal_cd8_ENCSR303SQG" "h3k27me3_normal_cd4_ENCSR002UMT" "h3k4me3_normal_cd4_ENCSR935ELX")

# Process each sample folder
for sample_id in "${sample_ids[@]}"; do
    input_folder="/cta/users/merve.kaftancioglu/alternative_scenario/seqs/$sample_id"
    output_folder="/cta/users/merve.kaftancioglu/alternative_scenario/trimmomatic_as/output/$sample_id"

    # Create output folder if it doesn't exist
    mkdir -p "$output_folder"

    # Process all .fastq.gz files in the specified folder
    for input_file in "$input_folder"/*.fastq.gz; do
        output_file="$output_folder/output_paired_$(basename "$input_file")"
        trimmomatic SE -threads 16 -phred33 -trimlog "$output_folder/trimmomatic_log.txt" "$input_file" "$output_file" HEADCROP:14
    done
done

```
