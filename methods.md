## Bash Scripts
### FastQC (samples divided into 31 files): 

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

### MultiQC (without bash):

```
module load miniconda/22.9
eval "$(conda shell.bash hook)"
conda activate multiqc
multiqc
```
<br>

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
multiqc -o /cta/users/merve.kaftancioglu/alternative_scenario/multiqc_data/h3k4me3_ms_cd8_ENCSR516CKJ/cta/users/merve.kaftancioglu/alternative_scenario/fastQC_as/output/h3k4me3_ms_cd8_ENCSR516CKJ/*_fastqc.zip
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
<br>

### Trimmomatic (divided into ~50 files):

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
<br>

### CutAdapt (didn't work):

```
#!/bin/bash
#SBATCH --job-name=cut
#SBATCH --nodes=1
#SBATCH --account=user
#SBATCH --qos=mid_investor
#SBATCH --partition=mid_investor
#SBATCH --cpus-per-task=8
#SBATCH -o %x-%j-slurm.out
#SBATCH -e %x-%j-slurm.err

# Module
module load cutadapt/4.0

# Define input and output directories
input_dir="/cta/users/merve.kaftancioglu/alternative_scenario/seqs/$sample_id"
output_dir="/cta/users/merve.kaftancioglu/alternative_scenario/cutadapt_M1/output/$sample_id"
mkdir -p $output_dir

# Define sample IDs
sample_ids=("h3k27ac_ms_b_ENCSR295QZX" "h3k27ac_ms_b_ENCSR538URI" "h3k27ac_ms_b_ENCSR617GMR" "h3k4me1_ms_b_ENCSR084FXT" "h3k4me1_ms_b_ENCSR480ITK" "h3k4me1_ms_b_ENCSR759DHN" "h3k4me1_ms_b_ENCSR900IAM" "h3k4me1_ms_b_ENCSR931JYC" "h3k4me3_ms_b_ENCSR260CRI" "h3k4me3_ms_b_ENCSR461QMZ" "h3k4me3_ms_b_ENCSR530AVK" "h3k4me3_ms_b_ENCSR713ZYF" "h3k4me3_ms_b_ENCSR923EGO" "h3k9me3_ms_b_ENCSR445DNH" "h3k9me3_ms_b_ENCSR486SZN" "h3k9me3_ms_b_ENCSR983MKX" "h3k27ac_ms_b_ENCSR364OIK" "h3k27ac_ms_b_ENCSR541PLO" "h3k36me3_ms_cd4_ENCSR471WZE" "h3k36me3_ms_cd4_ENCSR473TGK" "h3k4me1_ms_nk_ENCSR482UGV" "h3k4me1_ms_nk_ENCSR718LCI" "h3k4me1_ms_nk_ENCSR806HUT" "h3k4me3_ms_nk_ENCSR234BAD" "h3k4me3_ms_nk_ENCSR703TWY" "h3k4me3_ms_nk_ENCSR755ZGS" "h3k9me3_ms_nk_ENCSR061ATV" "h3k9me3_ms_nk_ENCSR611YIA" "h3k9me3_ms_nk_ENCSR667HUI" "h3k27ac_ms_nk_ENCSR246TTM" "h3k27ac_ms_nk_ENCSR469CFL" "h3k27ac_ms_nk_ENCSR746AIX" "h3k27me3_ms_b_ENCSR009ZRH" "h3k27me3_ms_b_ENCSR182NLA" "h3k27me3_ms_b_ENCSR272YVX" "h3k27me3_ms_b_ENCSR649FUX" "h3k27me3_ms_b_ENCSR842WWX" "h3k36me3_ms_b_ENCSR089VQL" "h3k36me3_ms_b_ENCSR119PSR" "h3k36me3_ms_b_ENCSR238WFK" "h3k36me3_ms_b_ENCSR987OPY" "h3k4me1_ms_cd4_ENCSR036YVP" "h3k4me1_ms_cd4_ENCSR043JIA" "h3k4me1_ms_cd4_ENCSR093VUP" "h3k4me1_ms_cd4_ENCSR124IGC" "h3k4me1_ms_cd4_ENCSR315MYP" "h3k4me1_ms_cd4_ENCSR458QEY" "h3k4me1_ms_cd4_ENCSR485FBT" "h3k4me1_ms_cd4_ENCSR641ZFV" "h3k4me1_ms_cd4_ENCSR687CQX" "h3k4me1_ms_cd8_ENCSR231XAP" "h3k4me1_ms_cd8_ENCSR572XTB" "h3k4me1_ms_cd8_ENCSR788TEF" "h3k4me1_ms_cd8_ENCSR815YZL" "h3k4me1_ms_cd8_ENCSR861FEC" "h3k4me1_ms_cd8_ENCSR940PHE" "h3k4me3_ms_cd4_ENCSR180NCM" "h3k4me3_ms_cd4_ENCSR341QLC" "h3k4me3_ms_cd4_ENCSR482TGI" "h3k4me3_ms_cd4_ENCSR486XJK" "h3k4me3_ms_cd4_ENCSR496LKR" "h3k4me3_ms_cd4_ENCSR603LTN" "h3k4me3_ms_cd4_ENCSR802MXQ" "h3k4me3_ms_cd4_ENCSR878YHM" "h3k4me3_ms_cd4_ENCSR954ZLD" "h3k4me3_ms_cd8_ENCSR231ZZH" "h3k4me3_ms_cd8_ENCSR278QHR" "h3k4me3_ms_cd8_ENCSR516CKJ" "h3k4me3_ms_cd8_ENCSR535YYH" "h3k4me3_ms_cd8_ENCSR741XAE" "h3k4me3_ms_cd8_ENCSR848XJL" "h3k9me3_ms_cd4_ENCSR057IZD" "h3k9me3_ms_cd4_ENCSR550DPT" "h3k9me3_ms_cd4_ENCSR677OEF" "h3k9me3_ms_cd4_ENCSR729ZXH" "h3k9me3_ms_cd4_ENCSR851RJV" "h3k9me3_ms_cd4_ENCSR919DFZ" "h3k9me3_ms_cd4_ENCSR953STZ" "h3k9me3_ms_cd4_ENCSR959VZU" "h3k9me3_ms_cd8_ENCSR101USF" "h3k9me3_ms_cd8_ENCSR354GNT" "h3k9me3_ms_cd8_ENCSR377QYB" "h3k9me3_ms_cd8_ENCSR733LCG" "h3k9me3_ms_cd8_ENCSR980KTW" "h3k27ac_ms_cd4_ENCSR200SSJ" "h3k27ac_ms_cd4_ENCSR322MTA" "h3k27ac_ms_cd4_ENCSR331WMS" "h3k27ac_ms_cd4_ENCSR350UKV" "h3k27ac_ms_cd4_ENCSR474PYR" "h3k27ac_ms_cd4_ENCSR520QDR" "h3k27ac_ms_cd4_ENCSR540XNK" "h3k27ac_ms_cd4_ENCSR705VSO" "h3k27ac_ms_cd4_ENCSR832UMM" "h3k27ac_ms_cd8_ENCSR078ATS" "h3k27ac_ms_cd8_ENCSR348YRH" "h3k27ac_ms_cd8_ENCSR458TOW" "h3k27ac_ms_cd8_ENCSR476IPR" "h3k27ac_ms_cd8_ENCSR787HDF" "h3k27ac_ms_cd8_ENCSR923JIB" "h3k27me3_ms_nk_ENCSR469QVG" "h3k27me3_ms_nk_ENCSR565WDW" "h3k36me3_ms_nk_ENCSR158VSE" "h3k36me3_ms_nk_ENCSR245KON" "h3k36me3_ms_nk_ENCSR530YDY" "h3k36me_ms_cd4_ENCSR532PXR" "h3k27me3_ms_cd4_ENCSR277XYX" "h3k27me3_ms_cd4_ENCSR526TNC" "h3k27me3_ms_cd4_ENCSR592EKF" "h3k27me3_ms_cd4_ENCSR613UFD" "h3k27me3_ms_cd4_ENCSR740SDR" "h3k27me3_ms_cd4_ENCSR779JLY" "h3k27me3_ms_cd4_ENCSR993CTA" "h3k27me3_ms_cd8_ENCSR116FVG" "h3k27me3_ms_cd8_ENCSR122JCM" "h3k27me3_ms_cd8_ENCSR216ZVA" "h3k27me3_ms_cd8_ENCSR284IKS" "h3k27me3_ms_cd8_ENCSR385BOZ" "h3k27me3_ms_cd8_ENCSR521SFR" "h3k36me3_ms_cd4_ENCSR276NGH" "h3k36me3_ms_cd4_ENCSR330CQU" "h3k36me3_ms_cd4_ENCSR482VIB" "h3k36me3_ms_cd4_ENCSR785PKM" "h3k36me3_ms_cd4_ENCSR865FTW" "h3k36me3_ms_cd4_ENCSR898VJE" "h3k36me3_ms_cd8_ENCSR239OWL" "h3k36me3_ms_cd8_ENCSR435MSB" "h3k36me3_ms_cd8_ENCSR631VIW" "h3k36me3_ms_cd8_ENCSR652WFO" "h3k36me3_ms_cd8_ENCSR757FGN" "h3k4me1_normal_b_ENCSR156BXM" "h3k4me3_normal_b_ENCSR791CAF" "h3k9me3_normal_b_ENCSR445LTM" "h3k27ac_normal_b_ENCSR685KZA" "h3k4me1_normal_nk_ENCSR277YKG" "h3k4me3_normal_nk_ENCSR394JFQ" "h3k9me3_normal_nk_ENCSR025UNZ" "h3k27ac_normal_nk_ENCSR977FMZ" "h3k27me3_normal_b_ENCSR589LHR" "h3k36me3_normal_b_ENCSR831AXK" "h3k4me1_normal_cd4_ENCSR102SOR" "h3k4me1_normal_cd8_ENCSR217SHH" "h3k4me3_normal_cd4_ENCSR537KJA" "h3k4me3_normal_cd8_ENCSR123ZAT" "h3k9me3_normal_cd4_ENCSR433EWI" "h3k9me3_normal_cd8_ENCSR294HTM" "h3k27ac_normal_cd4_ENCSR819NCZ" "h3k27ac_normal_cd8_ENCSR976RWL" "h3k27me3_normal_nk_ENCSR639NIG" "h3k36me3_normal_nk_ENCSR056GJY" "h3k27me3_normal_cd4_ENCSR068YVZ" "h3k27me3_normal_cd8_ENCSR720QRX" "h3k36me3_normal_cd4_ENCSR864OKB" "h3k36me3_normal_cd8_ENCSR303SQG" "h3k27me3_normal_cd4_ENCSR002UMT" "h3k4me3_normal_cd4_ENCSR935ELX")

for file in "${files[@]}"; do
    input_file="$input_dir/${file}.fastq.gz"
    output_file="$output_dir/${file}_processed.fastq.gz"
    cutadapt -b file:adapters.fasta -j 8 -q 15 -o -p "$output_file" "$input_file"
done

```
<br>

### STAR (didn't work):
(a)

```
#!/bin/bash
#SBATCH --job-name=cut
#SBATCH --nodes=1
#SBATCH --account=user
#SBATCH --qos=mid_investor
#SBATCH --partition=mid_investor
#SBATCH --cpus-per-task=16
#SBATCH -o /cta/users/merve.kaftancioglu/alternative_scenario/star_as/out-err/%x-%j-slurm.out
#SBATCH -e /cta/users/merve.kaftancioglu/alternative_scenario/star_as/out-err/%x-%j-slurm.err

# Load STAR module
module load star-2.7.0e-gcc-9.2.0-vynasg3

# Define paths (replace with your actual paths)
genome_fasta="/cta/users/merve.kaftancioglu/alternative_scenario/star_as/reference/human_g1k_v37.fasta"
star_index_dir="/cta/users/merve.kaftancioglu/star_as/index"

# Check if genome fasta file exists
if [ ! -f "$genome_fasta" ]; then
  echo "Error: Genome fasta file not found at: $genome_fasta"
  exit 1
fi

# Check if star_index_dir exists (and create if not)
if [ ! -d "$star_index_dir" ]; then
  echo "Creating star_genome_index directory: $star_index_dir"
  mkdir -p "$star_index_dir"
fi

# Define number of threads (adjust based on your system)
threads=8

# STAR command
star_cmd="STAR \
  --runThreadN $threads \
  --runMode genomeGenerate \
  --genomeDir $star_index_dir \
  --genomeFastaFiles $genome_fasta"

# Run STAR and handle potential errors
echo "Building STAR genome index..."
$star_cmd

if [ $? -ne 0 ]; then
  echo "Error: STAR failed to generate genome index."
  exit 1
fi

echo "STAR genome index created successfully in: $star_index_dir"
```

<br>

(b)

```
#!/bin/bash
#SBATCH --job-name=star_single
#SBATCH --nodes=1
#SBATCH --account=users
#SBATCH --qos=mid_mdbf
#SBATCH --partition=mid_mdbf
#SBATCH --cpus-per-task=16
#SBATCH -o /cta/users/merve.kaftancioglu/alternative_scenario/star_as/out-err/%x-%j-slurm.out
#SBATCH -e /cta/users/merve.kaftancioglu/alternative_scenario/star_as/out-err/%x-%j-slurm.err

# Define paths 
genome_fasta="/cta/users/merve.kaftancioglu/alternative_scenario/star_as/reference/human_g1k_v37.fasta"
star_index_dir="/cta/users/merve.kaftancioglu/alternative_scenario/star_as/index"

# Check if genome fasta file exists
if [ ! -f "$genome_fasta" ]; then
  echo "Error: Genome fasta file not found at: $genome_fasta"
  exit 1
fi

# Check if star_index_dir exists (and create if not)
if [ ! -d "$star_index_dir" ]; then
  echo "Creating star_genome_index directory: $star_index_dir"
  mkdir -p "$star_index_dir"
fi

# Define number of threads (adjust based on your system)
threads=8

# STAR command
star_cmd="STAR \
  --runThreadN $threads \
  --runMode genomeGenerate \
  --genomeDir $star_index_dir \
  --genomeFastaFiles $genome_fasta"

# Run STAR and handle potential errors
echo "Building STAR genome index..."
$star_cmd

if [ $? -ne 0 ]; then
  echo "Error: STAR failed to generate genome index."
  exit 1
fi

echo "STAR genome index created successfully in: $star_index_dir"

```

<br>

(c)

```
#!/bin/bash
#SBATCH --job-name=star_alignment
#SBATCH --nodes=1
#SBATCH --account=users
#SBATCH --qos=mid_mdbf
#SBATCH --partition=mid_mdbf
#SBATCH --cpus-per-task=16
#SBATCH -o /cta/users/merve.kaftancioglu/alternative_scenario/bowtie2_as/out-err/%x-%j-slurm.out
#SBATCH -e /cta/users/merve.kaftancioglu/alternative_scenario/bowtie2_as/out-err/%x-%j-slurm.err

# Load STAR module
module load star-2.7.0e-gcc-9.2.0-vynasg3

# Define directories
genome_dir="/cta/users/merve.kaftancioglu/alternative_scenario/star_as/index"  # Update with your genome directory
input_dir="/cta/users/merve.kaftancioglu/alternative_scenario/trimmomatic_as/output"  # Update with your input directory
output_dir="/cta/users/merve.kaftancioglu/alternative_scenario/star_as/output"  # Update with your output directory

# Run STAR for each folder within the input directory
for sample_dir in "$input_dir"/*; do
  # Check if it's a directory (avoid hidden files or non-directories)
  if [[ -d "$sample_dir" ]]; then
    # Extract sample ID from directory name (modify if needed)
    sample_id=$(basename "$sample_dir")

    # Define input and output paths with sample ID
    input_fastq="$sample_dir"/*.fastq.gz
    output_bam="$output_dir/$sample_id/$sample_id.bam"

    # Run STAR for the current sample
    STAR --genomeDir "$genome_dir" \
         --runThreadN 16 \  # Adjust threads based on your resource allocation
         --readFilesIn "$input_fastq" \  # No colon after flag name
         --outFileNamePrefix "$output_dir/$sample_id/" \  # Include trailing slash for proper path
         --outSAMtype BAM SortedByCoordinate \  # No colon after flag name
         --outSAMunmapped Within \
         --outSAMattributes Standard \
         --alignEndsType EndToEnd
  fi
done

```

<br>

### Bowtie2 (didn't work):

(a)

```
#!/bin/bash
#SBATCH --job-name=bowtie2_longread_alignment
#SBATCH --nodes=1
#SBATCH --account=users
#SBATCH --qos=mid_mdbf
#SBATCH --partition=mid_mdbf
#SBATCH --cpus-per-task=16
#SBATCH -o /cta/users/merve.kaftancioglu/alternative_scenario/bowtie2_as/out-err/%x-%j-slurm.out
#SBATCH -e /cta/users/merve.kaftancioglu/alternative_scenario/bowtie2_as/out-err/%x-%j-slurm.err

# Load Bowtie2 module
module load bowtie2-2.3.5.1-gcc-9.2.0-ae4dozk

# Define directories (consider adjusting paths)
input_dir="/cta/users/merve.kaftancioglu/alternative_scenario/trimmomatic_as/output"
output_dir="/cta/users/merve.kaftancioglu/alternative_scenario/bowtie2_as/longread_output"  # Separate output directory for long reads

# Path to reference genome FASTA file
reference_genome="/cta/users/merve.kaftancioglu/alternative_scenario/bowtie2_as/reference/human_g1k_v37.fasta"

# Create Bowtie2 index (if not already done)
bowtie2-build "$reference_genome" /cta/users/merve.kaftancioglu/alternative_scenario/bowtie2_as/index

# Loop through all FASTQ files in the input directory
for input_fastq in $(find "${input_dir}" -type f -name "*.fastq.gz" ! -name "*.*"); do

  sample_id="${input_fastq##*/}"  
  sample_id="${sample_id//.fastq.gz}" 

  echo "Processing ${sample_id}..."

  # Create output directory (if it doesn't exist)
  mkdir -p "${output_dir}/${sample_id}"

  # Output SAM file path
  output_sam="${output_dir}/${sample_id}/${sample_id}.sam"

  # Run Bowtie2 with the current FASTQ file, specifying output format (-S) as SAM and "--very-sensitive" flag
  bowtie2 --threads 16 -p 16 --trim3 30 --very-sensitive -x /cta/users/merve.kaftancioglu/alternative_scenario/bowtie2_as/index -U $input_fastq -S $output_sam

  # Optional: Check for successful completion (add this block if desired)
  if [ $? -eq 0 ]; then
    echo "Alignment for ${sample_id} completed successfully."
  else
    echo "Error occurred during alignment for ${sample_id}. Check slurm error logs for details."
  fi
done

```

<br>

(b)

```
#!/bin/bash
#SBATCH --job-name=bowtie2_alignment
#SBATCH --nodes=1
#SBATCH --account=users
#SBATCH --qos=mid_mdbf
#SBATCH --partition=mid_mdbf
#SBATCH --cpus-per-task=16
#SBATCH -o /cta/users/merve.kaftancioglu/alternative_scenario/bowtie2_as/out-err/%x-%j-slurm.out
#SBATCH -e /cta/users/merve.kaftancioglu/alternative_scenario/bowtie2_as/out-err/%x-%j-slurm.err

# Load Bowtie2 module
module load bowtie2-2.3.5.1-gcc-9.2.0-ae4dozk

# Define directories
input_dir="/cta/users/merve.kaftancioglu/alternative_scenario/trimmomatic_as/output"
output_dir="/cta/users/merve.kaftancioglu/alternative_scenario/bowtie2_as/output"

# Path to reference genome FASTA file
reference_genome="/cta/users/merve.kaftancioglu/alternative_scenario/bowtie2_as/reference/human_g1k_v37.fasta"

# Create Bowtie2 index (if not already done)
bowtie2-build "$reference_genome" /cta/users/merve.kaftancioglu/alternative_scenario/bowtie2_as/index

# Loop through all FASTQ files in the input directory using find
for input_fastq in $(find "${input_dir}" -type f -name "*.fastq.gz" ! -name "*.*"); do

  sample_id="${input_fastq##*/}"  
  sample_id="${sample_id//.fastq.gz}" 

  echo "Processing ${sample_id}..."

  # Create output directory (if it doesn't exist)
  mkdir -p "${output_dir}/${sample_id}"

  # Output SAM file path
  output_sam="${output_dir}/${sample_id}/${sample_id}.sam"

  # Run Bowtie2 with the current FASTQ file, specifying output format (-S) as SAM
  bowtie2 --threads 16 -p 16 --trim3 30 -x /cta/users/merve.kaftancioglu/alternative_scenario/bowtie2_as/index -U $input_fastq -S $output_sam
done

```
<br>

(c)

```
#!/bin/bash
#SBATCH --job-name=bowtie2_alignment
#SBATCH --nodes=1
#SBATCH --account=users
#SBATCH --qos=mid_mdbf
#SBATCH --partition=mid_mdbf
#SBATCH --cpus-per-task=16
#SBATCH -o /cta/users/merve.kaftancioglu/alternative_scenario/bowtie2_as/out-err/%x-%j-slurm.out
#SBATCH -e /cta/users/merve.kaftancioglu/alternative_scenario/bowtie2_as/out-err/%x-%j-slurm.err

# Load Bowtie2 module
module load bowtie2-2.3.5.1-gcc-9.2.0-ae4dozk

# Define directories
input_dir="/cta/users/merve.kaftancioglu/alternative_scenario/trimmomatic_as/output"
output_dir="/cta/users/merve.kaftancioglu/alternative_scenario/bowtie2_as/output"

# Path to reference genome FASTA file
reference_genome="/cta/users/merve.kaftancioglu/alternative_scenario/bowtie2_as/reference/human_g1k_v37.fasta"

# Create Bowtie2 index (if not already done)
bowtie2-build "$reference_genome" /cta/users/merve.kaftancioglu/alternative_scenario/bowtie2_as/index

# Loop through all FASTQ files in the input directory using find
for input_fastq in $(find "${input_dir}" -type f -name "*.fastq.gz" ! -name "*.*"); do

  sample_id="${input_fastq##*/}"  
  sample_id="${sample_id//.fastq.gz}" 

  echo "Processing ${sample_id}..."

  # Create output directory (if it doesn't exist)
  mkdir -p "${output_dir}/${sample_id}"

  # Output BAM file path
  output_bam="${output_dir}/${sample_id}/${sample_id}.bam"

  # Run Bowtie2 with the current FASTQ file
  bowtie2 --threads 16 -p 16 --trim3 30 -x reference_genome -U $input_fastq -S $output_bam
  # bowtie2 -x /cta/users/merve.kaftancioglu/alternative_scenario/bowtie2_as/index -U "$input_fastq" -S "$output_bam"
done

```

<br>

(d)

```
#!/bin/bash
#SBATCH --job-name=bowtie2_alignment
#SBATCH --nodes=1
#SBATCH --account=users
#SBATCH --qos=mid_mdbf
#SBATCH --partition=mid_mdbf
#SBATCH --cpus-per-task=16
#SBATCH -o /cta/users/merve.kaftancioglu/alternative_scenario/bowtie2_as/out-err/%x-%j-slurm.out
#SBATCH -e /cta/users/merve.kaftancioglu/alternative_scenario/bowtie2_as/out-err/%x-%j-slurm.err

# Load Bowtie2 module
module load bowtie2-2.3.5.1-gcc-9.2.0-ae4dozk

# Define directories
input_dir="/cta/users/merve.kaftancioglu/alternative_scenario/trimmomatic_as/output"
output_dir="/cta/users/merve.kaftancioglu/alternative_scenario/bowtie2_as/output"

# Path to reference genome FASTA file
reference_genome="/cta/users/merve.kaftancioglu/alternative_scenario/bowtie2_as/reference/human_g1k_v37.fasta"

# Create Bowtie2 index (if not already done)
bowtie2-build "$reference_genome" /cta/users/merve.kaftancioglu/alternative_scenario/bowtie2_as/index

# Loop through all FASTQ files in the input directory
for input_fastq in "${input_dir}"/*.fastq.gz; do
  # Extract sample ID from filename (assuming filenames start with sample ID)
  # Modify this part to match your specific filename format if needed
  sample_id="${input_fastq##*/}"  # Double ## removes everything before the last slash (/)
  sample_id="${sample_id//.fastq.gz}"  # Remove .fastq.gz extension

  echo "Processing ${sample_id}..."

  # Create output directory (if it doesn't exist)
  mkdir -p "${output_dir}/${sample_id}"

  # Output BAM file path
  output_bam="${output_dir}/${sample_id}/${sample_id}.bam"

  # Run Bowtie2 with the current FASTQ file
  bowtie2 -x /cta/users/merve.kaftancioglu/alternative_scenario/bowtie2_as/index -U "$input_fastq" -S "$output_bam"
done

```

<br>

(e)

```
#!/bin/bash
#SBATCH --job-name=bowtie2_alignment
#SBATCH --nodes=1
#SBATCH --account=users
#SBATCH --qos=mid_mdbf
#SBATCH --partition=mid_mdbf
#SBATCH --cpus-per-task=16
#SBATCH -o /cta/users/merve.kaftancioglu/alternative_scenario/bowtie2_as/out-err/%x-%j-slurm.out
#SBATCH -e /cta/users/merve.kaftancioglu/alternative_scenario/bowtie2_as/out-err/%x-%j-slurm.err

# Load Bowtie2 module
module load bowtie2-2.3.5.1-gcc-9.2.0-ae4dozk

# Define directories
input_dir="/cta/users/merve.kaftancioglu/alternative_scenario/trimmomatic_as/output"
output_dir="/cta/users/merve.kaftancioglu/alternative_scenario/bowtie2_as/output"

# Path to reference genome FASTA file
reference_genome="/cta/users/merve.kaftancioglu/alternative_scenario/bowtie2_as/reference/human_g1k_v37.fasta"

# Create Bowtie2 index (if not already done)
bowtie2-build "$reference_genome" /cta/users/merve.kaftancioglu/alternative_scenario/bowtie2_as/index

# Loop through all FASTQ files in the input directory using find
for input_fastq in $(find "${input_dir}" -type f -name "*.fastq.gz" ! -name "*.*"); do

  sample_id="${input_fastq##*/}"  
  sample_id="${sample_id//.fastq.gz}" 

  echo "Processing ${sample_id}..."

  # Create output directory (if it doesn't exist)
  mkdir -p "${output_dir}/${sample_id}"

  # Output SAM file path
  output_sam="${output_dir}/${sample_id}/${sample_id}.sam"

  # Run Bowtie2 with the current FASTQ file, specifying output format (-U) as SAM
  bowtie2 --threads 16 -p 16 --trim3 30 -x /cta/users/merve.kaftancioglu/alternative_scenario/bowtie2_as/index -U $input_fastq -U $output_sam
  # Note the change: -U $output_sam specifies the output SAM file
done

```

<br>

(f)

```
#!/bin/bash
#SBATCH --job-name=bowtie2_alignment
#SBATCH --nodes=1
#SBATCH --account=users
#SBATCH --qos=mid_mdbf
#SBATCH --partition=mid_mdbf
#SBATCH --cpus-per-task=16
#SBATCH -o /cta/users/merve.kaftancioglu/alternative_scenario/bowtie2_as/out-err/%x-%j-slurm.out
#SBATCH -e /cta/users/merve.kaftancioglu/alternative_scenario/bowtie2_as/out-err/%x-%j-slurm.err

# Load Bowtie2 module
module load bowtie2-2.3.5.1-gcc-9.2.0-ae4dozk

# Sample IDs 
sample_ids=("h3k27ac_ms_b_ENCSR295QZX" "h3k27ac_ms_b_ENCSR364OIK" "h3k27ac_ms_b_ENCSR538URI")

# Directories
input_dir="/cta/users/merve.kaftancioglu/alternative_scenario/trimmomatic_as/output"
output_dir="/cta/users/merve.kaftancioglu/alternative_scenario/bowtie2_as/output"
index_dir="/cta/users/merve.kaftancioglu/alternative_scenario/bowtie2_as/index"  

# Create Bowtie2 index 
bowtie2-build "$reference_genome" "$index_dir"

# Loop through samples
for sample_id in "${sample_ids[@]}"; do
  echo "Processing ${sample_id}..."

  # Create output directory 
  mkdir -p "${output_dir}/${sample_id}"

  # Combine input and output paths with sample ID
  input_fastq="${input_dir}/${sample_id}"/*.fastq.gz
  output_bam="${output_dir}/${sample_id}/${sample_id}.bam"

  # Run Bowtie2 with trimming
  bowtie2 --threads 16 -p 16 --trim3 30 -x "$index_dir" -U "$input_fastq" | samtools view -bS - > "$output_bam"
  
  # Convert BAM to SAM
  input_bam="$output_bam"
  output_sam="${output_dir}/${sample_id}/${sample_id}.sam"
  samtools view -h -o "$output_sam" "$input_b

  echo "Finished processing ${sample_id}"
done

echo "All samples processed!"

```

<br>

### BWA (divided into 6 files):


```

#!/bin/bash
#SBATCH --job-name=bwa_single
#SBATCH --nodes=1
#SBATCH --account=users
#SBATCH --qos=cuda
#SBATCH --partition=cuda
#SBATCH --cpus-per-task=16
#SBATCH -o /cta/users/merve.kaftancioglu/alternative_scenario/bwa_as/out-err/%x-%j-slurm.out
#SBATCH -e /cta/users/merve.kaftancioglu/alternative_scenario/bwa_as/out-err/%x-%j-slurm.err

# Load BWA module
module load bwa-0.7.17-gcc-9.2.0-xwkclqy

# Define directories (consider adjusting paths)
input_dir="/cta/users/merve.kaftancioglu/alternative_scenario/trimmomatic_as/output"
output_dir="/cta/users/merve.kaftancioglu/alternative_scenario/bwa_as/output"
reference_genome="/cta/users/merve.kaftancioglu/alternative_scenario/bwa_as/reference/human_g1k_v37.fasta"

# Check if BWA index exists
if [ ! -f "$reference_genome.bwt" ]; then
  echo "Error: BWA index not found. Please create the index using: bwa mem index $reference_genome"
  exit 1
fi

# Define sample IDs (modify as needed)
sample_ids=("h3k27ac_ms_b_ENCSR295QZX" "h3k27ac_ms_b_ENCSR538URI" "h3k27ac_ms_b_ENCSR617GMR" "h3k4me1_ms_b_ENCSR084FXT" "h3k4me1_ms_b_ENCSR480ITK" "h3k4me1_ms_b_ENCSR759DHN" "h3k4me1_ms_b_ENCSR900IAM" "h3k4me1_ms_b_ENCSR931JYC" "h3k4me3_ms_b_ENCSR260CRI" "h3k4me3_ms_b_ENCSR461QMZ" "h3k4me3_ms_b_ENCSR530AVK" "h3k4me3_ms_b_ENCSR713ZYF" "h3k4me3_ms_b_ENCSR923EGO" "h3k9me3_ms_b_ENCSR445DNH" "h3k9me3_ms_b_ENCSR486SZN" "h3k9me3_ms_b_ENCSR983MKX" "h3k27ac_ms_b_ENCSR364OIK" "h3k27ac_ms_b_ENCSR541PLO" "h3k36me3_ms_cd4_ENCSR471WZE" "h3k36me3_ms_cd4_ENCSR473TGK" "h3k4me1_ms_nk_ENCSR482UGV" "h3k4me1_ms_nk_ENCSR718LCI" "h3k4me1_ms_nk_ENCSR806HUT" "h3k4me3_ms_nk_ENCSR234BAD" "h3k4me3_ms_nk_ENCSR703TWY" "h3k4me3_ms_nk_ENCSR755ZGS" "h3k9me3_ms_nk_ENCSR061ATV" "h3k9me3_ms_nk_ENCSR611YIA" "h3k9me3_ms_nk_ENCSR667HUI" "h3k27ac_ms_nk_ENCSR246TTM" "h3k27ac_ms_nk_ENCSR469CFL" "h3k27ac_ms_nk_ENCSR746AIX" "h3k27me3_ms_b_ENCSR009ZRH" "h3k27me3_ms_b_ENCSR182NLA" "h3k27me3_ms_b_ENCSR272YVX" "h3k27me3_ms_b_ENCSR649FUX" "h3k27me3_ms_b_ENCSR842WWX" "h3k36me3_ms_b_ENCSR089VQL" "h3k36me3_ms_b_ENCSR119PSR" "h3k36me3_ms_b_ENCSR238WFK" "h3k36me3_ms_b_ENCSR987OPY" "h3k4me1_ms_cd4_ENCSR036YVP" "h3k4me1_ms_cd4_ENCSR043JIA" "h3k4me1_ms_cd4_ENCSR093VUP" "h3k4me1_ms_cd4_ENCSR124IGC" "h3k4me1_ms_cd4_ENCSR315MYP" "h3k4me1_ms_cd4_ENCSR458QEY" "h3k4me1_ms_cd4_ENCSR485FBT" "h3k4me1_ms_cd4_ENCSR641ZFV" "h3k4me1_ms_cd4_ENCSR687CQX" "h3k4me1_ms_cd8_ENCSR231XAP" "h3k4me1_ms_cd8_ENCSR572XTB" "h3k4me1_ms_cd8_ENCSR788TEF" "h3k4me1_ms_cd8_ENCSR815YZL" "h3k4me1_ms_cd8_ENCSR861FEC" "h3k4me1_ms_cd8_ENCSR940PHE" "h3k4me3_ms_cd4_ENCSR180NCM" "h3k4me3_ms_cd4_ENCSR341QLC" "h3k4me3_ms_cd4_ENCSR482TGI" "h3k4me3_ms_cd4_ENCSR486XJK" "h3k4me3_ms_cd4_ENCSR496LKR" "h3k4me3_ms_cd4_ENCSR603LTN" "h3k4me3_ms_cd4_ENCSR802MXQ" "h3k4me3_ms_cd4_ENCSR878YHM" "h3k4me3_ms_cd4_ENCSR954ZLD" "h3k4me3_ms_cd8_ENCSR231ZZH" "h3k4me3_ms_cd8_ENCSR278QHR" "h3k4me3_ms_cd8_ENCSR516CKJ" "h3k4me3_ms_cd8_ENCSR535YYH" "h3k4me3_ms_cd8_ENCSR741XAE" "h3k4me3_ms_cd8_ENCSR848XJL" "h3k9me3_ms_cd4_ENCSR057IZD" "h3k9me3_ms_cd4_ENCSR550DPT" "h3k9me3_ms_cd4_ENCSR677OEF" "h3k9me3_ms_cd4_ENCSR729ZXH" "h3k9me3_ms_cd4_ENCSR851RJV" "h3k9me3_ms_cd4_ENCSR919DFZ" "h3k9me3_ms_cd4_ENCSR953STZ" "h3k9me3_ms_cd4_ENCSR959VZU" "h3k9me3_ms_cd8_ENCSR101USF" "h3k9me3_ms_cd8_ENCSR354GNT" "h3k9me3_ms_cd8_ENCSR377QYB" "h3k9me3_ms_cd8_ENCSR733LCG" "h3k9me3_ms_cd8_ENCSR980KTW" "h3k27ac_ms_cd4_ENCSR200SSJ" "h3k27ac_ms_cd4_ENCSR322MTA" "h3k27ac_ms_cd4_ENCSR331WMS" "h3k27ac_ms_cd4_ENCSR350UKV" "h3k27ac_ms_cd4_ENCSR474PYR" "h3k27ac_ms_cd4_ENCSR520QDR" "h3k27ac_ms_cd4_ENCSR540XNK" "h3k27ac_ms_cd4_ENCSR705VSO" "h3k27ac_ms_cd4_ENCSR832UMM" "h3k27ac_ms_cd8_ENCSR078ATS" "h3k27ac_ms_cd8_ENCSR348YRH" "h3k27ac_ms_cd8_ENCSR458TOW" "h3k27ac_ms_cd8_ENCSR476IPR" "h3k27ac_ms_cd8_ENCSR787HDF" "h3k27ac_ms_cd8_ENCSR923JIB" "h3k27me3_ms_nk_ENCSR469QVG" "h3k27me3_ms_nk_ENCSR565WDW" "h3k36me3_ms_nk_ENCSR158VSE" "h3k36me3_ms_nk_ENCSR245KON" "h3k36me3_ms_nk_ENCSR530YDY" "h3k36me_ms_cd4_ENCSR532PXR" "h3k27me3_ms_cd4_ENCSR277XYX" "h3k27me3_ms_cd4_ENCSR526TNC" "h3k27me3_ms_cd4_ENCSR592EKF" "h3k27me3_ms_cd4_ENCSR613UFD" "h3k27me3_ms_cd4_ENCSR740SDR" "h3k27me3_ms_cd4_ENCSR779JLY" "h3k27me3_ms_cd4_ENCSR993CTA" "h3k27me3_ms_cd8_ENCSR116FVG" "h3k27me3_ms_cd8_ENCSR122JCM" "h3k27me3_ms_cd8_ENCSR216ZVA" "h3k27me3_ms_cd8_ENCSR284IKS" "h3k27me3_ms_cd8_ENCSR385BOZ" "h3k27me3_ms_cd8_ENCSR521SFR" "h3k36me3_ms_cd4_ENCSR276NGH" "h3k36me3_ms_cd4_ENCSR330CQU" "h3k36me3_ms_cd4_ENCSR482VIB" "h3k36me3_ms_cd4_ENCSR785PKM" "h3k36me3_ms_cd4_ENCSR865FTW" "h3k36me3_ms_cd4_ENCSR898VJE" "h3k36me3_ms_cd8_ENCSR239OWL" "h3k36me3_ms_cd8_ENCSR435MSB" "h3k36me3_ms_cd8_ENCSR631VIW" "h3k36me3_ms_cd8_ENCSR652WFO" "h3k36me3_ms_cd8_ENCSR757FGN" "h3k4me1_normal_b_ENCSR156BXM" "h3k4me3_normal_b_ENCSR791CAF" "h3k9me3_normal_b_ENCSR445LTM" "h3k27ac_normal_b_ENCSR685KZA" "h3k4me1_normal_nk_ENCSR277YKG" "h3k4me3_normal_nk_ENCSR394JFQ" "h3k9me3_normal_nk_ENCSR025UNZ" "h3k27ac_normal_nk_ENCSR977FMZ" "h3k27me3_normal_b_ENCSR589LHR" "h3k36me3_normal_b_ENCSR831AXK" "h3k4me1_normal_cd4_ENCSR102SOR" "h3k4me1_normal_cd8_ENCSR217SHH" "h3k4me3_normal_cd4_ENCSR537KJA" "h3k4me3_normal_cd8_ENCSR123ZAT" "h3k9me3_normal_cd4_ENCSR433EWI" "h3k9me3_normal_cd8_ENCSR294HTM" "h3k27ac_normal_cd4_ENCSR819NCZ" "h3k27ac_normal_cd8_ENCSR976RWL" "h3k27me3_normal_nk_ENCSR639NIG" "h3k36me3_normal_nk_ENCSR056GJY" "h3k27me3_normal_cd4_ENCSR068YVZ" "h3k27me3_normal_cd8_ENCSR720QRX" "h3k36me3_normal_cd4_ENCSR864OKB" "h3k36me3_normal_cd8_ENCSR303SQG" "h3k27me3_normal_cd4_ENCSR002UMT" "h3k4me3_normal_cd4_ENCSR935ELX")

# Loop through all sample subdirectories
for sample_dir in "$input_dir"/*; do
  # Check if directory name matches any sample ID
  if [[ " ${sample_ids[@]} " =~ " $(basename "$sample_dir") " ]]; then
    echo "Processing sample: $(basename "$sample_dir")"

    # Create output directory for the sample (if it doesn't exist)
    mkdir -p "$output_dir/$(basename "$sample_dir")"

    # Find FASTQ files within the sample directory
    # Adjust the pattern below if your file names differ
    fastq_files=( "$sample_dir"/*.fastq.gz )

    # Check if any FASTQ files found
    if [[ ${#fastq_files[@]} -eq 0 ]]; then
      echo "Warning: No FASTQ files found in $sample_dir"
      continue # Skip to next iteration if no FASTQ files
    fi

    # BWA mem command with paired-end or single-end handling
    if [[ ${#fastq_files[@]} -eq 2 ]]; then
      # Paired-end reads (assuming *_1.fastq.gz and *_2.fastq.gz naming)
      bwa mem $reference_genome ${fastq_files[0]} ${fastq_files[1]} > "$output_dir/$(basename "$sample_dir").sam"
    else
      # Single-end read (assuming only one FASTQ file)
      # Extract only the sample name from the directory path
      sample_name=$(basename "$sample_dir" /)
      bwa mem $reference_genome ${fastq_files[0]} > "$output_dir/$sample_name.sam"
    fi
  fi
done

```

<br>

### nf-core ChIP Seq (didn't work):

(a)

```

#!/bin/bash
#SBATCH --job-name=nf-core
#SBATCH --nodes=1
#SBATCH --account=users
#SBATCH --qos=cuda
#SBATCH --partition=cuda
#SBATCH --cpus-per-task=16
#SBATCH -o /cta/users/merve.kaftancioglu/alternative_scenario/processed_seqs/out-err/%x-%j-slurm.out
#SBATCH -e /cta/users/merve.kaftancioglu/alternative_scenario/processed_seqs/out-err/%x-%j-slurm.err

# Load necessary modules (if applicable for your cluster)
module load docker

# Specify the Docker image (replace with nf-core/chipseq version you downloaded)
docker_image=nf-core/chipseq:latest

# Navigate to the directory containing your Slurm script and pipeline files
cd /cta/users/merve.kaftancioglu/alternative_scenario/processed_seqs/all_data

# Create symbolic links for all fastq files in the current directory
for file in *.fastq; do ln -s "$file" all_fastq_files; done

# Run the nextflow command within the container
docker run --rm -v /cta/users/merve.kaftancioglu/alternative_scenario/trimmomatic_as/output

# Create a directory to hold symbolic links (soft links) to your fastq files (optional)
# This step is optional, you can comment it out if you don't want to create a separate directory
# mkdir all_fastq_files

# Loop through all subdirectories under processed_seqs
for sample_dir in */ ; do
  # Navigate to the current subdirectory
  cd "$sample_dir"

  # Loop through all fastq.gz files in the current subdirectory
  for file in *.fastq.gz; do
    # Create a symbolic link (optional, uncomment if using the directory)
    # ln -s "$file" ../all_fastq_files  # Uncomment if using all_fastq_files directory
    
    # Construct the full path to the fastq file (alternative to symbolic links)
    full_fastq_path=$(pwd)/"$file"
    
    # You can echo the full path for verification (optional)
    # echo "Found fastq: $full_fastq_path"
  done
  
```

<br>

(b)

```
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
# SET /P blacklist_path=Enter the path to your blacklist BED file (optional): 
SET /P output_path=/cta/users/merve.kaftancioglu/alternative_scenario/nf-core/chipseq-2.0.0/output: 
SET /P profile=singularity:3.8.4  

REM Build the command with user inputs
SET command=nextflow run nf-core/chipseq -r 2.0.0 

#REM Add options based on data format (only for paired-end)
#IF "%data_format%"=="paired-end" (
#  SET command=%command% --paired
)

SET command=%command% --input %input_path% --genome %genome_path% --outdir %output_path% -profile %profile%

#REM Add optional blacklist path (uncomment if needed)
#IF NOT "%blacklist_path%"=="" (
#  SET command=%command% --blacklist %blacklist_path%
#)

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

```

<br>

(c)

```
#!/bin/bash
#SBATCH --job-name=bwa_single
#SBATCH --nodes=1
#SBATCH --account=users
#SBATCH --qos=cuda
#SBATCH --partition=cuda
#SBATCH --cpus-per-task=16
#SBATCH -o /cta/users/merve.kaftancioglu/alternative_scenario/processed_seqs/out-err/%x-%j-slurm.out
#SBATCH -e /cta/users/merve.kaftancioglu/alternative_scenario/processed_seqs/out-err/%x-%j-slurm.err

# Load BWA module
module load rsync-3.1.3-gcc-9.2.0-mrwcim2

# Define directories (consider adjusting paths)
input_dir="/cta/users/merve.kaftancioglu/alternative_scenario/processed_seqs/all_data"
output_dir="/cta/users/merve.kaftancioglu/alternative_scenario/bwa_as/output"
reference_genome="/cta/users/merve.kaftancioglu/alternative_scenario/bwa_as/reference/human_g1k_v37.fasta"

# Check if BWA index exists
if [ ! -f "$reference_genome.bwt" ]; then
  echo "Error: BWA index not found. Please create the index using: bwa mem index $reference_genome"
  exit 1
fi

# Define sample IDs (modify as needed)
sample_ids=("h3k27ac_ms_b_ENCSR295QZX" "h3k27ac_ms_b_ENCSR538URI" "h3k27ac_ms_b_ENCSR617GMR" "h3k4me1_ms_b_ENCSR084FXT" "h3k4me1_ms_b_ENCSR480ITK" "h3k4me1_ms_b_ENCSR759DHN" "h3k4me1_ms_b_ENCSR900IAM" "h3k4me1_ms_b_ENCSR931JYC" "h3k4me3_ms_b_ENCSR260CRI" "h3k4me3_ms_b_ENCSR461QMZ" "h3k4me3_ms_b_ENCSR530AVK" "h3k4me3_ms_b_ENCSR713ZYF" "h3k4me3_ms_b_ENCSR923EGO" "h3k9me3_ms_b_ENCSR445DNH" "h3k9me3_ms_b_ENCSR486SZN" "h3k9me3_ms_b_ENCSR983MKX" "h3k27ac_ms_b_ENCSR364OIK" "h3k27ac_ms_b_ENCSR541PLO" "h3k36me3_ms_cd4_ENCSR471WZE" "h3k36me3_ms_cd4_ENCSR473TGK" "h3k4me1_ms_nk_ENCSR482UGV" "h3k4me1_ms_nk_ENCSR718LCI" "h3k4me1_ms_nk_ENCSR806HUT" "h3k4me3_ms_nk_ENCSR234BAD" "h3k4me3_ms_nk_ENCSR703TWY" "h3k4me3_ms_nk_ENCSR755ZGS" "h3k9me3_ms_nk_ENCSR061ATV" "h3k9me3_ms_nk_ENCSR611YIA" "h3k9me3_ms_nk_ENCSR667HUI" "h3k27ac_ms_nk_ENCSR246TTM")

# Loop through all sample subdirectories
for sample_dir in "$input_dir"/*; do
  # Check if directory name matches any sample ID
  if [[ " ${sample_ids[@]} " =~ " $(basename "$sample_dir") " ]]; then
    echo "Processing sample: $(basename "$sample_dir")"

    # Create output directory for the sample (if it doesn't exist)
    mkdir -p "$output_dir/$(basename "$sample_dir")"

    # Find FASTQ files within the sample directory
    # Adjust the pattern below if your file names differ
    fastq_files=( "$sample_dir"/*.fastq.gz )

    # Check if any FASTQ files found
    if [[ ${#fastq_files[@]} -eq 0 ]]; then
      echo "Warning: No FASTQ files found in $sample_dir"
      continue # Skip to next iteration if no FASTQ files
    fi

    # BWA mem command with paired-end or single-end handling
    if [[ ${#fastq_files[@]} -eq 2 ]]; then
      # Paired-end reads (assuming *_1.fastq.gz and *_2.fastq.gz naming)
      bwa mem $reference_genome ${fastq_files[0]} ${fastq_files[1]} > "$output_dir/$(basename "$sample_dir").sam"
    else
      # Single-end read (assuming only one FASTQ file)
      # Extract only the sample name from the directory path
      sample_name=$(basename "$sample_dir" /)
      bwa mem $reference_genome ${fastq_files[0]} > "$output_dir/$sample_name.sam"
    fi
  fi
done

```

<br>

### samtools (sam-bam transition):

```

#!/bin/bash
#SBATCH --job-name=sam_to_bam
#SBATCH --nodes=1
#SBATCH --account=users
#SBATCH --qos=cuda
#SBATCH --partition=cuda
#SBATCH --cpus-per-task=16
#SBATCH -o /cta/users/merve.kaftancioglu/alternative_scenario/samtools_as/out-err/%x-%j-slurm.out
#SBATCH -e /cta/users/merve.kaftancioglu/alternative_scenario/samtools_as/out-err/%x-%j-slurm.err

# Load samtools module
module load samtools-1.9-gcc-9.2.0-w7pulwi

# Define directories (consider adjusting paths)
input_dir="/cta/users/merve.kaftancioglu/alternative_scenario/bwa_as/output"
output_dir="/cta/users/merve.kaftancioglu/alternative_scenario/samtools_as/output"

# Loop through all subdirectories in the input directory
for subdir in $input_dir/*; do
  # Check if it's a directory (avoid hidden folders)
  if [ -d "$subdir" ]; then
    # Get the sample name from the subdirectory name
    sample_name=$(basename "$subdir")

    # Look for a .sam file within the subdirectory
    sam_file="$subdir/$sample_name.sam"

    # Check if the SAM file exists
    if [ -f "$sam_file" ]; then
      # Create a folder with the sample name in the output directory
      mkdir -p "$output_dir/$sample_name"

      # Convert SAM to BAM using samtools view and place it in the new folder
      samtools view -bS "$sam_file" > "$output_dir/$sample_name/$sample_name.bam"

      # Check for conversion errors (optional)
      if [ $? -ne 0 ]; then
        echo "Error converting $sam_file to BAM. See slurm error logs for details."
      fi
    fi
  fi
done

echo "Done!")
```
<br>

### picard (didn't work)

```
#!/bin/bash
#SBATCH --job-name=picard
#SBATCH --nodes=1
#SBATCH --account=users
#SBATCH --qos=cuda
#SBATCH --partition=cuda
#SBATCH --cpus-per-task=16
#SBATCH -o /cta/users/merve.kaftancioglu/alternative_scenario/picard_as/out-err/%x-%j-slurm.out
#SBATCH -e /cta/users/merve.kaftancioglu/alternative_scenario/picard_as/out-err/%x-%j-slurm.err

# Load BWA module

# Load required modules
module load picard samtools

# Define input and output file variables
MAIN_DIR="/cta/users/merve.kaftancioglu/alternative_scenario/samtools_as/output"  # Replace with path to main directory
REALIGNED_BAMS_DIR="$MAIN_DIR/realigned_bams"  # Output directory for realigned BAMs
TEMP_DIR="$MAIN_DIR/temp"  # Optional temporary directory (create if needed)

# **Step 1: Mark duplicates on each individual BAM**

# Loop through subdirectories
for subdir in $MAIN_DIR/*; do
  # Check if subdirectory is a directory (avoid hidden files)
  if [[ -d "$subdir" ]]; then
    # Extract subdirectory name
    subdir_name=$(basename "$subdir")

    # Find BAM files with wildcard
    bam_file="$subdir/$subdir_name*.bam"

    # Check if BAM file exists
    if [[ -f "$bam_file" ]]; then
      # Mark duplicates and place output in subdirectory
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

# **Step 2: Merge marked duplicate BAMs (optional)**

# List of temporary marked duplicate BAMs
marked_dups_bams="$TEMP_DIR/$SAMPLE_NAME*.marked_dups.bam"

# Merge marked duplicate BAMs
merge_bam_cmd="picard MergeSamFiles \
  O=$REALIGNED_BAMS_DIR/$SAMPLE_NAME.merged.marked_dups.bam \
  $marked_dups_bams"

echo "Merging marked duplicate BAMs..."
$merge_bam_cmd

# **Step 3: Re-mark duplicates on the merged BAM**

# Re-mark duplicates on the merged BAM
re_mark_dups_cmd="picard MarkDuplicates \
  I=$REALIGNED_BAMS_DIR/$SAMPLE_NAME.merged.marked_dups.bam \
  O=$REALIGNED_BAMS_DIR/$SAMPLE_NAME.realigned.bam \
  M=$REALIGNED_BAMS_DIR/$SAMPLE_NAME.realigned_metrics.txt \
  CREATE_INDEX=true"

echo "Re-marking duplicates on merged BAM..."
$re_mark_dups_cmd

# **Optional Step: Clean up temporary files (if used)**
# rm -rf $TEMP_DIR  # Uncomment to remove temporary files

echo "Done!"

# **Step 3: Re-mark duplicates on merged BAMs (optional)**

# If you performed merging, adjust commands to work with merged BAMs.

# **Optional Cleanup (if using temporary directory)**
# rm -rf $TEMP_DIR  # Uncomment to remove temporary files

echo "Done!"

```

<br>

### Docker (used locally, containers that was tried):


```

docker pull nfcore/base

docker pull pegi3s/igv:latest

docker pull staphb/multiqc

docker pull nanozoo/multiqc

docker pull zavolab/multiqc


```

<br>

### nextflow (recently added to Tosun, didn't work):

```

module load miniconda/22.9
eval "$(conda shell.bash hook)"
conda activate nextflow
nextflow

```

<br>

### singularity (used with nf-core, didn't work):

```

module spider singularity/3.8.4

```

<br> 

### snakemake (used with chocolate ChIP, didnt work):


```

snakemake-5.23.0-gcc-9.2.0-awyfooy

```
