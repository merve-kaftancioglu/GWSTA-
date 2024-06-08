# GWSTA- Project Report

# Epigenetic differences of Retromer proteins for immune cells in Multiple Sclerosis

## Merve Kaftancioglu

## Introduction

This project investigates potential epigenetic alterations in immune cells of Multiple Sclerosis (MS) patients, specifically focusing on recycling proteins called retromers in CD4+ T cells, CD8+ T cells, B cells, and NK cells.
As the most common chronic inflammatory disease of the central nervous system (CNS), it is estimated that MS is affecting over 2 million people globally. Although its cause is likely a complex combination of genetic and environmental factors, it is certainly known that the immunopathological nature of the disease triggers the neural degeneration. Adaptive immune system components, T and B cells, attack the myelin sheath that covers nerve fibers and thus lead to various neurological symptoms: Specific types of T cells, particularly Th1 and Th17, along with CD8+ and regulatory T cells create an unusual immune response. On the other hand, B cells and their antibodies are also involved in further damage in the inflammatory process (Filippi et al., 2018a; Haki et al., 2024; Jakimovski et al., 2024a; Pomales, 1990; Reich et al., 2018).

The complexity of the disease's underlying causes inevitably bring a possibility of epigenetic mechanisms in minds. Indeed, histone modifications might influence MS related genes in three ways: 1. making them more susceptible to MS: Certain genes linked to MS risk might be regulated by histone modifications, 2. changing normal mechanisms of autoimmune system: Histone modifications might influence the differentiation of several types of T cells and their functionality, and 3. regulating monocyte/macrophage functions: Histone methylation might be involved in the polarization of monocytes/macrophages towards a pro-inflammatory or anti-inflammatory state. As studies indicate, abnormal epigenetic motifs in the immune cells of MS patients definitely highlight crucial roles of epigenetics in this disease (Gacias & Casaccia, n.d.; Li et al., 2017; Ruhrmann, 2017).

Retromers are protein complexes which are found in eukaryotic cells and they are in charge of recycling other proteins. Retromer complexes help movement of proteins from endosomes either back to the Golgi network or the cell membrane. This regulation prevents proteins from being broken down and allows them to be used again and again. Due to this function, retromers are essential for continuation of various cell processes, particularly the maintenance of neuronal functions. Interestingly, aberrant functioning of retromers has been linked to several neurodegenerative diseases such as Alzheimer’s and Parkinson(Belenkaya et al., 2008; Cui et al., 2018, 2019; Gallon & Cullen, 2015; Wang & Bellen, 2015).

Here, the aim of the project is to enlighten the possible epigenetic changes of immune cells in MS, with a particular focus on retromer proteins, which is hoped to bring a new angle to investigate this disease. The ENCODE data that has been gathered actually belongs to a research (NIH Grant No: UM1HG009390) which was conducted by the members of Bradley Bernstein's lab in Broad Institute. For this study, Mint ChIP Sequencing has been preferred since it enables to barcode multiple samples together, which has a definite advantage over regular ChIP-Seq considering its large scope. Among many immune cell samples, CD4+ T cells, CD8+ T cells, and B cells are chosen for this study. Also, to see its unique role among the innate components, NK cells are also included. Histone modifications like H3K27ac, H3K27me3, H3K36me3, H3K4me1, H3K4me3, and H3K9me3 are decided to be the scope of the modifications. A comparison of these histone modifications in MS patients versus healthy subjects is expected to reveal significant differences and to provide compelling evidence for the involvement of epigenetic dynamics in the development of this disease (Abascal et al., 2020; Schreiber et al., 2023; Walter et al., 2019).

## Results

Computational analysis for this project was limited due to time constraints and technical issues. Consequently, conclusive results haven't been generated in a general context. However, the stepwise progress is presented here (A). Note that, to expedite the pipeline, samples were divided into smaller tasks based on processing time and run in parallel. 

After determination of the topic and realization of the time factor, beside downloading the raw versions of the chosen data, already processed versions of the samples which were provided in their respective ENCODE pages are also downloaded parallelly. There were two reasons behind this step, first while continuing conventional way of analyzing, utilizing this data and automatized pipelines such as nf-core to have a tangible result at the end, second, making comparisons of the processed and the personal data, and by that way, understanding the different options and settings of the tools from experienced researches in some kind of a reverse engineering. Unfortunately, this attempt also failed, but its unsuccessful results will also be shared (B), as well.
(A) 1st FastQC analysis: At the total 153 samples were evaluated following their download process. Other than having some CG content warnings in a few samples, the data quality was top notch. Randomly chosen samples revealed a potential adapter region of 15-20bp, which might affect further analysis.

Trimming step: Trimmomatic and Cutadapt were tested. Initially, specific adapter lists were used to make more exact trims, but mismatches with actual adapters or scripting errors occurred. Consequently, Trimmomatic was employed with a ~20bp trim based on random sample evaluation.

2nd FastQC analysis and MultiQC: FastQC data was evaluated with MultiQC due to the critical nature of this step.  The partial nature of the samples necessitated running all 4732 gz files together. Important features like per base sequence content and quality, per sequence quality,  per sequence GC content, sequence duplication levels, overrepresented sequences, adapter content and status check heatmap were checked. As mentioned earlier, other than some high GC content samples appeared. Also, based on adapter content analyze results, it can be deduced that ~20bp should have been extended to 25bp to have a clearer analyzing process in the further steps. Other than this, the overall data quality was good based on the provided result report. 

Alignment Step: For this job, three different software were used: BWA, Bowtie2 and STAR. Hisat2 and NextGenMap were also considered, but it was learnt that algorithm of Hisat2 was better for shorter sequences, and NextGenMap was not provided in the local slurm. As mentioned above, separated files were a problem in directory detections, thus, sometimes wildcard option and for loops did not agreed with the tool settings. This was especially apparent with STAR trials. On the other hand, although Bowtie2 was a more preferable option since it also provided trimming function, several batch script trials were unsuccessful. At the end, BWA mem was used to generate sam files. An important challenge BWA threw was unpredictable nature of parallel runs, sometimes when a specific sample did not progress, it was tried to be added as another task group script, yet, after a while stalling samples decided to progress as well, and at the and same sample began to run parallelly in two different scripts. If they used to produce a second sam file this would not be considered as an obstacle, but somehow the data began to add up to the same file, and thus, unusually large files began to appear. For that reason, the BWA scripts were run a second time, and waited till the end of job queue to re-run undone samples in the batches. 
samtools: This time, sam files were transferred into bam files. Fortunately, this process took significantly less time compared to the other steps.

IGV visualization: Again, the size of the sample data and the time limitation became a problem for visualization of all the bam data. Some selected samples were tried to generate some images, however with 64GB size bam files, the waiting period became a tedious issue.

(B) As mentioned above, besides providing their whole raw data, protocols and analysis processes, Bernstein Lab also added several data types which have already been processed, such as bam, bigWig, bigBed narrowPeak, and bed narrowPeak. This data gives a great chance to evaluate choices of the experienced researchers for the analyzes, an information which can be gained only after a lot of experience and digestion of hundreds of articles. Since these kind of tips and tricks only can be understood after comparison with personal analyses, they are downloaded manually while the other scripts were running in the cluster. Ready-to-use pipeline structures such as nf-core are also invaluable tools for faster and better results. However, their highly dependent nature for other tools makes this process a bit more complicated. For example, nf-core chipseq pipeline is utilizing nextflow, a workflow system which was not provided in Tosun. For that reason, as an alternative, snakemake was tried. Moreover, nf-core also needed virtual and independent computational containers such as Docker and Singularity. Earlier experiences with Docker leaded an assumption to complete tasks much earlier, unfortunately containers with on the current Linux system (ubuntu jammy) did not work, after pulling them and run, images immediately got disconnected. When plans with local Docker were fruitless, Singularity was tried. This time, a ready chipseq pipeline called chocolate was planned to be adapted, however for it, the required Singularity version for that was not compatible with the one HPC provided since it was very old. At the end, nextflow was requested to the Tosun, but it only came on the day of the deadline of the report.
 
![Figure 1: General idea of the project](https://github.com/merve-kaftancioglu/GWSTA-/assets/160722990/96af95ac-678d-4792-91e4-0e5a9c2abe4a)

Figure 1: General idea of the project

![Figure 2: Pipeline representation (modified from nf-core chipseq page)](https://github.com/merve-kaftancioglu/GWSTA-/assets/160722990/db916ca4-df60-4473-b4f2-a9875d6b2595)

Figure 2: Pipeline representation (modified from nf-core chipseq page)
 
![Figure 3: Used software and data list from the cluster](https://github.com/merve-kaftancioglu/GWSTA-/assets/160722990/b786d462-fcbe-470d-beaf-d94f158fa0a7)

Figure 3: Used software and data list from the cluster

![Figure 4: Total adapter content plot from MultiQC](https://github.com/merve-kaftancioglu/GWSTA-/assets/160722990/9416f711-a63f-404c-9c30-fd76d393f129)

Figure 4: Total adapter content plot from MultiQC

## Discussion

The development of the project process was a quite informative, specifically in terms of the nature of bioinformatic analyses, computational experience, and inevitable consequences of human factors. During the research process many online platforms were used, and various researchers guided and suggested different approaches. That’s why, after the end of semester, full complication of the project is planned, and noted suggestions are also hoped to contribute overall understanding of the pipelines and tools. For example, footprinting analysis on main components of retromer proteins like VPS35 (Ch. 16-11.2) and VPS26 (Ch. 10-22.1) can be further investigated for their transcription factor bindings via HOMER. Genrich was another tool which was highly recommended for peak detections. From the readings, one can gather that Galaxy is a great tool definitely worth to discover. 

One of the greatest challenges should be noted is the batch script generation. Due to differences of the tools, some code patterns do not work, and relying solely on Python or R knowledge is a serious limitation. In this stage, several AI services such as Gemini, Copilot became a vital contributor for correcting and improving the scripts.

The theoretical background of tools and their performances requires some feedback from experienced people, as well as trials with the alternatives at the same time to compare their strong and weak features. In this matter, the development of nf-core is very useful while enabling alternatives and supporting to make personal changes on the general pipeline, although relying on many other tools. From a personal point of view, the active and easy to researchers who contribute progress of nf-core are actually helping the newbie bioinformatians by giving solid structural idea about sequencing technologies instead of making them to reinvent the tire. Furthermore, as reproducibility in bioinformatics research is turning into a serious issue, it can be deduced that well-structured pipelines on artificial environments will be more favorable in the near future.

Although troubleshooting processes turned into very informative lessons, they can be very time consuming if the problematic part cannot be determined, or at least estimated. In these kind of cases, guiding the right way without fully addressing the solution/issue is the essence for understanding and interpreting analyses. 

Aside from time constraints and technical issues that prevented conclusive results, outcomes of this project can be states as: successful data quality evaluation using FastQC and MultiQC, understanding building different script structures for BWA, Bowtie2, STAR, and gaining more experience with containerization tools. 

### Last edit (8.6.2024): The full sample bash job (nf-core chipseq) is working more than 8 days, however, due to cluster speed it is still going on and it might take more than 15 days if this situation continues. That's why, I used sakura for cell and histone type based small jobs. So far, h3km3me1 b cells and cd4 cells are completed (cd4s ended after 12pm, but I wanted to compare these two batches, so, last modification of the report is done on 8th of June, 5am) Partially downloaded output report: 

pdf: [MultiQC_Report_h3k4me1_b.pdf](https://github.com/user-attachments/files/15746230/MultiQC_Report_h3k4me1_b.pdf)

igv broad peaks (bigWig files) for VPS35 region (ch16: 46,656,130 - 46,689,519), one of the most imporant component of the retromer complex: ![h3k4me1_b_igv](https://github.com/merve-kaftancioglu/GWSTA-/assets/160722990/2a1569c1-1ae2-4b3a-b45e-9cbc01c725ec)


## Materials & Methodology

To keep an organized structure for the project data, each and every sample was labelled, and they are kept separately in folders, thus, for loops were utilized for input and output directories. Different numbers of scripts were created for parallel run based on the process speed of the tool. For example, for FastQC 31 batches were made, compared to BWA which required only 5. 
The versions of the tools that are used: 

System: Linux Ubuntu jammy release cran40

Tools: FastQC v 0.11.7

MultiQC (miniconda v. 22.9)

Trimmomatic v. 0.38

Cutadapt v. 4.0

STAR v. 2.7

Bowtie2 v. 2.2.3.5.1

BWA v. 0.7.17

nf-core ChIP-seq 2.0

samtools v. 1.9

picard v. 2.20.8

Docker v. 26.1.2

Nextflow (miniconda v. 22.9)

Singularity v. 3.8.4

Snakemake v. 5.23

IGV v. 2.17.3

## References 

Abascal, F., Acosta, R., Addleman, N. J., Adrian, J., Afzal, V., Aken, B., Akiyama, J. A., Jammal, O. Al, Amrhein, H., Anderson, S. M., Andrews, G. R., Antoshechkin, I., Ardlie, K. G., Armstrong, J., Astley, M., Banerjee, B., Barkal, A. A., Barnes, I. H. A., Barozzi, I., … Zimmerman, J. (2020). Expanded encyclopaedias of DNA elements in the human and mouse genomes. Nature, 583(7818), 699–710. https://doi.org/10.1038/s41586-020-2493-4

Ahmed, T. (2022). Neural stem cell engineering for the treatment of multiple sclerosis. Biomedical Engineering Advances, 4, 100053. https://doi.org/10.1016/j.bea.2022.100053

Andersen, M. H. (2017). Anti-regulatory T cells. In Seminars in Immunopathology (Vol. 39, Issue 3, pp. 317–326). Springer Verlag. https://doi.org/10.1007/s00281-016-0593-x

Asashima, H., Axisa, P. P., Pham, T. H. G., Longbrake, E. E., Ruff, W. E., Lele, N., Cohen, I., Raddassi, K., Sumida, T. S., & Hafler, D. A. (2022). Impaired TIGIT expression on B cells drives circulating follicular helper T cell expansion in multiple sclerosis. Journal of Clinical Investigation, 132(20). https://doi.org/10.1172/JCI156254

Attfield, K. E., Jensen, L. T., Kaufmann, M., Friese, M. A., & Fugger, L. (2022). The immunology of multiple sclerosis. In Nature Reviews Immunology (Vol. 22, Issue 12, pp. 734–750). Nature Research. https://doi.org/10.1038/s41577-022-00718-z

Barnett, M. H., Henderson, A. P. D., & Prineas, J. W. (2006). The macrophage in MS: Just a scavenger after all? Pathology and pathogenesis of the acute MS lesion. In Multiple Sclerosis (Vol. 12, Issue 2, pp. 121–132). https://doi.org/10.1191/135248506ms1304rr

Belenkaya, T. Y., Wu, Y., Tang, X., Zhou, B., Cheng, L., Sharma, Y. V., Yan, D., Selva, E. M., & Lin, X. (2008). The Retromer Complex Influences Wnt Secretion by Recycling Wntless from Endosomes to the Trans-Golgi Network. Developmental Cell, 14(1), 120–131. https://doi.org/10.1016/j.devcel.2007.12.003

Bjornevik, K., Cortese, M., Healy, B. C., Kuhle, J., Mina, M. J., Leng, Y., Elledge, S. J., Niebuhr, D. W., Scher, A. I., Munger, K. L., & Ascherio, A. (2022a). MULTIPLE SCLEROSIS Longitudinal analysis reveals high prevalence of Epstein-Barr virus associated with multiple sclerosis. In Science (Vol. 375). https://www.science.org

Boccazzi, M., Raffaele, S., & Fumagalli, M. (2022). Not only myelination: the immune-inflammatory functions of oligodendrocytes. Neural Regeneration Research, 17(12), 2661–2663. https://doi.org/10.4103/1673-5374.342678
Cheru, N., Hafler, D. A., & Sumida, T. S. (2023). Regulatory T cells in peripheral tissue tolerance and diseases. In Frontiers in Immunology (Vol. 14). Frontiers Media S.A. https://doi.org/10.3389/fimmu.2023.1154575

Chipseq: Introduction. chipseq: Introduction. (n.d.). https://nf-co.re/chipseq/2.0.0

Chomyk, A. M., Volsko, C., Tripathi, A., Deckard, S. A., Trapp, B. D., Fox, R. J., & Dutta, R. (2017). DNA methylation in demyelinated multiple sclerosis hippocampus. Scientific Reports, 7(1). https://doi.org/10.1038/s41598-017-08623-5

Chu, Y., Dai, E., Li, Y., Han, G., Pei, G., Ingram, D. R., Thakkar, K., Qin, J. J., Dang, M., Le, X., Hu, C., Deng, Q., Sinjab, A., Gupta, P., Wang, R., Hao, D., Peng, F., Yan, X., Liu, Y., … Wang, L. (2023). Pan-cancer T cell atlas links a cellular stress response state to immunotherapy resistance. Nature Medicine, 29(6), 1550–1562. https://doi.org/10.1038/s41591-023-02371-y

Cui, Y., Carosi, J. M., Yang, Z., Ariotti, N., Kerr, M. C., Parton, R. G., Sargeant, T. J., & Teasdale, R. D. (2019). Retromer has a selective function in cargo sorting via endosome transport carriers. Journal of Cell Biology, 218(2), 615–631. https://doi.org/10.1083/jcb.201806153

Cui, Y., Yang, Z., & Teasdale, R. D. (2018). The functional roles of retromer in Parkinson’s disease. In FEBS Letters (Vol. 592, Issue 7, pp. 1096–1112). Wiley Blackwell. https://doi.org/10.1002/1873-3468.12931
Erra Díaz, F., Dantas, E., & Geffner, J. (2018). Unravelling the interplay between extracellular acidosis and immune cells. In Mediators of Inflammation (Vol. 2018). Hindawi Limited. https://doi.org/10.1155/2018/1218297

Eva, L., Pleș, H., Covache-Busuioc, R. A., Glavan, L. A., Bratu, B. G., Bordeianu, A., Dumitrascu, D. I., Corlatescu, A. D., & Ciurea, A. V. (2023b). A Comprehensive Review on Neuroimmunology: Insights from Multiple Sclerosis to Future Therapeutic Developments. In Biomedicines (Vol. 11, Issue 9). Multidisciplinary Digital Publishing Institute (MDPI). https://doi.org/10.3390/biomedicines11092489

Fathallah, S., Abdellatif, A., & Saadeldin, M. (2023). Unleashing nature’s potential and limitations: Exploring molecular targeted pathways and safe alternatives for the treatment of multiple sclerosis (Review). Medicine International, 3(5). https://doi.org/10.3892/mi.2023.102

Filippi, M., Bar-Or, A., Piehl, F., Preziosa, P., Solari, A., Vukusic, S., & Rocca, M. A. (2018a). Multiple sclerosis. Nature Reviews Disease Primers, 4(1). https://doi.org/10.1038/s41572-018-0041-4

Freedman, M. S., Bar-Or, A., Atkins, H. L., Karussis, D., Frassoni, F., Lazarus, H., Scolding, N., Slavin, S., Le Blanc, K., & Uccelli, A. (2010). The therapeutic potential of mesenchymal stem cell transplantation as a treatment for multiple sclerosis: Consensus report of the international MSCT study group. Multiple Sclerosis, 16(4), 503–510. https://doi.org/10.1177/1352458509359727

Gacias, M., & Casaccia, P. (n.d.). EPIGENETIC MECHANISMS IN MULTIPLE SCLEROSIS HHS Public Access.

Gallon, M., & Cullen, P. J. (2015). Retromer and sorting nexins in endosomal sorting. Biochemical Society Transactions, 43, 33–47. https://doi.org/10.1042/BST20140290

Guo, L., Li, X., Gould, T., Wang, Z. Y., & Cao, W. (2023). T cell aging and Alzheimer’s disease. In Frontiers in Immunology (Vol. 14). Frontiers Media S.A. https://doi.org/10.3389/fimmu.2023.1154699

Haghani, A., Li, C. Z., Robeck, T. R., Zhang, J., Lu, A. T., Ablaeva, J., Acosta-Rodríguez, V. A., Adams, D. M., Alagaili, A. N., Almunia, J., Aloysius, A., Amor, N. M., Ardehali, R., Arneson, A., Scott Baker, C., Banks, G., Belov, K., Bennett, N. C., Black, P., … Horvath, S. (2023). DNA methylation networks underlying mammalian traits. In Science (Vol. 381). https://www.science.org

Haki, M., Al-Biati, H. A., Al-Tameemi, Z. S., Ali, I. S., & Al-Hussaniy, H. A. (2024). Review of multiple sclerosis: Epidemiology, etiology, pathophysiology, and treatment. Medicine (United States), 103(8), E37297. https://doi.org/10.1097/MD.0000000000037297

Hammack, B. N., Fung, K., Hunsucker, S. W., Duncan, M. W., Burgoon, M. P., Owens, G. P., Gilden, D. H., & Hammack, B. N. (2004b). Proteomic analysis of multiple sclerosis cerebrospinal fluid. In Multiple Sclerosis (Vol. 10).

He, H., Hu, Z., Xiao, H., Zhou, F., & Yang, B. (2018). The tale of histone modifications and its role in multiple sclerosis. In Human Genomics (Vol. 12, Issue 1). BioMed Central Ltd. https://doi.org/10.1186/s40246-018-0163-5

Hecht, V., Dong, K., Rajesh, S., Shpilker, P., Wekhande, S., & Shoresh, N. (2023). Analyzing histone ChIP-seq data with a bin-based probability of being signal. PLoS Computational Biology, 19(10 October). https://doi.org/10.1371/journal.pcbi.1011568

Jakimovski, D., Bittner, S., Zivadinov, R., Morrow, S. A., Benedict, R. H., Zipp, F., & Weinstock-Guttman, B. (2024a). Multiple sclerosis. In The Lancet (Vol. 403, Issue 10422, pp. 183–202). Elsevier B.V. https://doi.org/10.1016/S0140-6736(23)01473-3

Jones, T., Naslavsky, N., & Caplan, S. (2020). Eps15 Homology Domain Protein 4 (EHD4) is required for Eps15 Homology Domain Protein 1 (EHD1)-mediated endosomal recruitment and fission. PLoS ONE, 15(9 September 2020). https://doi.org/10.1371/journal.pone.0239657

Kaaij, L. J. T., Mokry, M., Zhou, M., Musheev, M., Geeven, G., Melquiond, A. S. J., de Jesus Domingues, A. M., de Laat, W., Niehrs, C., Smith, A. D., & Ketting, R. F. (2016). Enhancers reside in a unique epigenetic environment during early zebrafish development. Genome Biology, 17(1). https://doi.org/10.1186/s13059-016-1013-1

Kornek, B., & Lassmann, H. (2003). Neuropathology of multiple sclerosis - New concepts. Brain Research Bulletin, 61(3), 321–326. https://doi.org/10.1016/S0361-9230(03)00095-9

Kosa, P., Barbour, C., Varosanec, M., Wichman, A., Sandford, M., Greenwood, M., & Bielekova, B. (2022). Molecular models of multiple sclerosis severity identify heterogeneity of pathogenic mechanisms. Nature Communications, 13(1). https://doi.org/10.1038/s41467-022-35357-4

Kouli, A., Jensen, M., Papastavrou, V., Scott, K. M., Kolenda, C., Parker, C., Solim, I. H., Camacho, M., Martin-Ruiz, C., & Williams-Gray, C. H. (2021). T lymphocyte senescence is attenuated in Parkinson’s disease. Journal of Neuroinflammation, 18(1). https://doi.org/10.1186/s12974-021-02287-9

Li, H. (n.d.). Which human reference genome to use?. Sitewide ATOM. https://lh3.github.io/2017/11/13/which-human-reference-genome-to-use 

Li, X., Xiao, B., & Chen, X. S. (2017). DNA Methylation: a New Player in Multiple Sclerosis. In Molecular Neurobiology (Vol. 54, Issue 6, pp. 4049–4059). Humana Press Inc. https://doi.org/10.1007/s12035-016-9966-3

Lincoln, M. R., Axisa, P. P., & Hafler, D. A. (2020). Epigenetic fine-mapping: identification of causal mechanisms for autoimmunity. In Current Opinion in Immunology (Vol. 67, pp. 50–56). Elsevier Ltd. https://doi.org/10.1016/j.coi.2020.09.002

Lippi, A., & Krisko, A. (2023). Protein aggregation: A detrimental symptom or an adaptation mechanism? In Journal of Neurochemistry. John Wiley and Sons Inc. https://doi.org/10.1111/jnc.15955

Lucca, L. E., Lerner, B. A., Park, C., Debartolo, D., Harnett, B., Kumar, V. P., Ponath, G., Raddassi, K., Huttner, A., Hafler, D. A., & Pitt, D. (2020). Differential expression of the T-cell inhibitor TIGIT in glioblastoma and MS. Neurology: Neuroimmunology and NeuroInflammation, 7(3). https://doi.org/10.1212/NXI.0000000000000712

Martin, M., Vermeiren, S., Bostaille, N., Eubelen, M., Spitzer, D., Vermeersch, M., Profaci, C. P., Pozuelo, E., Toussay, X., Raman-Nair, J., Tebabi, P., America, M., de Groote, A., Sanderson, L. E., Cabochette, P., Germano, R. F., Torres, D., Boutry, S., de 
Kerchove, A., … Vanhollebeke, B. (n.d.). Engineered Wnt ligands enable blood-brain barrier repair in neurological disorders.

Noble, W., & Hanger, D. P. (2023). Trimming away tau in neurodegeneration. Science, 381(6656), 377–378. https://doi.org/10.1126/science.adj0256

Park, P. J. (2009). ChIP-seq: Advantages and challenges of a maturing technology. In Nature Reviews Genetics (Vol. 10, Issue 10, pp. 669–680). https://doi.org/10.1038/nrg2641

Pomales, Z. A. (1990). Multiple sclerosis: a review. In Boletín de la Asociación Médica de Puerto Rico (Vol. 82, Issue 2, pp. 67–71). https://doi.org/10.32628/ijsrst218131

Pompura, S. L., Hafler, D. A., & Dominguez-Villar, M. (2022). Fatty Acid Metabolism and T Cells in Multiple Sclerosis. In Frontiers in Immunology (Vol. 13). Frontiers Media S.A. https://doi.org/10.3389/fimmu.2022.869197

Raab, J. R., & Kamakaka, R. T. (2010). Insulators and promoters: Closer than we think. In Nature Reviews Genetics (Vol. 11, Issue 6, pp. 439–446). https://doi.org/10.1038/nrg2765

Reich, D. S., Lucchinetti, C. F., & Calabresi, P. A. (2018). Multiple Sclerosis. New England Journal of Medicine, 378(2), 169–180. https://doi.org/10.1056/nejmra1401483

Reitz, C. (2015). The role of the retromer complex in aging-related neurodegeneration: a molecular and genomic review. In Molecular Genetics and Genomics (Vol. 290, Issue 2, pp. 413–427). Springer Verlag. https://doi.org/10.1007/s00438-014-0939-9

Roy, R., Ramamoorthy, S., Shapiro, B. D., Kaileh, M., Hernandez, D., Sarantopoulou, D., Arepalli, S., Boller, S., Singh, A., Bektas, A., Kim, J., Moore, A. Z., Tanaka, T., McKelvey, J., Zukley, L., Nguyen, C., Wallace, T., Dunn, C., Wersto, R., … Sen, R. (2021). DNA methylation signatures reveal that distinct combinations of transcription factors specify human immune cell epigenetic identity. Immunity, 54(11), 2465-2480.e5. https://doi.org/10.1016/j.immuni.2021.10.001

Rudick, R. A., Fisher, E., Lee, J.-C., Duda, J. T., & Simon, J. (n.d.). Brain atrophy in relapsing multiple sclerosis: relationship to relapses, EDSS, and treatment with interferon b-1a. www.arnoldpublishers.com/journals

Ruhrmann, S. (2017). INVOLVEMENT OF EPIGENETIC MECHANISMS IN DISEASE INHERITANCE AND PATHOGEGENESIS OF MULTIPLE SCLEROSIS (MS) WITH A FOCUS ON GENOMIC IMPRINTING AND DNA METHYLATION IN CD4 + T CELLS.

Schattgen, S. A., & Thomas, P. G. (2021). T RH cells, helpers making an impact in their local community Downloaded from. Sci. Immunol, 6, 2886. https://doi.org/10.1101/2020.02.28.970400

Schreiber, J., Boix, C., wook Lee, J., Li, H., Guan, Y., Chang, C. C., Chang, J. C., Hawkins-Hooker, A., Schölkopf, B., Schweikert, G., Carulla, M. R., Canakoglu, A., Guzzo, F., Nanni, L., Masseroli, M., Carman, M. J., Pinoli, P., Hong, C., Yip, K. Y., … Kundaje, A. (2023). The ENCODE Imputation Challenge: a critical assessment of methods for cross-cell type imputation of epigenomic profiles. Genome Biology, 24(1). https://doi.org/10.1186/s13059-023-02915-y

Sendinc, E., & Shi, Y. (2023). RNA m6A methylation across the transcriptome. In Molecular Cell (Vol. 83, Issue 3, pp. 428–441). Cell Press. https://doi.org/10.1016/j.molcel.2023.01.006

Shang, Z., Sun, W., Zhang, M., Xu, L., Jia, X., Zhang, R., & Fu, S. (2020). Identification of key genes associated with multiple sclerosis based on gene expression data from peripheral blood mononuclear cells. PeerJ, 2020(2). https://doi.org/10.7717/peerj.8357

Steinman, L. (1996a). Multiple Sclerosis: A Coordinated Review Immunological Attack against Myelin in the Central Nervous System. In Cell (Vol. 85).

Sun, Z., Kwon, J.-S., Ren, Y., Chen, S., Cates, K., Lu, X., Walker, C. K., Karahan, H., Sviben, S., Fitzpatrick, J. A. J., Valdez, C., Houlden, H., Karch, C. M., Bateman, R. J., Sato, C., Mennerick, S. J., Diamond, M. I., Kim, J., Tanzi, R. E., … Yoo, A. S. (2023). Endogenous recapitulation of Alzheimer’s disease neuropathology through human 3D direct neuronal reprogramming. BioRxiv : The Preprint Server for Biology. https://doi.org/10.1101/2023.05.24.542155

Tiwari, S. (2017). Recent Advancement in Methodology for Understanding Epigenetic Modifications. Journal of Clinical Epigenetics, 03(03). https://doi.org/10.21767/2472-1158.100055

van Galen, P., Viny, A. D., Ram, O., Ryan, R. J. H., Cotton, M. J., Donohue, L., Sievers, C., Drier, Y., Liau, B. B., Gillespie, S. M., Carroll, K. M., Cross, M. B., Levine, R. L., & Bernstein, B. E. (2016). A Multiplexed System for Quantitative Comparisons of Chromatin Landscapes. Molecular Cell, 61(1), 170–180. https://doi.org/10.1016/j.molcel.2015.11.003

Varadarajan, S. G., Hunyara, J. L., Hamilton, N. R., Kolodkin, A. L., & Huberman, A. D. (2022). Central nervous system regeneration. In Cell (Vol. 185, Issue 1, pp. 77–94). Elsevier B.V. https://doi.org/10.1016/j.cell.2021.10.029

Walter, L. D., Galen, V., Bernstein, B. E., & Epstein, C. B. (2019). Mint-ChIP3: A low-input ChIP-seq protocol using multiplexed chromatin and T7 amplification. https://doi.org/10.17504/protocols.io.wbefaje

Wang, S., & Bellen, H. J. (2015). The retromer complex in development and disease. Development (Cambridge), 142(14), 2392–2396. https://doi.org/10.1242/dev.123737

Weiner, H. L., & Cohen, J. A. (2002). Treatment of multiple sclerosis with cyclophosphamide: critical review of clinical and immunologic effects. https://doi.org/10.1191/1352458502ms790o

Wong, W. K., Yin, B., Lam, C. Y. K., Huang, Y., Yan, J., Tan, Z., & Wong, S. H. D. (2022). The Interplay Between Epigenetic Regulation and CD8+ T Cell Differentiation/Exhaustion for T Cell Immunotherapy. In Frontiers in Cell and Developmental Biology (Vol. 9). Frontiers Media S.A. https://doi.org/10.3389/fcell.2021.783227

Wu, Y. L., Lin, Z. J., Li, C. C., Lin, X., Shan, S. K., Guo, B., Zheng, M. H., Li, F., Yuan, L. Q., & Li, Z. hong. (2023). Epigenetic regulation in metabolic diseases: mechanisms and advances in clinical study. In Signal Transduction and Targeted Therapy (Vol. 8, Issue 1). Springer Nature. https://doi.org/10.1038/s41392-023-01333-7

Yong, X., Zhao, L., Deng, W., Sun, H., Zhou, X., Mao, L., Hu, W., Shen, X., Sun, Q., Billadeau, D. D., Xue, Y., & Jia, D. (2020). Mechanism of cargo recognition by retromerlinked SNX-BAR proteins. PLoS Biology, 18(3). https://doi.org/10.1371/journal.pbio.3000631

Yoshida, T. M., Wang, A., & Hafler, D. A. (2022). Basic principles of neuroimmunology. In Seminars in Immunopathology (Vol. 44, Issue 5, pp. 685–695). Springer Science and Business Media Deutschland GmbH. https://doi.org/10.1007/s00281-022-00951-7

You, Z., Wang, L., He, H., Wu, Z., Zhang, X., Xue, S., Xu, P., Hong, Y., Xiong, M., Wei, W., & Chen, Y. (2023). Mapping of clonal lineages across developmental stages in human neural differentiation. Cell Stem Cell, 30(4), 473-487.e9. https://doi.org/10.1016/j.stem.2023.02.007

YouTube. (2024, April 25). Biyoenformatik Analizlerde Tekrarlanabilirlik: “nextflow ve NF-core” workshop’u - phd. Kübra Narcı. YouTube. https://www.youtube.com/watch?v=AqNmIkoQrNo

Zhang, F., Aschenbrenner, D., Yoo, J. Y., & Zuo, T. (2022). The gut mycobiome in health, disease, and clinical applications in association with the gut bacterial microbiome assembly. In The Lancet Microbe (Vol. 3, Issue 12, pp. e969–e983). Elsevier Ltd. https://doi.org/10.1016/S2666-5247(22)00203-8

Zhang, X., Liu, T., Hou, X., Zhou, Z., Zhang, F., Ma, H., Wu, X., & Jiang, J. (2023). Exosomes secreted by mesenchymal stem cells delay brain aging b upregulating SIRT1 expression. Scientific Reports, 13(1). https://doi.org/10.1038/s41598-023-40543-5

Zierfuss, B., Larochelle, C., & Prat, A. (2024). Blood–brain barrier dysfunction in multiple sclerosis: causes, consequences, and potential effects of therapies. In The Lancet Neurology (Vol. 23, Issue 1, pp. 95–109). Elsevier Ltd. https://doi.org/10.1016/S1474-4422(23)00377-0
 
### Supplemental Data

Details of the samples: https://docs.google.com/spreadsheets/d/1cLsejNY6V_bZS2QRhoWMrBAu4JIZdUGq/edit#gid=2018365031

Full MultiQC analysis report: https://drive.google.com/file/d/1SPVfamqntS4RtfoqQ0XMC3Ut8fVsBFFn/view?usp=drive_link

pdf version of the report: https://drive.google.com/file/d/11abFgQmwUsYuG9pNh7RxAraQb3OWmNUS/view?usp=sharing
