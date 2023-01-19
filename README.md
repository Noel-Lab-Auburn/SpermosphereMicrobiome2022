# SpermosphereMicrobiome2022
Everything data analysis for the spermosphere microbiome of cotton and soybean!

File tree for this project. It contains all data output from HPC to run within R
R scripts are found within the R scripts directory


Description of R scripts 
1. [Preprocessing](R_Scripts/Preprocessing.R) - run this script first since it generates input data for the rest of the    scirpts 
2. [Diversity](R_Scripts/Diversity.R) - general microbial diversity statistics 
3. [Differential Abundance](R_Scripts/diff_abund.R) - differential abundance analysis 
4. [Network](R_Scripts/NetworkResultsPlot.R)- generating and plotting networks. 
5. [Transmission Analysis](R_Scripts/transmission_analysis.R) - Not actually in the         manuscript, but attempts to find OTUs only present in certain samples. this study was not really designed for this analysis, but still fun to look through! 

```bash
.
└── SpermosphereMicrobiome2022
    ├── Bacteria
    │   ├── 16s_taxonomy.csv
    │   ├── Bacteria_spermosphere_CSS_112922.rds
    │   ├── Bacteria_spermosphere_nonnorm_112922.rds
    │   ├── otu_table_16s.csv
    │   ├── otus_16s.fasta
    │   ├── otus_16s.tre
    │   ├── otus_16s_midpoint.tre
    │   └── spermospheremetadata.csv
    ├── Differential_Abundance
    │   ├── differential_abund_121422.rds
    │   └── differential_abund_alloutput_121422.rds
    ├── Figures
    │   ├── Main Figures
    │   │   ├── PDF
    │   │   │   ├── Figure1_method.pdf
    │   │   │   ├── Figure2_imbibitiondiversity.pdf
    │   │   │   ├── Figure3_composition_top_20.pdf
    │   │   │   ├── Figure4_PCOA.pdf
    │   │   │   ├── Figure5_Diff_abundance2.pdf
    │   │   │   └── Figure6_Networks2.pdf
    │   │   └── PNG
    │   │       ├── Figure1_method.png
    │   │       ├── Figure2_imbibitiondiversity.png
    │   │       ├── Figure3_composition_top_20.png
    │   │       ├── Figure4_PCOA.png
    │   │       ├── Figure5_Diff_abundance2.png
    │   │       └── Figure6_Networks2.png
    │   └── SupplementalFigures
    │       ├── PDF
    │       │   ├── SupplementalFig1_PreliminaryLogCFU.pdf
    │       │   ├── SupplementalFig2_SequencingPerformance_120522.pdf
    │       │   ├── SupplementalFig3.pdf
    │       │   └── SupplementalFig4_composition_top_20_fungi.pdf
    │       └── PNG
    │           ├── SupplementalFig1_PreliminaryLogCFU.png
    │           ├── SupplementalFig2_SequencingPerformance_120522.png
    │           ├── SupplementalFig3.png
    │           └── SupplementalFig4_composition_top_20_fungi.png
    ├── Fungi
    │   ├── Fungi_CSSNorm_083022.rds
    │   ├── Fungi_spermosphere_unedited_083022.rds
    │   ├── METADATA.csv
    │   ├── OTU_Table.csv
    │   ├── fungal_taxonomy_DADA2_NBC.csv
    │   └── otus_R1.fasta
    ├── Networks
    │   ├── cotton_spieceasi_network.rds
    │   ├── netCompare_cotton_soil.rds
    │   ├── netCompare_cottonvssoybean.rds
    │   ├── netCompare_soybean_soil.rds
    │   ├── netConstruct_cottonvsoil.rds
    │   ├── netConstruct_cottonvssoybean.rds
    │   ├── netConstruct_soybeanvssoil.rds
    │   ├── netanalyse_cottonvsoil.rds
    │   ├── netanalyse_cottonvssoybean.rds
    │   ├── netanalyse_soybeanvssoil.rds
    │   ├── soil_spieceasi_network.rds
    │   └── soybean_spieceasi_network.rds
    ├── README.md
    ├── R_Scripts
    │   ├── Diversity.R
    │   ├── NetworkResultsPlot.R
    │   ├── Preprocessing.R
    │   ├── diff_abund.R
    │   └── transmission_analysis.R
    ├── SpermosphereMicrobiome2022.Rproj
    └── Tables
        ├── Main Tables
        │   ├── Table1_PERMANOVA_Prok_Split_120622.xlsx
        │   ├── Table2_CentralityNetworkHub.xlsx
        │   └── Table3_HubTaxa.xlsx
        └── Supplemental
            ├── SupplementalTable1_GlobalPERMANOVA_112922.xlsx
            ├── SupplementalTable2_diffabund.csv
            └── SupplementalTable3_NetworkProperties.xlsx
```
