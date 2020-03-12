Untitled
================
Kendall Lin
3/5/2020

``` r
library(bio3d)
```

    ## Warning: package 'bio3d' was built under R version 3.6.3

``` r
aligned_file <- read.fasta("Icandoit")
aligned_matrix <- seqidentity(aligned_file, normalize=TRUE, similarity=FALSE, ncore=1, nseg.scale=1)
png(filename="heatmap.png", width = 700, height = 700, pointsize = 14)
heatmap(aligned_matrix, margins= c(12,12))
```

``` r
matrixsum <- rowSums(aligned_matrix)
matrixsum
```

    ##                 Platypus           Tammar_Wallaby                   Wombat 
    ##                    8.926                    8.595                    9.426 
    ##   Long-tailed_Chinchilla          Florida_Manatee African_Savanna_Elephant 
    ##                    9.132                    9.885                   10.192 
    ##           Golden_Hamster                    Human      Great_Roundleaf_Bat 
    ##                    9.710                   10.106                    9.914 
    ##      African_Clawed_Frog               Coelacanth       Tibetan_Ground_Tit 
    ##                    8.861                    8.956                    9.751 
    ##              Green_Anole   Central_Bearded_Dragon 
    ##                    9.047                    9.375

``` r
ok<- consensus(aligned_file)
boots <- toString(ok$seq)
toast<- gsub(", ", "", boots)
toast
```

    ## [1] "MALW---LPLLALLA---P----A-VNQHLCGSHLVEALYLVCGERGFFY-PK-RR--E-P-V-------G--------L-----E----KRGIV-QCC---CSLYQLENYCN"

``` r
#eleseq
```

``` r
eleseq <- c("MALWTRLLPLLALLAVGAPPPARAFVNQHLCGSHLVEALYLVCGERGFFYTPKTRREVEDTQVGEVELGTGLQPFPAEAPKQKRGIVEQCCTGVCSLYQLENYCN")
eleblast <- blast.pdb(eleseq, database = "pdb", time.out = NULL, chain.single=TRUE)
```

    ##  Searching ... please wait (updates every 5 seconds) RID = 6NHHG8UA016 
    ##  .
    ##  Reporting 100 hits

``` r
selectdata<- eleblast$raw$subjectid
annotate <- pdb.annotate(selectdata, anno.terms = c("structureId", "experimentalTechnique", "resolution", "source"))
```

    ## Warning in pdb.annotate(selectdata, anno.terms = c("structureId",
    ## "experimentalTechnique", : ids should be standard 4 character PDB-IDs: trying
    ## first 4 characters...

``` r
Evalue <- eleblast$hit.tbl$evalue
Identity <- eleblast$hit.tbl$identity
been <- cbind(annotate, Evalue , Identity )
eleblast$hit.tbl
```

    ##         queryid subjectids identity alignmentlength mismatches gapopens q.start
    ## 1   Query_56859     6B3Q_a   79.091             110         18        1       1
    ## 2   Query_56859     2KQP_A   76.744              86         15        1      25
    ## 3   Query_56859     6PXV_D   76.543              81         12        1      25
    ## 4   Query_56859     1EFE_A   66.667              81          6        2      25
    ## 5   Query_56859     1ZEI_A   58.025              81          6        1      25
    ## 6   Query_56859     5WDM_A   60.494              81          8        2      25
    ## 7   Query_56859     6INS_E   59.259              81          2        1      25
    ## 8   Query_56859     5WBT_A   59.259              81          9        2      25
    ## 9   Query_56859     1SJU_A   56.790              81          4        1      25
    ## 10  Query_56859     2JZQ_A   59.259              81          9        2      25
    ## 11  Query_56859     2KXK_B   96.875              32          1        0      25
    ## 12  Query_56859     2LGB_B  100.000              31          0        0      25
    ## 13  Query_56859     5MHD_B   96.774              31          1        0      25
    ## 14  Query_56859     1AI0_B  100.000              30          0        0      25
    ## 15  Query_56859     1IZA_B   96.667              30          1        0      25
    ## 16  Query_56859     4UNE_B   96.667              30          1        0      25
    ## 17  Query_56859     2MVD_B   96.667              30          1        0      25
    ## 18  Query_56859     1APH_B   96.667              30          1        0      25
    ## 19  Query_56859     1HLS_B   96.667              30          1        0      25
    ## 20  Query_56859     1ZEG_B   96.667              30          1        0      25
    ## 21  Query_56859     1B9E_B   96.667              30          1        0      25
    ## 22  Query_56859     1UZ9_B  100.000              29          0        0      25
    ## 23  Query_56859     1QIY_B   96.667              30          1        0      25
    ## 24  Query_56859     3V1G_B   96.667              30          1        0      25
    ## 25  Query_56859     5HPR_B   96.667              30          1        0      25
    ## 26  Query_56859     3JSD_B   96.667              30          1        0      25
    ## 27  Query_56859     1W8P_B   93.333              30          2        0      25
    ## 28  Query_56859     2M2M_B   96.667              30          1        0      25
    ## 29  Query_56859     3ZQR_B   96.667              30          1        0      25
    ## 30  Query_56859     2M2O_B   96.667              30          1        0      25
    ## 31  Query_56859     1HIQ_B   96.667              30          1        0      25
    ## 32  Query_56859     2WRX_B   96.667              30          1        0      25
    ## 33  Query_56859     3BXQ_B   96.667              30          1        0      25
    ## 34  Query_56859     4UNG_B   96.667              30          1        0      25
    ## 35  Query_56859     1HIT_B   96.667              30          1        0      25
    ## 36  Query_56859     1MHI_B   96.552              29          1        0      26
    ## 37  Query_56859     4UNH_B   96.667              30          1        0      25
    ## 38  Query_56859     6TYH_B   96.667              30          1        0      25
    ## 39  Query_56859     6CK2_B   93.333              30          2        0      25
    ## 40  Query_56859     2INS_B   96.552              29          1        0      26
    ## 41  Query_56859     6GV0_B   93.333              30          2        0      25
    ## 42  Query_56859     4A7E_B   93.333              30          2        0      25
    ## 43  Query_56859     1QIY_D   96.552              29          1        0      26
    ## 44  Query_56859     5AZZ_B   93.333              30          2        0      25
    ## 45  Query_56859     5BPO_B   93.333              30          2        0      25
    ## 46  Query_56859     4NIB_B   93.333              30          2        0      25
    ## 47  Query_56859     4P65_B   93.333              30          2        0      25
    ## 48  Query_56859     2N2V_B   93.333              30          2        0      25
    ## 49  Query_56859     1LPH_B   93.333              30          2        0      25
    ## 50  Query_56859     1JCO_B   93.333              30          2        0      25
    ## 51  Query_56859     3U4N_B   96.552              29          1        0      25
    ## 52  Query_56859     1HO0_A   93.103              29          2        0      25
    ## 53  Query_56859     5BQQ_B   96.429              28          1        0      25
    ## 54  Query_56859     1HTV_B  100.000              27          0        0      25
    ## 55  Query_56859     1SJT_B   93.103              29          2        0      25
    ## 56  Query_56859     1MHJ_B   96.667              30          0        1      25
    ## 57  Query_56859     4EFX_B   96.429              28          1        0      25
    ## 58  Query_56859     2N2X_B   90.000              30          3        0      25
    ## 59  Query_56859     2MPG_B   90.000              30          3        0      25
    ## 60  Query_56859     1K3M_B   90.000              30          3        0      25
    ## 61  Query_56859     1A7F_B   93.103              29          2        0      25
    ## 62  Query_56859     3ZS2_B   86.667              30          4        0      25
    ## 63  Query_56859     1T1P_B   86.667              30          4        0      25
    ## 64  Query_56859     1T1K_B   86.667              30          4        0      25
    ## 65  Query_56859     2HHO_B   86.667              30          4        0      25
    ## 66  Query_56859     1T1Q_B   86.667              30          4        0      25
    ## 67  Query_56859     2KJU_B   86.667              30          4        0      25
    ## 68  Query_56859     2L1Y_B   86.667              30          4        0      25
    ## 69  Query_56859     1HIS_B  100.000              25          0        0      25
    ## 70  Query_56859     1BZV_B  100.000              25          0        0      25
    ## 71  Query_56859     2HH4_B   86.667              30          4        0      25
    ## 72  Query_56859     2K9R_B   86.667              30          4        0      25
    ## 73  Query_56859     2MLI_B   86.667              30          4        0      25
    ## 74  Query_56859     2H67_B   86.667              30          4        0      25
    ## 75  Query_56859     2K91_B   86.667              30          4        0      25
    ## 76  Query_56859     2WS4_B  100.000              25          0        0      25
    ## 77  Query_56859     2MPI_B   86.667              30          4        0      25
    ## 78  Query_56859     1HUI_B   89.286              28          3        0      26
    ## 79  Query_56859     1SDB_B  100.000              23          0        0      27
    ## 80  Query_56859     1DEI_B  100.000              23          0        0      25
    ## 81  Query_56859     1IGL_A   33.333              78         30        1      27
    ## 82  Query_56859     6UM2_B   33.333              78         30        1      27
    ## 83  Query_56859     3LRI_A   33.708              89         39        2      18
    ## 84  Query_56859     2QIU_A   90.909              22          2        0      84
    ## 85  Query_56859     1BQT_A   36.364              77         31        2      28
    ## 86  Query_56859     5L3M_A   35.443              79         29        2      27
    ## 87  Query_56859     6RVA_X   35.065              77         32        2      28
    ## 88  Query_56859     2BN3_A   95.238              21          1        0      85
    ## 89  Query_56859     5L3N_A   35.443              79         29        2      27
    ## 90  Query_56859     5MHD_A   90.476              21          2        0      85
    ## 91  Query_56859     2KXK_A   90.476              21          2        0      85
    ## 92  Query_56859     1A7F_A   90.476              21          2        0      85
    ## 93  Query_56859     1TGR_A   35.065              77         22        2      28
    ## 94  Query_56859     1APH_A   90.476              21          2        0      85
    ## 95  Query_56859     6TYH_A   90.476              21          2        0      85
    ## 96  Query_56859     1XW7_A   85.714              21          3        0      85
    ## 97  Query_56859     2JUM_A   85.714              21          3        0      85
    ## 98  Query_56859     1K3M_A   85.714              21          3        0      85
    ## 99  Query_56859     2WBY_D  100.000              19          0        0      25
    ## 100 Query_56859     1JCA_A   85.714              21          3        0      85
    ##     q.end s.start s.end   evalue bitscore positives mlog.evalue pdb.id    acc
    ## 1     105       1   110 2.15e-50    155.0     81.82   114.36379 6B3Q_a 6B3Q_a
    ## 2     105       1    86 1.82e-41    132.0     79.07    93.80715 2KQP_A 2KQP_A
    ## 3     105       1    74 1.94e-39    126.0     77.78    89.13813 6PXV_D 6PXV_D
    ## 4     105       1    60 6.72e-31    104.0     70.37    69.47505 1EFE_A 1EFE_A
    ## 5     105       1    53 2.13e-25     90.5     60.49    56.80851 1ZEI_A 1ZEI_A
    ## 6     105       1    57 3.60e-25     90.1     64.20    56.28369 5WDM_A 5WDM_A
    ## 7     105       1    50 5.32e-25     89.7     60.49    55.89315 6INS_E 6INS_E
    ## 8     105       1    57 1.07e-23     86.7     61.73    52.89180 5WBT_A 5WBT_A
    ## 9     105       1    50 3.41e-23     85.1     58.02    51.73274 1SJU_A 1SJU_A
    ## 10    105       1    57 6.04e-23     84.7     61.73    51.16105 2JZQ_A 2JZQ_A
    ## 11     56       1    32 2.04e-17     70.1    100.00    38.43100 2KXK_B 2KXK_B
    ## 12     55       1    31 5.21e-17     68.9    100.00    37.49337 2LGB_B 2LGB_B
    ## 13     55       1    31 3.72e-16     66.6     96.77    35.52764 5MHD_B 5MHD_B
    ## 14     54       1    30 3.92e-16     66.6    100.00    35.47527 1AI0_B 1AI0_B
    ## 15     54       1    30 1.02e-15     65.5    100.00    34.51897 1IZA_B 1IZA_B
    ## 16     54       1    30 1.27e-15     65.5    100.00    34.29976 4UNE_B 4UNE_B
    ## 17     54       1    30 1.77e-15     65.1    100.00    33.96780 2MVD_B 2MVD_B
    ## 18     54       1    30 1.87e-15     65.1     96.67    33.91284 1APH_B 1APH_B
    ## 19     54       1    30 1.99e-15     65.1    100.00    33.85064 1HLS_B 1HLS_B
    ## 20     54       1    30 2.08e-15     64.7     96.67    33.80641 1ZEG_B 1ZEG_B
    ## 21     54       1    30 2.22e-15     64.7     96.67    33.74127 1B9E_B 1B9E_B
    ## 22     53       1    29 2.29e-15     64.7    100.00    33.71022 1UZ9_B 1UZ9_B
    ## 23     54       1    30 2.77e-15     64.7    100.00    33.51993 1QIY_B 1QIY_B
    ## 24     54       1    30 3.15e-15     64.3     96.67    33.39137 3V1G_B 3V1G_B
    ## 25     54       1    30 3.75e-15     64.3     96.67    33.21702 5HPR_B 5HPR_B
    ## 26     54       1    30 4.23e-15     63.9     96.67    33.09657 3JSD_B 3JSD_B
    ## 27     54       1    30 4.42e-15     63.9    100.00    33.05264 1W8P_B 1W8P_B
    ## 28     54       1    30 4.52e-15     63.9     96.67    33.03026 2M2M_B 2M2M_B
    ## 29     54       1    30 5.04e-15     63.9     96.67    32.92137 3ZQR_B 3ZQR_B
    ## 30     54       1    30 5.04e-15     63.9     96.67    32.92137 2M2O_B 2M2O_B
    ## 31     54       1    30 5.44e-15     63.9     96.67    32.84500 1HIQ_B 1HIQ_B
    ## 32     54       1    30 5.50e-15     63.9     96.67    32.83403 2WRX_B 2WRX_B
    ## 33     54       1    30 5.62e-15     63.9     96.67    32.81244 3BXQ_B 3BXQ_B
    ## 34     54       1    30 6.40e-15     63.5     96.67    32.68248 4UNG_B 4UNG_B
    ## 35     54       1    30 1.01e-14     63.2     96.67    32.22624 1HIT_B 1HIT_B
    ## 36     54       2    30 1.05e-14     63.2     96.55    32.18740 1MHI_B 1MHI_B
    ## 37     54       1    30 1.17e-14     63.2     96.67    32.07919 4UNH_B 4UNH_B
    ## 38     54       1    30 1.19e-14     62.8     96.67    32.06224 6TYH_B 6TYH_B
    ## 39     54       1    30 1.22e-14     62.8     96.67    32.03734 6CK2_B 6CK2_B
    ## 40     54       1    29 1.26e-14     62.8     96.55    32.00508 2INS_B 2INS_B
    ## 41     54       1    30 1.45e-14     62.8     96.67    31.86463 6GV0_B 6GV0_B
    ## 42     54       1    30 1.57e-14     62.8     93.33    31.78512 4A7E_B 4A7E_B
    ## 43     54       2    30 1.93e-14     62.4    100.00    31.57867 1QIY_D 1QIY_D
    ## 44     54       1    30 2.37e-14     62.4     93.33    31.37330 5AZZ_B 5AZZ_B
    ## 45     54       1    30 2.40e-14     62.0     93.33    31.36072 5BPO_B 5BPO_B
    ## 46     54       1    30 3.87e-14     61.6     93.33    30.88294 4NIB_B 4NIB_B
    ## 47     54       1    30 4.36e-14     61.6     93.33    30.76372 4P65_B 4P65_B
    ## 48     54       1    30 5.08e-14     61.2     93.33    30.61088 2N2V_B 2N2V_B
    ## 49     54       1    30 5.67e-14     61.2     93.33    30.50100 1LPH_B 1LPH_B
    ## 50     54       1    30 6.05e-14     61.2     93.33    30.43613 1JCO_B 1JCO_B
    ## 51     53       1    29 6.60e-14     61.2     96.55    30.34912 3U4N_B 3U4N_B
    ## 52     53       1    29 7.28e-14     60.8     93.10    30.25106 1HO0_A 1HO0_A
    ## 53     52       1    28 1.37e-13     60.1     96.43    29.61880 5BQQ_B 5BQQ_B
    ## 54     51       1    27 1.42e-13     60.1    100.00    29.58295 1HTV_B 1HTV_B
    ## 55     53       1    29 1.63e-13     60.1     93.10    29.44503 1SJT_B 1SJT_B
    ## 56     54       1    29 1.80e-13     60.1     96.67    29.34582 1MHJ_B 1MHJ_B
    ## 57     52       1    28 3.81e-13     58.9     96.43    28.59598 4EFX_B 4EFX_B
    ## 58     54       1    30 6.09e-13     58.5     90.00    28.12696 2N2X_B 2N2X_B
    ## 59     54       1    30 6.30e-13     58.5     90.00    28.09306 2MPG_B 2MPG_B
    ## 60     54       1    30 7.25e-13     58.5     90.00    27.95260 1K3M_B 1K3M_B
    ## 61     53       1    29 1.06e-12     58.2     93.10    27.57275 1A7F_B 1A7F_B
    ## 62     54       1    30 2.38e-12     57.0     90.00    26.76392 3ZS2_B 3ZS2_B
    ## 63     54       1    30 2.80e-12     57.0     86.67    26.60140 1T1P_B 1T1P_B
    ## 64     54       1    30 2.89e-12     57.0     86.67    26.56976 1T1K_B 1T1K_B
    ## 65     54       1    30 3.84e-12     56.6     86.67    26.28555 2HHO_B 2HHO_B
    ## 66     54       1    30 4.18e-12     56.6     86.67    26.20071 1T1Q_B 1T1Q_B
    ## 67     54       1    30 5.80e-12     56.2     86.67    25.87316 2KJU_B 2KJU_B
    ## 68     54       1    30 6.06e-12     56.2     86.67    25.82931 2L1Y_B 2L1Y_B
    ## 69     49       1    25 7.82e-12     55.8    100.00    25.57434 1HIS_B 1HIS_B
    ## 70     49       1    25 7.88e-12     55.8    100.00    25.56669 1BZV_B 1BZV_B
    ## 71     54       1    30 8.04e-12     55.8     86.67    25.54659 2HH4_B 2HH4_B
    ## 72     54       1    30 9.26e-12     55.5     86.67    25.40532 2K9R_B 2K9R_B
    ## 73     54       1    30 9.26e-12     55.5     86.67    25.40532 2MLI_B 2MLI_B
    ## 74     54       1    30 1.07e-11     55.5     86.67    25.26078 2H67_B 2H67_B
    ## 75     54       1    30 1.13e-11     55.5     86.67    25.20622 2K91_B 2K91_B
    ## 76     49       1    25 1.13e-11     55.5    100.00    25.20622 2WS4_B 2WS4_B
    ## 77     54       1    30 1.78e-11     54.7     86.67    24.75182 2MPI_B 2MPI_B
    ## 78     53       2    29 2.26e-11     54.7     89.29    24.51307 1HUI_B 1HUI_B
    ## 79     49       1    23 2.37e-10     52.0    100.00    22.16296 1SDB_B 1SDB_B
    ## 80     47       1    23 4.59e-10     51.2    100.00    21.50197 1DEI_B 1DEI_B
    ## 81    104       5    60 1.06e-08     48.5     46.15    18.36241 1IGL_A 1IGL_A
    ## 82    104       6    61 1.10e-08     48.5     46.15    18.32537 6UM2_B 6UM2_B
    ## 83    104       4    74 1.98e-08     48.1     47.19    17.73758 3LRI_A 3LRI_A
    ## 84    105       1    22 2.96e-08     46.6     95.45    17.33549 2QIU_A 2QIU_A
    ## 85    104       3    61 3.10e-08     47.4     49.35    17.28928 1BQT_A 1BQT_A
    ## 86    104       5    62 4.57e-08     47.0     46.84    16.90117 5L3M_A 5L3M_A
    ## 87    104       4    62 6.07e-08     46.6     49.35    16.61732 6RVA_X 6RVA_X
    ## 88    105       1    21 7.73e-08     45.4     95.24    16.37557 2BN3_A 2BN3_A
    ## 89    104       5    62 8.10e-08     46.2     45.57    16.32882 5L3N_A 5L3N_A
    ## 90    105       1    21 1.51e-07     44.7     95.24    15.70599 5MHD_A 5MHD_A
    ## 91    105       1    21 1.80e-07     44.7     95.24    15.53031 2KXK_A 2KXK_A
    ## 92    105       1    21 2.11e-07     44.3     95.24    15.37141 1A7F_A 1A7F_A
    ## 93    104       3    51 3.05e-07     44.7     46.75    15.00295 1TGR_A 1TGR_A
    ## 94    105       1    21 4.03e-07     43.5     90.48    14.72433 1APH_A 1APH_A
    ## 95    105       1    21 5.06e-07     43.5     95.24    14.49673 6TYH_A 6TYH_A
    ## 96    105       1    21 7.05e-07     43.1     95.24    14.16507 1XW7_A 1XW7_A
    ## 97    105       1    21 7.85e-07     42.7     90.48    14.05758 2JUM_A 2JUM_A
    ## 98    105       1    21 7.94e-07     42.7     90.48    14.04618 1K3M_A 1K3M_A
    ## 99     43       1    19 9.15e-07     42.7    100.00    13.90434 2WBY_D 2WBY_D
    ## 100   105       1    21 9.33e-07     42.7     90.48    13.88486 1JCA_A 1JCA_A

``` r
chosen_id <- c("6B3Q_a", "1ZEI_A", "1APH_B")
annotate2 <- pdb.annotate(chosen_id, anno.terms = c("structureId", "experimentalTechnique", "resolution", "source" ))
```

    ## Warning in pdb.annotate(chosen_id, anno.terms = c("structureId",
    ## "experimentalTechnique", : ids should be standard 4 character PDB-IDs: trying
    ## first 4 characters...

``` r
cleantable <- been[chosen_id,]
cleantable
```

    ##        structureId experimentalTechnique resolution       source   Evalue
    ## 6B3Q_a        6B3Q   ELECTRON MICROSCOPY        3.7 Homo sapiens 2.15e-50
    ## 1ZEI_A        1ZEI     X-RAY DIFFRACTION        1.9   Sus scrofa 2.13e-25
    ## 1APH_B        1APH     X-RAY DIFFRACTION        2.0   Bos taurus 1.87e-15
    ##        Identity
    ## 6B3Q_a   79.091
    ## 1ZEI_A   58.025
    ## 1APH_B   96.667

``` r
file.name <- get.pdb("1aph_B")
```

    ## Warning in get.pdb("1aph_B"): ./1aph.pdb exists. Skipping download

``` r
insulin <- read.pdb(file.name)
```

    ##    PDB has ALT records, taking A only, rm.alt=TRUE

``` r
prot <- trim.pdb(insulin, "protein")
lig <- trim.pdb(insulin, "ligand")
write.pdb(prot, file="insulin_protein.pdb")
write.pdb(lig, file="insulin_ligand.pdb")
```
