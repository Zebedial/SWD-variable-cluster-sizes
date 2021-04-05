This repository includes R code for reproducing Figure 1, 2, 3, 4, 5, and Table 2, 3 in the article "Impact of unequal cluster sizes for GEE analyses of stepped wedge cluster randomized trials with binary outcomes" (under review), and Web Figure 1-50, and Web Table 1-23 in its supplementary material.

For questions or comments about the code, please contact Zibo Tian at <zibo.tian@yale.edu>.

I. Supporting Files: These supporting files are sourced in the main files that reproduce the simulation results in the manuscript.

1) NEX_functions.R = basic functions used to calculate variance of treatment effect estimator when the true correlation structure is nested exchangeable;

2) ED_functions.R = basic functions used to calculate variance of treatment effect estimator when the true correlation structure is exponential decay;

II. Main Files (Part 1): These main files are used to reproduce the simulation results and store the data into appropriate tables. Please make sure to change the directories in the files of interest in order to specify the path for saving the simulation results (.csv file).

3) RECV_nex.R = reproduce simulation results used to plot Figures like Figure 2 in the main article, which illustrate the relationship between relative efficiency (RE) and coefficient of variation (CV), when the true correlation structure is nested exhangeable.  Data used to reproduce Figure 2 and Figure 3 in the main article can be reproduced by runing related chunks of code in this file.

4) RECV_ed.R = reproduce simulation results used to plot RE vs CV figures when the true correlation structure is exponential decay.   

5) ICC_nex.R = reproduce simulation results used to plot Figures like Figure 4 in the main article, which illustrate the relationship between RE and alpha1/alpha0 or the impact of intraclass correlation coefficients (ICCs), when the true correlation structure is nested exhangeable and no within-cluster imbalance is introduced.  Data used to reproduce Figure 4 and Figure 5 in the main article can be reproduced by runing related chunks of code in this file.

6) ICC_nex_var.R = reproduce simulation results used to plot median RE vs alpha1/alpha0 figures when the true correlation structure is nested exhangeable and different patterns of within-cluster imbalance are introduced.  

7) ICC_ed.R = reproduce simulation results used to plot median RE vs rho figures when the true correlation structure is exponential decay.  

8) ICC_nex_cmean.R = reproduce simulation results used to plot additional median RE vs alpha1/alpha0 figures mentioned in Section 5.5 (illustrating the impact of the mean cluster-period size) when the true correlation structure is nested exhangeable.  

9) ICC_ed_cmean.R = reproduce simulation results used to plot additional median RE vs rho figures which are the counterparts to the figures mentioned in Section 5.5 (illustrating the impact of the mean cluster-period size) when the true correlation structure is exponential decay. 

10) Table_NEX_period.R = reproduce simulation results used to construct tables like Table 2 in the main article which illustrate the impact of number of period on RE, when the true correlation structure is nested exhangeable. Table 2 in the main article can be reproduced by runing related chunks of code in this file.

11) Table_ED_period.R = reproduce simulation results used to construct counterparts to tables like Table 2 in the main article when the true correlation structure is exponential decay. 

12) Table_NEX_nuisance_fixed.R = reproduce simulation results used to construct tables mentioned in Section 5.6 (illustrating the impact of treatment effect, prevalence trend, and baseline prevalence on RE), when the true correlation structure is nested exhangeable and no within-cluster imbalance is introduced.

13) Table_NEX_nuisance_pattern4.R = reproduce simulation results used to construct tables mentioned in Section 5.6 when the true correlation structure is nested exhangeable and pattern 4 (randomly permuted) within-cluster imbalance is introduced.

14) Table_ED_nuisance_fixed.R = reproduce simulation results used to construct counterparts to the tables produced in "Table_NEX_nuisance_fixed.R" when the true correlation structure is exponential decay.

15) Table_ED_nuisance_pattern4.R = reproduce simulation results used to construct counterparts to the tables produced in "Table_NEX_nuisance_pattern4.R" when the true correlation structure is exponential decay.

16) SWD_MC_Sample_Size_ttest.R = reproduce Table 3 in the main article that summarizes the application of sample size calculation algorithm developed in Section 7 of the main article.

17) EQcompare_nex.R = reproduce simulation results used to construct figures mentioned in Section 8 (illustrating the relationship between newly defined RE and alpha1/alpha0 under equal cluster sizes), when the true correlation structure is nested exhangeable.

18) EQcompare_ed.R = reproduce simulation results used to construct figures mentioned in Section 8 (illustrating the relationship between newly defined RE and rho under equal cluster sizes), when the true correlation structure is exponential decay.

III. Main Files (Part 2): These main files are used to reproduce the figures showed in the main article or the supplementary materials. Please make sure to change the directories in the files of interest to the location that stores the data sets obtained from Part 1 above. Please also make sure to change the directories in the files of interest in order to specify the path for saving the output figures.

19) Trial diagram_WEPT.R = reproduce Figure 1 in the main article.

20) Fig_REvsCV_NEX_fixed.R = reproduce Figure 2 in the main article and Web Figure 12. Need data produced by "RECV_nex.R".

21) Fig_REvsCV_NEX_pattern1.R = reproduce Web Figure 1. Need data produced by "RECV_nex.R".

22) Fig_REvsCV_NEX_pattern2.R = reproduce Web Figure 2. Need data produced by "RECV_nex.R".

23) Fig_REvsCV_NEX_pattern3.R = reproduce Web Figure 3. Need data produced by "RECV_nex.R".

24) Fig_REvsCV_NEX_pattern4.R = reproduce Figure 3 in the main article. Need data produced by "RECV_nex.R".

25) Fig_ICC_nex_fixed.R = reproduce Figure 4 and 5 in the main article. Need data produced by "ICC_nex.R".

26) Fig_ICC_nex_pattern1.R = reproduce Web Figure 4 and 8. Need data produced by "ICC_nex_var.R".

27) Fig_ICC_nex_pattern2.R = reproduce Web Figure 5 and 9. Need data produced by "ICC_nex_var.R".

28) Fig_ICC_nex_pattern3.R = reproduce Web Figure 6 and 10. Need data produced by "ICC_nex_var.R".

29) Fig_ICC_nex_pattern4.R = reproduce Web Figure 7 and 11. Need data produced by "ICC_nex_var.R".

30) Fig_cmean_nex_fixed.R = reproduce Web Figure 13, 14, 17, and 18. Need data produced by "ICC_nex_cmean.R".

31) Fig_cmean_nex_pattern4.R = reproduce Web Figure 15, 16, 19, and 20. Need data produced by "ICC_nex_cmean.R".

32) Fig_REvsCV_ED_fixed.R = reproduce Web Figure 21 and 22. Need data produced by "RECV_ed.R".

33) Fig_REvsCV_ED_pattern1.R = reproduce Web Figure 23. Need data produced by "RECV_ed.R".

34) Fig_REvsCV_ED_pattern2.R = reproduce Web Figure 24. Need data produced by "RECV_ed.R".

35) Fig_REvsCV_ED_pattern3.R = reproduce Web Figure 25. Need data produced by "RECV_ed.R".

36) Fig_REvsCV_ED_pattern4.R = reproduce Web Figure 26. Need data produced by "RECV_ed.R".

37) Fig_ICC_ed_fixed.R = reproduce Web Figure 27 and 32 in the main article. Need data produced by "ICC_ed.R".

38) Fig_ICC_ed_pattern1.R = reproduce Web Figure 28 and 33. Need data produced by "ICC_ed.R".

39) Fig_ICC_ed_pattern2.R = reproduce Web Figure 29 and 34. Need data produced by "ICC_ed.R".

40) Fig_ICC_ed_pattern3.R = reproduce Web Figure 30 and 35. Need data produced by "ICC_ed.R".

41) Fig_ICC_ed_pattern4.R = reproduce Web Figure 31 and 36. Need data produced by "ICC_ed.R".

42) Fig_cmean_ed_fixed.R = reproduce Web Figure 37, 38, 41, and 42. Need data produced by "ICC_ed_cmean.R".

43) Fig_cmean_ed_pattern4.R = reproduce Web Figure 39, 40, 43, and 44. Need data produced by "ICC_ed_cmean.R".

44) Fig_EQcompare_NEX.R = reproduce Web Figure 45. Need data produced by "EQcompare_nex.R".

45) Fig_EQcompare_ED.R = reproduce Web Figure 46. Need data produced by "EQcompare_ed.R".

46) Fig_REvsCV_NEX_fixed_Gbound.R = reproduce Web Figure 47 and 48. Need data produced by "RECV_nex.R".

47) Fig_REvsCV_ED_fixed_Gbound.R = reproduce Web Figure 49 and 50. Need data produced by "RECV_ed.R".


IV. Software 

Analyses were conducted with R, version 3.6.2 (https://www.r-project.org/)

V. R commands for the installation of R packages (Producing Figures)

install.packages(c("ggplot2", "ggpubr")) 