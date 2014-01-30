The following files are the supplementary files of "Inferring
Time-Delayed Causal Gene Network using Time-series Expression Data" by
Leung-Yau Lo (lylo@cse.cuhk.edu.hk), Prof. Kwong-Sak Leung
(ksleung@cse.cuhk.edu.hk) and Prof. Kin-Hong Lee
(khlee@cse.cuhk.edu.hk) from the Chinese University of Hong Kong.

clinde_supp.ps:
	The supplementary figures and tables in Postscript format.
clinde_supp.pdf:
	The supplementary figures and tables in PDF format.

programs/
	infer.R:
		Contains the main function (infer_grn) for CLINDE.
	pcor.R:
		By Kim, S-H. and Yi, S. (2007) Understanding relationship between sequence and functional evolution in yeast proteins . Genetica, In press. Downloaded from http://www.yilab.gatech.edu/pcor.html. Functions to calculate partial correlation.
	sim.R:
		functions for simulating the synthetic networks and gene expression data.
	test.R:
		when provided with the values of test.n, test.p1, test.p2, test.m, test.max.n, test.max.parents, test.max.delay, test.n.points, test.is.acyclic, test.is.self, test.method and test.pruning, will randomly generate synthetic network, gene expression and infer the gene network from the generated expression, and calculate the accuracy of the inferred network.
	test_real.R:
		provides the function (test1) to test the inference performance on real expression data. When the expression data has only one segment, segs should be set to c(L), where L is the number of time points.
	real_data.RData:
		R data file containing the relevant expression data and true network of the IRMA dataset.
	exp_data/kegg_meiosis_exp.RData:
		R data file containing the relevant expression data and true network of the KEGG (cell-cycle and meiosis) network dataset.
	kegg_meiosis_sn.RData
		R data file containing the relevant expression data and true network of the subnetworks of the KEGG (cell-cycle and meiosis) network dataset.

results/
	test_grn_res_2012_m_01_d_01_h_23_m_39_pcor_medians.csv:
	test_grn_res_2012_m_01_d_03_h_10_m_05_pcor_medians.csv:
	test_grn_res_2012_m_01_d_03_h_10_m_07_pcor_medians.csv:
	test_grn_res_2012_m_01_d_03_h_10_m_06_pcor_medians.csv:
		The above 4 files contain the results on the synthetic data for various parameter settings.
	test_grn_res_2012_m_01_d_02_h_01_m_37_pcor_real.csv:
		Contains the results on the IRMA dataset.
	test_grn_res_2012_m_04_d_30_h_00_m_22_pcor_real_sn.csv:
		Contains the results on the KEGG subnetworks.

comparison.rar:
	contains the testing shell scripts, R scripts for testing our algorithm CLINDE with TD-ARACNE and Banjo on synthetic data.
	also contains the testing scripts for segmented data on our algorithm CLINDE.

