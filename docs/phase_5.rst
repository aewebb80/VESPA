**************************************
Phase 5: selection analysis assessment
**************************************

CodeML results assessment function
==================================
The ‘codeml_reader’ function is designed to parse the complex codeML directory structure and create simplified results for inexperienced users. This is achieved by incorporating in-house software ‘CreateSummaryReport.pl’ written by Dr. Thomas Walsh [Walsh, 2013] to produce the majority of the codeML results. In addition to automating ‘CreateSummaryReport.pl’, ‘codeml_reader’ produces supplementary output files (Figure 12) and specialized MSAs that are designed to aid in the detection of false positives (Figure 13). If the user specifies a branch-label table (Section 2.9.6) ‘codeml_reader’ will produce codeML MSAs, these MSAs are characterized by the addition of i) the putative positively selected sites, and ii) the codons/amino acids that are positively selected in the respective lineage/s.
usr$ python vespa.py codeml_reader –input=USR_INPUT
Supported file format(s): ‘input’: VESPA formatted codeML standard output.

Figure 12: Sample supplementary output file created by ‘codeml_reader’

Figure 12 Legend
The supplementary output file includes information for each site-specific and branch-specific model of codeML. The following information is provided for each model: the tree tested; the type of model (i.e. site-specific or branch-specific) being tested; number of free parameters in the ω distribution that are estimated by codeML, the initial ω value used by codeML; the resulting log likelihood (lnL) of the analysis; the resulting model of the likelihood ratio test (LRT); the parameter estimates of codeML; if positive selection was detected; and the positively selected sites (if positive selection was detected).

Figure 13: Sample specialized MSA created by ‘codeml_reader’

Figure 13 Legend
The specialized MSA shown above includes data on the location of positively selected codons or residues. Depending on the type of model being explored, the MSA will include additional information. For all models (site-specific or branch-specific), the header ‘PS_Sites’ indicates the position of the positively selected codons (shown as NNN) or residues (shown as X). For branch-specific, the characters under positive selection are shown for each relevant lineage using the header ‘PS_Characters’ followed the by the lineage of interest (i.e. PS_Characters|Chimp above).

