# Semi-Automation of the Ancestral Sequence Reconstruction Workflow
# This is a program that constructs rough estimates of ancestral sequence reconstruction from an amino acid sequence.
# This program requires instalation of the Bioservieces module
# Written by James VanAntwerp from September 2020 through May 2023 - contact  vanantj @ udel . edu
# Written by Pattrick Finneran, Menten AI, Palo Alto, California, United States of America
# Written for the Woldring Lab, Michigan State University in East Lansing,
# Michigan, USA.

import sys
from Bio.Blast import NCBIWWW
import xml.etree.ElementTree as ET
import os
import argparse

Directory_Structure = '''The program will make the following directory structure, where ROOT is the directory the script is placed in:
    ROOT/
    |AP-LASR.py                               This script
    |--Input.fasta                            The input fasta file, if you used one.
    |--ASR/...................................The directory containing all work from the most recent run. The software WILL overwrite old runs.
    | |--BlastP_Results.fasta                       Fasta file containing all BlastP results.
    | |--BlastP_XML_                                The raw XML returned by NCBI's BlastP, if you want to find more information about the BlastP resutls.
    | |--HitsInfo.csv                               Data about the BlastP hits - Sequence IDs and notes.
    | |--Final_Sequences.fasta                      The set of modern sequences used for the ASR and the foundation of all the later calculation.
    | |--Final_Tree.treefile                        The final tree which is best for reading, as it has the most information in one place.
    | |--Concesus_Ancestors_with_Gaps.fasta         The set of likely ancestors, aligned and with gaps where the binary gap analysis predicts them.
    | |--Concesus_Ancestors_without_Gaps.fasta      The set of likely ancestors, unaligned and without gaps.
    | |--IQTree_Phylo/..............................The directory for IQTree files from the phylogeny.
    | | |--Phylo.*                                      All IQTree files.
    | | |--Supports.*                                   Data about the confidence of tree topology - .csv holds all data, .txt is summary.
    | | |--UFB_Confidences.png                          Histogram of UFB node support values (if made).
    | | |--SHaLRT_Confidences.png                       Histogram of SHaLRT node support values (if made).
    | |--IQTree_ASR/................................The directory for IQTree files from ASR.
    | | |--ASR.*                                        All IQTree files.
    | | |--Ancestral_Sequence_Confidences.*             Data about the confidence of ASR prediction - .csv holds all data, .txt is summary, .png is a histogram (if made).
    | |--IQTree_Binary/.............................The directory for IQTree files from gap analysis.
    | | |--Binary_Alignment.fasta                       A binary alignment of modern sequences - gap or not.
    | | |--Binary.*                                     All IQTree files.
    | |--DNA_Libraries/.............................The directory for IQTree files from gap analysis.
    | | |--Library_Size_Information.csv                 Information on the size of generated libraries.
    | | |--Library_**%_Cutoff.fasta                     Fasta of DNA for all high-confidence ancestors, with degenerate codons coding for all AAs with more than a certian likelihood.
    | |--Sequence_Supplement/.......................The directory for the files created by adding additional sequence data to poorly supported regions of the tree.
    | | |--*__supplement_HitsInfo.csv                   Same as HitInfo.csv, but for each supplement.
    | | |--*_BlastP_Results.fasta                       Each top 50 hits from sequences that are supplemented.
    | |--Confidence_Heatmaps/.......................The directory for IQTree files from gap analysis. Only generated after running the script again with the MakeFigures option.
    | | |--*_Confidences.pdf                            A heatmap of the confidence values for each ancestral position.
    '''

Help_String = '''This software automates the process of generating combinatorial protein libraries from ancestral sequence reconstruction.
This software can run in one of four modes.
     Assistance mode will show this message, software prerequistes, and an example directory structure.
     ASR mode will take in a user sequence and conduct ASR, then generate combinatorial protein libraries of ancestral proteins. Command-line parameters should be used to specify the directory where output is stored and the executables for CD-Hit, MAFFT, and IQ-Tree. The desired final dataset size can be specified, and sequence supplementation can be turned on here too.
        option -i specifies the name of the input Fasta file, or can take a raw protein sequence. This is the only mandatory input.
        option -n specifies the name of the output directory. The default is "ASR".
        option -s specifies the desired maximum size of the final dataset of modern homologs. This allows users some control over the time AP-LASR will take to run, and the detail/quality of the ASR. The default is 500.
        option -supplement will turn on supplementation at the provided similarity cutoff (entered as a decimal). If you do turn on supplementation, there is not a default value but we reccomend users specify something between 0.7 and 0.85. 
        option -iqtree allows users to specify an executable for IQTree other than the default ("iqtree").
        option -cdhit allows users to specify an executable for CH-Hit other than the default ("cd-hit").
        option -mafft allows users to specify an executable for MAFFT other than the default ("mafft").
        option -MSI allows users to specify the number of times the function "Post_MAFFT_Processing" iterates over the raw alignment to generate the final dataset. Each iteation will remove a large number of sequences, more if the different input sequecnes are dissimilar. Default is 2.
     MakeFigures mode will generate figures from ASR data of a previous ASR run - specify the directory where the results of interest are stored.
        option -n specifies the name of the directory where the ASR results are stored. The default is "ASR".
     RemakeLibrareies mode will generate new combinatorial libraries from ancestral proteins of a previous ASR run with a specified threshold confidence. This will allow generation of a library with a different size than what was generated with the default threshold confidences. Specify the directory where the ASR results of interest are stored.
        option -n specifies the name of the directory where the ASR results are stored. The default is "ASR".
        option -t specifies the new threshold, expressed as a decimal. 
Entering just the option "--help" or "-h" will give a less verbose version this description. '''

Software_Prerequistes = '''This AP-LASR script makes use of external software, which must be set up before it can be used. See the GitHub page for this project (github.com/jjvanantwerp/Automated-ASR) for download links.
This software requires IQTree, CD-Hit, and MAFFT to be downloaded and set up with their standard executable names added to the path.
It also requires the python modules BioPython (Bio) and xml to be installed. Optionally, the python module matplotlib can be installed for production of figures from the data.
To make the pdf confidence heatmaps, the python modules pandas and seaborn are also needed.'''

# This dictionary provides the amino acid encoded for by every codon
Codon_to_AA = {
    'ata': 'I', 'atc': 'I', 'att': 'I', 'atg': 'M',
    'aca': 'T', 'acc': 'T', 'acg': 'T', 'act': 'T',
    'aac': 'N', 'aat': 'N', 'aaa': 'K', 'aag': 'K',
    'agc': 'S', 'agt': 'S', 'aga': 'R', 'agg': 'R',
    'cta': 'L', 'ctc': 'L', 'ctg': 'L', 'ctt': 'L',
    'cca': 'P', 'ccc': 'P', 'ccg': 'P', 'cct': 'P',
    'cac': 'H', 'cat': 'H', 'caa': 'Q', 'cag': 'Q',
    'cga': 'R', 'cgc': 'R', 'cgg': 'R', 'cgt': 'R',
    'gta': 'V', 'gtc': 'V', 'gtg': 'V', 'gtt': 'V',
    'gca': 'A', 'gcc': 'A', 'gcg': 'A', 'gct': 'A',
    'gac': 'D', 'gat': 'D', 'gaa': 'E', 'gag': 'E',
    'gga': 'G', 'ggc': 'G', 'ggg': 'G', 'ggt': 'G',
    'tca': 'S', 'tcc': 'S', 'tcg': 'S', 'tct': 'S',
    'ttc': 'F', 'ttt': 'F', 'tta': 'L', 'ttg': 'L',
    'tac': 'Y', 'tat': 'Y', 'tgc': 'C', 'tgt': 'C',
    'taa': ' stop ', 'tag': ' stop ', 'tga': ' stop ', 'tgg': 'W'
}
# This dictionary provides the best codon for every amino acid in E Coli
AA_to_Codon_Ecoli = {
    'A': 'gcc', 'R': 'cgt', 'N': 'aac', 'D': 'gat', 'B': 'aac',
    'C': 'tgc', 'E': 'gaa', 'Q': 'cag', 'Z': 'cag', 'G': 'ggc',
    'H': 'cat', 'I': 'att', 'L': 'ctg', 'K': 'aaa', 'M': 'atg',
    'F': 'ttt', 'P': 'ccg', 'S': 'agc', 'T': 'acc', 'W': 'tgg',
    'Y': 'tat', 'V': 'gtg', ' stop ': 'taa'
}
# This dictionary provides the best codon for every amino acid in humans
AA_to_Codon_Human = {
    'A': 'gcc', 'R': 'cgg', 'N': 'aac', 'D': 'gac', 'B': 'aac',
    'C': 'tgc', 'E': 'gag', 'Q': 'cag', 'Z': 'cag', 'G': 'ggc',
    'H': 'cac', 'I': 'atc', 'L': 'ctg', 'K': 'aag', 'M': 'atg',
    'F': 'ttc', 'P': 'ccc', 'S': 'agc', 'T': 'acc', 'W': 'tgg',
    'Y': 'tac', 'V': 'gtc', ' stop ': 'tga'
}
# This dictionary provides the mixed-base code for every combination of 2-4 nucleic acids
Mixed_Bases_lookup = {
    'a':'a','c':'c','g':'g','t':'t',

    'ag':'r', 'ga':'r',
    'ct':'y', 'tc':'y',
    'ac':'m', 'ca':'m',
    'gt':'k', 'tg':'k',
    'gc':'s', 'cg':'s',
    'at':'w', 'ta':'w',

    'act':'h', 'atc':'h', 'cat':'h', 'cta':'h', 'tac':'h', 'tca':'h',
    'gct':'b', 'gtc':'b', 'ctg':'b', 'cgt':'b', 'tgc':'b', 'tcg':'b',
    'acg':'v', 'agc':'v', 'cag':'v', 'cga':'v', 'gca':'v', 'gac':'v',
    'agt':'d', 'atg':'d', 'gat':'d', 'gta':'d', 'tga':'d', 'tag':'d',

    'acgt':'n','actg':'n','agct':'n','agtc':'n','atcg':'n','atgc':'n',
    'cagt':'n','catg':'n','gact':'n','gatc':'n','tacg':'n','tagc':'n',
    'cgat':'n','ctag':'n','gcat':'n','gtac':'n','tcag':'n','tgac':'n',
    'cgta':'n','ctga':'n','gcta':'n','gtca':'n','tcga':'n','tgca':'n'
    }
# This dictionary provides the best degenerate codon for every combination of two amino acids for Humans
AA_Pair_lookup_Human = { 
    'AC':'ksc', 'AD':'gmc', 'AE':'gma', 'AF':'kyc', 'AG':'gsc', 'AH':'smc', 'AI':'ryc', 'AK':'rma', 'AL':'syc', 'AM':'ryg', 'AN':'rmc', 'AP':'scc', 'AQ':'sma', 'AR':'rsa', 'AS':'kcc', 'AT':'rcc', 'AV':'gyc', 'AW':'ksg', 'AY':'kmc', 
    'CA':'ksc', 'CD':'krc', 'CE':'krs', 'CF':'tkc', 'CG':'kgc', 'CH':'yrc', 'CI':'wkc', 'CK':'wrs', 'CL':'ykc', 'CM':'wks', 'CN':'wrc', 'CP':'ysc', 'CQ':'yrs', 'CR':'ygc', 'CS':'tsc', 'CT':'wsc', 'CV':'kkc', 'CW':'tgs', 'CY':'trc', 
    'DA':'gmc', 'DC':'krc', 'DE':'gas', 'DF':'kwc', 'DG':'grc', 'DH':'sac', 'DI':'rwc', 'DK':'ras', 'DL':'swc', 'DM':'rws', 'DN':'rac', 'DP':'smc', 'DQ':'sas', 'DR':'src', 'DS':'rrc', 'DT':'rmc', 'DV':'gwc', 'DW':'krs', 'DY':'kac', 
    'EA':'gma', 'EC':'krs', 'ED':'gas', 'EF':'kws', 'EG':'grg', 'EH':'sas', 'EI':'rwa', 'EK':'rag', 'EL':'swg', 'EM':'rwg', 'EN':'ras', 'EP':'smg', 'EQ':'sag', 'ER':'rrg', 'ES':'kmg', 'ET':'rmg', 'EV':'gwg', 'EW':'rrg', 'EY':'kas', 
    'FA':'kyc', 'FC':'tkc', 'FD':'kwc', 'FE':'kws', 'FG':'kkc', 'FH':'ywc', 'FI':'wtc', 'FK':'wwm', 'FL':'ytc', 'FM':'wts', 'FN':'wwc', 'FP':'yyc', 'FQ':'yws', 'FR':'ykc', 'FS':'tyc', 'FT':'wyc', 'FV':'ktc', 'FW':'tks', 'FY':'twc', 
    'GA':'gsc', 'GC':'kgc', 'GD':'grc', 'GE':'grg', 'GF':'kkc', 'GH':'src', 'GI':'rkc', 'GK':'rra', 'GL':'skc', 'GM':'rrg', 'GN':'rrc', 'GP':'ssc', 'GQ':'grg', 'GR':'sgg', 'GS':'rgc', 'GT':'rsc', 'GV':'gkc', 'GW':'kgg', 'GY':'krc', 
    'HA':'smc', 'HC':'yrc', 'HD':'sac', 'HE':'sas', 'HF':'ywc', 'HG':'src', 'HI':'mwc', 'HK':'mas', 'HL':'cwc', 'HM':'mws', 'HN':'mac', 'HP':'cmc', 'HQ':'cas', 'HR':'crc', 'HS':'mrc', 'HT':'mmc', 'HV':'swc', 'HW':'yrs', 'HY':'yac', 
    'IA':'ryc', 'IC':'wkc', 'ID':'rwc', 'IE':'rwa', 'IF':'wtc', 'IG':'rkc', 'IH':'mwc', 'IK':'awa', 'IL':'mtc', 'IM':'ats', 'IN':'awc', 'IP':'myc', 'IQ':'mya', 'IR':'aka', 'IS':'akc', 'IT':'ayc', 'IV':'rtc', 'IW':'wks', 'IY':'wwc', 
    'KA':'rma', 'KC':'wrs', 'KD':'ras', 'KE':'rag', 'KF':'wwm', 'KG':'rra', 'KH':'mas', 'KI':'awa', 'KL':'mwa', 'KM':'awg', 'KN':'aas', 'KP':'mma', 'KQ':'maa', 'KR':'arg', 'KS':'ars', 'KT':'ama', 'KV':'rwa', 'KW':'wrg', 'KY':'was', 
    'LA':'syc', 'LC':'ykc', 'LD':'swc', 'LE':'swg', 'LF':'ytc', 'LG':'skc', 'LH':'csc', 'LI':'mtc', 'LK':'mwa', 'LM':'mtg', 'LN':'mwc', 'LP':'cyc', 'LQ':'cyg', 'LR':'ckg', 'LS':'tyr', 'LT':'myc', 'LV':'stg', 'LW':'tkg', 'LY':'ywc', 
    'MA':'ryg', 'MC':'wks', 'MD':'rws', 'ME':'rwg', 'MF':'wts', 'MG':'rrg', 'MH':'mws', 'MI':'ats', 'MK':'awg', 'ML':'mtg', 'MN':'aws', 'MP':'myg', 'MQ':'mwg', 'MR':'akg', 'MS':'aks', 'MT':'ayg', 'MV':'rtg', 'MW':'wkg', 'MY':'wws', 
    'NA':'rmc', 'NC':'wrc', 'ND':'rac', 'NE':'ras', 'NF':'wwc', 'NG':'rrc', 'NH':'mac', 'NI':'awc', 'NK':'aas', 'NL':'mwc', 'NM':'aws', 'NP':'mmc', 'NQ':'mas', 'NR':'ars', 'NS':'arc', 'NT':'amc', 'NV':'rwc', 'NW':'wrs', 'NY':'wac', 
    'PA':'scc', 'PC':'ysc', 'PD':'smc', 'PE':'smg', 'PF':'yyc', 'PG':'ssc', 'PH':'cmc', 'PI':'myc', 'PK':'mma', 'PL':'cyc', 'PM':'myg', 'PN':'mmc', 'PQ':'cmr', 'PR':'csc', 'PS':'ycc', 'PT':'mcc', 'PV':'syc', 'PW':'ysg', 'PY':'ymc', 
    'QA':'sma', 'QC':'yrs', 'QD':'sas', 'QE':'sag', 'QF':'yws', 'QG':'grg', 'QH':'cas', 'QI':'mya', 'QK':'maa', 'QL':'cyg', 'QM':'mwg', 'QN':'mas', 'QP':'cmr', 'QR':'crg', 'QS':'yma', 'QT':'mma', 'QV':'swg', 'QW':'yrg', 'QY':'yas', 
    'RA':'rsa', 'RC':'ygc', 'RD':'src', 'RE':'rrg', 'RF':'ykc', 'RG':'sgg', 'RH':'crc', 'RI':'aka', 'RK':'arg', 'RL':'ckg', 'RM':'akg', 'RN':'ars', 'RP':'csc', 'RQ':'crg', 'RS':'ags', 'RT':'asc', 'RV':'rkc', 'RW':'ygg', 'RY':'yrc', 
    'SA':'kcc', 'SC':'tsc', 'SD':'rrc', 'SE':'kmg', 'SF':'tyc', 'SG':'rgc', 'SH':'mrc', 'SI':'akc', 'SK':'ars', 'SL':'tyr', 'SM':'aks', 'SN':'arc', 'SP':'ycc', 'SQ':'yma', 'SR':'ags', 'ST':'asc', 'SV':'kyc', 'SW':'tsg', 'SY':'tmc', 
    'TA':'rcc', 'TC':'wsc', 'TD':'rmc', 'TE':'rmg', 'TF':'wyc', 'TG':'rsc', 'TH':'mmc', 'TI':'ayc', 'TK':'ama', 'TL':'myc', 'TM':'ayg', 'TN':'amc', 'TP':'mcc', 'TQ':'mma', 'TR':'asc', 'TS':'asc', 'TV':'ryc', 'TW':'wsg', 'TY':'wmc', 
    'VA':'gyc', 'VC':'kkc', 'VD':'gwc', 'VE':'gwg', 'VF':'ktc', 'VG':'gkc', 'VH':'swc', 'VI':'rtc', 'VK':'rwa', 'VL':'stg', 'VM':'rtg', 'VN':'rwc', 'VP':'syc', 'VQ':'swg', 'VR':'rkc', 'VS':'kyc', 'VT':'ryc', 'VW':'kkg', 'VY':'kwc', 
    'WA':'ksg', 'WC':'tgs', 'WD':'krs', 'WE':'rrg', 'WF':'tks', 'WG':'kgg', 'WH':'yrs', 'WI':'wks', 'WK':'wrg', 'WL':'tkg', 'WM':'wkg', 'WN':'wrs', 'WP':'ysg', 'WQ':'yrg', 'WR':'ygg', 'WS':'tsg', 'WT':'wsg', 'WV':'kkg', 'WY':'trs', 
    'YA':'kmc', 'YC':'trc', 'YD':'kac', 'YE':'kas', 'YF':'twc', 'YG':'krc', 'YH':'yac', 'YI':'wwc', 'YK':'was', 'YL':'ywc', 'YM':'wws', 'YN':'wac', 'YP':'ymc', 'YQ':'yas', 'YR':'yrc', 'YS':'tmc', 'YT':'wmc', 'YV':'kwc', 'YW':'trs',
}
# This dictionary provides the best degenerate codon for every combination of two amino acids for EColi
AA_Pair_lookup_EColi = { 
    'AC':'ksc', 'AD':'gmc', 'AE':'gma', 'AF':'kyc', 'AG':'gsc', 'AH':'smc', 'AI':'ryc', 'AK':'rma', 'AL':'syc', 'AM':'ryg', 'AN':'rmc', 'AP':'scc', 'AQ':'sma', 'AR':'rsa', 'AS':'kcc', 'AT':'rcc', 'AV':'gyc', 'AW':'ksg', 'AY':'kmc', 
    'CA':'ksc', 'CD':'krc', 'CE':'krs', 'CF':'tkc', 'CG':'kgc', 'CH':'yrc', 'CI':'wkc', 'CK':'wrs', 'CL':'ykc', 'CM':'wks', 'CN':'wrc', 'CP':'ysc', 'CQ':'yrs', 'CR':'ygc', 'CS':'tsc', 'CT':'wsc', 'CV':'kkc', 'CW':'tgs', 'CY':'trc', 
    'DA':'gmc', 'DC':'krc', 'DE':'gas', 'DF':'kwc', 'DG':'grc', 'DH':'sac', 'DI':'rwc', 'DK':'ras', 'DL':'swc', 'DM':'rws', 'DN':'rac', 'DP':'smc', 'DQ':'sas', 'DR':'src', 'DS':'rrc', 'DT':'rmc', 'DV':'gwc', 'DW':'krs', 'DY':'kac', 
    'EA':'gma', 'EC':'krs', 'ED':'gas', 'EF':'kws', 'EG':'grg', 'EH':'sas', 'EI':'rwa', 'EK':'rag', 'EL':'swg', 'EM':'rwg', 'EN':'ras', 'EP':'smg', 'EQ':'sag', 'ER':'rrg', 'ES':'kmg', 'ET':'rmg', 'EV':'gwg', 'EW':'rrg', 'EY':'kas', 
    'FA':'kyc', 'FC':'tkc', 'FD':'kwc', 'FE':'kws', 'FG':'kkc', 'FH':'ywc', 'FI':'wtc', 'FK':'wwm', 'FL':'ytc', 'FM':'wts', 'FN':'wwc', 'FP':'yyc', 'FQ':'yws', 'FR':'ykc', 'FS':'tyc', 'FT':'wyc', 'FV':'ktc', 'FW':'tks', 'FY':'twc', 
    'GA':'gsc', 'GC':'kgc', 'GD':'grc', 'GE':'grg', 'GF':'kkc', 'GH':'src', 'GI':'rkc', 'GK':'rra', 'GL':'skc', 'GM':'rrg', 'GN':'rrc', 'GP':'ssc', 'GQ':'grg', 'GR':'sgg', 'GS':'rgc', 'GT':'rsc', 'GV':'gkc', 'GW':'kgg', 'GY':'krc', 
    'HA':'smc', 'HC':'yrc', 'HD':'sac', 'HE':'sas', 'HF':'ywc', 'HG':'src', 'HI':'mwc', 'HK':'mas', 'HL':'cwc', 'HM':'mws', 'HN':'mac', 'HP':'cmc', 'HQ':'cas', 'HR':'crc', 'HS':'mrc', 'HT':'mmc', 'HV':'swc', 'HW':'yrs', 'HY':'yac', 
    'IA':'ryc', 'IC':'wkc', 'ID':'rwc', 'IE':'rwa', 'IF':'wtc', 'IG':'rkc', 'IH':'mwc', 'IK':'awa', 'IL':'mtc', 'IM':'ats', 'IN':'awc', 'IP':'myc', 'IQ':'mya', 'IR':'aka', 'IS':'akc', 'IT':'ayc', 'IV':'rtc', 'IW':'wks', 'IY':'wwc', 
    'KA':'rma', 'KC':'wrs', 'KD':'ras', 'KE':'rag', 'KF':'wwm', 'KG':'rra', 'KH':'mas', 'KI':'awa', 'KL':'mwa', 'KM':'awg', 'KN':'aas', 'KP':'mma', 'KQ':'maa', 'KR':'arg', 'KS':'ars', 'KT':'ama', 'KV':'rwa', 'KW':'wrg', 'KY':'was', 
    'LA':'syc', 'LC':'ykc', 'LD':'swc', 'LE':'swg', 'LF':'ytc', 'LG':'skc', 'LH':'csc', 'LI':'mtc', 'LK':'mwa', 'LM':'mtg', 'LN':'mwc', 'LP':'cyc', 'LQ':'cyg', 'LR':'ckg', 'LS':'tyr', 'LT':'myc', 'LV':'stg', 'LW':'tkg', 'LY':'ywc', 
    'MA':'ryg', 'MC':'wks', 'MD':'rws', 'ME':'rwg', 'MF':'wts', 'MG':'rrg', 'MH':'mws', 'MI':'ats', 'MK':'awg', 'ML':'mtg', 'MN':'aws', 'MP':'myg', 'MQ':'mwg', 'MR':'akg', 'MS':'aks', 'MT':'ayg', 'MV':'rtg', 'MW':'wkg', 'MY':'wws', 
    'NA':'rmc', 'NC':'wrc', 'ND':'rac', 'NE':'ras', 'NF':'wwc', 'NG':'rrc', 'NH':'mac', 'NI':'awc', 'NK':'aas', 'NL':'mwc', 'NM':'aws', 'NP':'mmc', 'NQ':'mas', 'NR':'ars', 'NS':'arc', 'NT':'amc', 'NV':'rwc', 'NW':'wrs', 'NY':'wac', 
    'PA':'scc', 'PC':'ysc', 'PD':'smc', 'PE':'smg', 'PF':'yyc', 'PG':'ssc', 'PH':'cmc', 'PI':'myc', 'PK':'mma', 'PL':'cyc', 'PM':'myg', 'PN':'mmc', 'PQ':'cmr', 'PR':'csc', 'PS':'yct', 'PT':'mcc', 'PV':'syc', 'PW':'ysg', 'PY':'ymc', 
    'QA':'sma', 'QC':'yrs', 'QD':'sas', 'QE':'sag', 'QF':'yws', 'QG':'grg', 'QH':'cas', 'QI':'mya', 'QK':'maa', 'QL':'cyg', 'QM':'mwg', 'QN':'mas', 'QP':'cmr', 'QR':'crg', 'QS':'yma', 'QT':'mma', 'QV':'swg', 'QW':'yrg', 'QY':'yas', 
    'RA':'rsa', 'RC':'ygc', 'RD':'src', 'RE':'rrg', 'RF':'ykc', 'RG':'sgg', 'RH':'crc', 'RI':'aka', 'RK':'arg', 'RL':'ckg', 'RM':'akg', 'RN':'ars', 'RP':'csc', 'RQ':'crg', 'RS':'ags', 'RT':'asc', 'RV':'rkc', 'RW':'ygg', 'RY':'yrc', 
    'SA':'kcc', 'SC':'tsc', 'SD':'rrc', 'SE':'kmg', 'SF':'tyc', 'SG':'rgc', 'SH':'mrc', 'SI':'akc', 'SK':'ars', 'SL':'tyr', 'SM':'aks', 'SN':'arc', 'SP':'yct', 'SQ':'yma', 'SR':'ags', 'ST':'asc', 'SV':'kyc', 'SW':'tsg', 'SY':'tmc', 
    'TA':'rcc', 'TC':'wsc', 'TD':'rmc', 'TE':'rmg', 'TF':'wyc', 'TG':'rsc', 'TH':'mmc', 'TI':'ayc', 'TK':'ama', 'TL':'myc', 'TM':'ayg', 'TN':'amc', 'TP':'mcc', 'TQ':'mma', 'TR':'asc', 'TS':'asc', 'TV':'ryc', 'TW':'wsg', 'TY':'wmc', 
    'VA':'gyc', 'VC':'kkc', 'VD':'gwc', 'VE':'gwg', 'VF':'ktc', 'VG':'gkc', 'VH':'swc', 'VI':'rtc', 'VK':'rwa', 'VL':'stg', 'VM':'rtg', 'VN':'rwc', 'VP':'syc', 'VQ':'swg', 'VR':'rkc', 'VS':'kyc', 'VT':'ryc', 'VW':'kkg', 'VY':'kwc', 
    'WA':'ksg', 'WC':'tgs', 'WD':'krs', 'WE':'rrg', 'WF':'tks', 'WG':'kgg', 'WH':'yrs', 'WI':'wks', 'WK':'wrg', 'WL':'tkg', 'WM':'wkg', 'WN':'wrs', 'WP':'ysg', 'WQ':'yrg', 'WR':'ygg', 'WS':'tsg', 'WT':'wsg', 'WV':'kkg', 'WY':'trs', 
    'YA':'kmc', 'YC':'trc', 'YD':'kac', 'YE':'kas', 'YF':'twc', 'YG':'krc', 'YH':'yac', 'YI':'wwc', 'YK':'was', 'YL':'ywc', 'YM':'wws', 'YN':'wac', 'YP':'ymc', 'YQ':'yas', 'YR':'yrc', 'YS':'tmc', 'YT':'wmc', 'YV':'kwc', 'YW':'trs',
}
# A reverse lookup for degerate bases
Degenerate_Base_lookup = {
    'a': 'a', 'c': 'c','g': 'g','t': 't',
    'r': 'ag','y': 'ct','m': 'ac','k': 'gt','s': 'cg','w': 'at',
    'h': 'act','b': 'cgt','v': 'acg','d': 'agt',
    'n': 'acgt'
    }
# An amino-acid key for IQ-Tree *.state files.
AA_key = [
    'A','R','N','D','C','Q','E','G','H','I','L','K','M','F','P','S','T','W','Y','V'
    ]


# Reads in fasta file and renames certain sequences based on forbidden
# characters in IQ Tree as needed
def fasta2dict(fasta_path, return_dict={}):
    # Read in the file and prepare some variables
    with open(fasta_path, 'r') as infile:
        fastafile = infile.readlines()
    working_sequence = ''
    key = None
    for line in fastafile:  # For every line in the input fasta
        if line[0] == '>':  # Check if it's a name
            if working_sequence != '':  # If we already have a working sequence, another name indicates we're done. Otherwise record the name
                if len(working_sequence) >= 0:
                    return_dict[key] = working_sequence
            key = line[1:].rstrip().replace(':', '_')
            key = key.replace('|', '_')
            key = line[1:].rstrip()
            working_sequence = ''  # clear the working sequence
        else:  # If the line isn't a name, it's part of the sequence.
            if not all([char for char in working_sequence if (
                    char.isalpha() or char == '-')]):
                raise ValueError(
                    f"The provided file ({fasta_path}) was not a fasta file.")
            working_sequence = working_sequence + line.rstrip()
    # Finally, clean up anything left in the working sequence before returning
    # the dictionary.
    if working_sequence != '':
        if not all([char for char in working_sequence if (
                char.isalpha() or char == '-')]):
            raise ValueError(
                f"The provided file ({fasta_path}) was not a fasta file.")
        else:
            if len(working_sequence) >= 0:
                return_dict[key] = working_sequence
    if len(return_dict) == 0:
        raise ValueError(
            f"The provided file ({fasta_path}) was not a fasta file.")
    return return_dict


# Saves dictionary to fasta where the dictionary key is the protein name
# and the value is the sequence
def dict2fasta(fasta_d, fpath):
    with open(fpath, 'w+') as out_fasta:
        for key, seq in fasta_d.items():
            out_fasta.write(f'>{key}\n{seq}\n')


def Is_Valid_AA(AA):  # Is the argurment a valid amino acid or list of amino acids
    if isinstance(
            AA, str):  # This block lets us evaluate strings of amino acids
        return all([True for residue in AA if (
            (residue in AA_key) or (residue == '-'))])
    if isinstance(
            AA,
            list) and isinstance(
            AA[0],
            str):  # This block lets us evaluate lists of amino acids
        return all((i in AA_key) or (i == '-') for i in AA)
    if isinstance(
            AA,
            list) and isinstance(
            AA[0],
            list):  # This block lets us evaluate lists of lists of amino acids
        return all(all(((i in AA_key) or (i == '-'))
                   for i in lst) for lst in AA)
    else:
        raise ValueError("A bad type was evaluated as an amino acid list....")


def Is_Valid_Codon(codon):  # Is the argurment a valid codon
    dna_leters = ['a','c','t','g','r','y','m','k','s','w','h','b','v','d','n']
    if isinstance(codon, str):
        return (((len(codon) % 3) == 0) and all(
            i in dna_leters for i in codon))
    if isinstance(codon, list):
        return all(((len(c) == 3) and all(i in dna_leters for i in c))
                   for c in codon)


# Interact with BlastP, and record the XML
def NCBI_to_XML(dirname, sequence, hits=2000, expect_value=0.30, seq_name=''):
    # For the given sequence, we will run a BlastP search and parse the XML
    # return to make a multi-fasta file
    blast_result_handle = NCBIWWW.qblast(
        "blastp",
        "nr",
        sequence,
        hitlist_size=hits,
        expect=expect_value)
    XML_Name = f"BlastP_XML_{seq_name}"
    with open(f"./{dirname}/{XML_Name}", "w+") as fout:
        fout.write(blast_result_handle.read())
    blastp_xml = ET.parse(f"./{dirname}/{XML_Name}")
    return (blastp_xml)  # Returns the ElementTree type from the xml library


# Parse the BlastP XML - record Fasta
def Parse_BlastP_XML(dirname, blastp_xml, sequence, sequence_name=None):
    Fasta_Dict = {}
    return_name = "BlastP_Results.fasta"
    # If no sequence_name is identified (Meaning this is the search with the
    # user sequence)
    if sequence_name is None:
        if isinstance(sequence, str):
            with open(f"./{dirname}/HitsInfo.csv", "a+") as fout:
                fout.write(f"Hit ID,Hit Description,Hit Sequence\n")
                # Parsing the XML object, looking for hits
                for hit in blastp_xml.findall(
                        './BlastOutput_iterations/Iteration/Iteration_hits/Hit'):
                    # I've tried to also add the Hit_accession, but I can't
                    # access that form the XML for some reason
                    hitid = (hit.find('Hit_id')).text
                    hitdef = (hit.find('Hit_def')).text
                    hitaccession = (
                        hit.find('Hit_accession')).text.replace(
                        ".", "_")
                    seq = (hit.find('Hit_hsps/Hsp/Hsp_hseq')).text
                    # If the sequence doesn't have unknowns amino acids (annoying) then record it.
                    # The optional second method also removes exceptionally short or long sequences - be sure to synch with the code ~13 lines below
                    # if (("X" not in seq) and
                    # (len(seq)<((1+length_cutoff)*User_Sequence_Length)) and
                    # (len(seq)>((1-length_cutoff)*User_Sequence_Length))):
                    if ("X" not in seq):
                        fout.write(f"{hitid},{hitdef},{seq}\n")
                        Fasta_Dict[hitaccession] = seq
            with open(f"{dirname}/{return_name}", "a+") as blastp_file:
                for key, F_D_Sequence in Fasta_Dict.items():
                    # if
                    # (len(Sequence)<((1+length_cutoff)*User_Sequence_Length))
                    # and
                    # (len(Sequence)>((1-length_cutoff)*User_Sequence_Length)):
                    if len(F_D_Sequence) > 0.5 * len(sequence) and len(F_D_Sequence) < 1.5 * \
                            len(sequence) and Is_Valid_AA(F_D_Sequence):
                        # We remove all gaps, because CD-Hit cannot handle
                        # gaps.
                        blastp_file.write(
                            f">{key}\n{F_D_Sequence.replace('-','')}\n")
            # Modify the BlastP return to remove duplicate sequences.
            Remove_Duplicate_Sequences_FASTA(dirname, f"{return_name}")
            with open(f"{dirname}/{return_name}", "a+") as blastp_file:
                blastp_file.write(f">User_Sequence\n{sequence}\n")
        # If we've been given multiple sequences, we have a dict for sequence
        # and a list of xmls for blastp_xml
        elif isinstance(sequence, dict):
            for xml in blastp_xml:
                with open(f"./{dirname}/HitsInfo.csv", "a+") as fout:
                    # .csv header
                    fout.write(f"Hit ID,Hit Description,Hit Sequence\n")
                    # Parsing the XML objects, looking for hits
                    for hit in xml.findall(
                            './BlastOutput_iterations/Iteration/Iteration_hits/Hit'):
                        # I've tried to also add the Hit_accession, but I can't
                        # access that form the XML for some reason
                        hitid = (hit.find('Hit_id')).text
                        hitdef = (hit.find('Hit_def')).text
                        hitaccession = (
                            hit.find('Hit_accession')).text.replace(
                            ".", "_")
                        seq = (hit.find('Hit_hsps/Hsp/Hsp_hseq')).text
                        # If the sequence doesn't have unknowns amino acids (annoying) then record it.
                        # The optional second method also removes exceptionally short or long sequences - be sure to synch with the code ~13 lines below
                        # if (("X" not in seq) and
                        # (len(seq)<((1+length_cutoff)*User_Sequence_Length))
                        # and
                        # (len(seq)>((1-length_cutoff)*User_Sequence_Length))):
                        if ("X" not in seq):
                            fout.write(f"{hitid},{hitdef},{seq}\n")
                            Fasta_Dict[hitaccession] = seq
                for key, seq in sequence.items(
                ):  # Ensure all user seqeucnes are in final fasta file
                    if key not in Fasta_Dict:
                        Fasta_Dict.update({key: seq})
                with open(f"{dirname}/{return_name}", "a+") as blastp_file:
                    for key, F_D_Sequence in Fasta_Dict.items():
                        # if
                        # (len(Sequence)<((1+length_cutoff)*User_Sequence_Length))
                        # and
                        # (len(Sequence)>((1-length_cutoff)*User_Sequence_Length)):
                        if len(F_D_Sequence) < 10000 and len(
                                F_D_Sequence) > 10 and Is_Valid_AA(F_D_Sequence):
                            # We remove all gaps, because CD-Hit cannot handle
                            # gaps.
                            blastp_file.write(
                                f"\n>{key}\n{F_D_Sequence.replace('-','')}")
            # Modify the BlastP return to remove duplicate sequences.
            Remove_Duplicate_Sequences_FASTA(dirname, return_name)
        return (return_name)
    # If a sequence_name has been provided, this means we're doing Supplement
    # searches so our output directory structure should be different.
    elif isinstance(sequence_name, str):
        with open(f"{dirname}/{sequence_name}_supplement.csv", "a+") as fout:  # DIFFERENT FROM ABOVE
            fout.write(f"Hit ID,Hit Description,Hit Sequence\n")
            # Parsing the XML object, looking for hits
            for hit in blastp_xml.findall(
                    './BlastOutput_iterations/Iteration/Iteration_hits/Hit'):
                # I've tried to also add the Hit_accession, but I can't access
                # that form the XML for some reason
                hitid = (hit.find('Hit_id')).text
                hitdef = (hit.find('Hit_def')).text
                hitaccession = (
                    hit.find('Hit_accession')).text.replace(
                    ".", "_")
                seq = (hit.find('Hit_hsps/Hsp/Hsp_hseq')).text
                # If the sequence doesn't have unknowns amino acids (annoying) then record it.
                # The optional second method also removes exceptionally short or long sequences - be sure to synch with the code ~13 lines below
                # if (("X" not in seq) and
                # (len(seq)<((1+length_cutoff)*User_Sequence_Length)) and
                # (len(seq)>((1-length_cutoff)*User_Sequence_Length))):
                if ("X" not in seq):
                    fout.write(f"{hitid},{hitdef},{seq}\n")
                    Fasta_Dict[hitaccession] = seq
        # DIFFERENT FROM ABOVE
        with open(f"{dirname}/{sequence_name}_{return_name}", "a+") as blastp_file:
            # DIFFERENT FROM ABOVE
            blastp_file.write(f">{sequence_name}\n{sequence}\n")
            for key, F_D_Sequence in Fasta_Dict.items():
                if Is_Valid_AA(F_D_Sequence):
                    # if
                    # (len(Sequence)<((1+length_cutoff)*User_Sequence_Length))
                    # and
                    # (len(Sequence)>((1-length_cutoff)*User_Sequence_Length)):
                    # We remove all gaps, because CD-Hit cannot handle gaps.
                    blastp_file.write(
                        f">{key}\n{F_D_Sequence.replace('-','')}\n")
        # Modify the BlastP return to remove duplicate sequences. #DIFFERENT
        # FROM ABOVE
        Remove_Duplicate_Sequences_FASTA(
            dirname, f"{sequence_name}_{return_name}", True)
        return (f"{sequence_name}_{return_name}")
    else:
        print("sequence_name was not a string type")
        raise ValueError("sequence_name was not a string type")


# This function takes an amino acid sequence, submitts a BlastP search,
# and records the result in a fasta file
def BlastP(dirname, sequence, hits=2000, expect_value=0.2, sequence_name=None):
    if not (os.path.isdir(dirname)):
        os.mkdir(dirname)
    if isinstance(sequence, str):
        sequence = sequence.replace('-', '')
        if bool([char for char in sequence if (char not in AA_key)]
                ):  # Empty list (all chars in AA_key) comes up true
            print(
                f"Invalid sequence submitted for BlastP search. Contains the invalid charecter(s) {[char for char in sequence if (char not in AA_key)]}")
            raise ValueError("Invalid sequence submitted for BlastP search")
    elif isinstance(sequence, dict):
        for key, item in sequence.items():
            item = item.replace('-', '')
            sequence.update({key: item})
            if bool([char for char in item if (char not in AA_key)]):
                print(
                    f"Invalid sequence submitted for BlastP search ({key}). Contains the invalid charecter(s) {[char for char in item if (char not in AA_key)]}")
                raise ValueError(
                    f"Invalid sequence submitted for BlastP search ({key})")
    try:
        print("Acessing the NCBI database....")
        if isinstance(sequence, str):
            # No need to wait for an API call if this is just a re-run and we
            # already have the xmls available
            if not os.path.exists(f"{dirname}/BlastP_XML"):
                blastp_xml = NCBI_to_XML(dirname, sequence, hits, expect_value)
            else:
                blastp_xml = ET.parse(f"./{dirname}/BlastP_XML")
        elif isinstance(sequence, dict):
            xmls = []  # List of ElementTree objects
            for key, item in sequence.items():
                if not os.path.exists(
                        f"{dirname}/BlastP_XML_{key.replace('.','_')}"):
                    xmls.append(
                        NCBI_to_XML(
                            dirname,
                            item,
                            hits,
                            expect_value,
                            key.replace(
                                '.',
                                '_')))
                else:
                    xmls.append(
                        ET.parse(f"./{dirname}/BlastP_XML_{key.replace('.','_')}"))
    except BaseException:
        print("There was an error fetching the BlastP results")
        raise RuntimeError("There was an error fetching the BlastP results.")
    try:
        # Now, we parse the XML object and make a multi-fasta file. We also write the fasta file which is the result of our BlastP search.
        # The output behavior must be different if this the user sequence or a
        # Supplement search though.
        if sequence_name is None and isinstance(
                sequence, str):  # A primary search with one user sequence
            if not os.path.exists(f"{dirname}/BlastP_Results.fasta"):
                return_string = Parse_BlastP_XML(dirname, blastp_xml, sequence)
                return (return_string)
            else:
                return ("BlastP_Results.fasta")
        # A primary search with multiple user sequences
        elif sequence_name is None and isinstance(sequence, dict):
            if not os.path.exists(f"{dirname}/BlastP_Results.fasta"):
                return_string = Parse_BlastP_XML(dirname, xmls, sequence)
                return (return_string)
            else:
                return ("BlastP_Results.fasta")
        elif isinstance(sequence_name, str):  # A supplement search
            if not os.path.exists(
                    f"{dirname}/{sequence_name}_BlastP_Results.fasta"):
                try:
                    return_string = Parse_BlastP_XML(
                        dirname, blastp_xml, sequence, sequence_name)
                    return (return_string)
                except BaseException:
                    return (False)
            else:
                return (f"{sequence_name}_BlastP_Results.fasta")
        else:
            print("Variable sequence_name was not a string type")
            raise ValueError("Variable sequence_name was not a string type")
    except BaseException:
        print("There was an error recording the BlastP Results")
        raise RuntimeError("There was an error recording the BlastP Results")


def Sequence_Processing(dirname, finname, sequence):
    Fasta_Dict = {}
    if isinstance(sequence, dict):
        with open(f"{dirname}/{finname}", 'a') as fin:
            fin.write("\n")
            for key, item in sequence.items():
                fin.write(f">{key}\n{item}\n")
    else:
        try:
            # Alignment is needed for the Hamming Distance
            os.system(
                f"{MAFFT_Executable} {dirname}/{finname} > {dirname}/Early_Alignment_temp.fasta")
        except BaseException:
            raise RuntimeError("There was an error running MAFFT.")
        Fasta_Dict = fasta2dict(f"{dirname}/Early_Alignment_temp.fasta", {})
        os.remove(f"{dirname}/Early_Alignment_temp.fasta")
        # We have now computed the hamming distance of all sequences.
        Hamming_Dict = Fasta_Dict_Hamming(
            Fasta_Dict, Fasta_Dict["User_Sequence"])
        for key, item in Hamming_Dict.items():
            # If a given sequence has less than 50% similarity with the user
            # sequence, remove it.
            if (item / len(Fasta_Dict["User_Sequence"])
                    ) > 0.5 and key != "User_Sequence":
                Fasta_Dict.pop(key)
    for key, item in Fasta_Dict.items():
        # We now need to remove all the gaps in all the sequences to use CD-Hit
        Fasta_Dict.update({key: item.replace("-", "")})
    if supplement_Cutoff != 0:  # If we are doing supplementation
        # This one line supplements poorly supported areas of the tree, and
        # updates the dictionary with those additional sequences.
        Fasta_Dict.update(
            fasta2dict(f"{dirname}/{Supplement_Sequences(dirname,Fasta_Dict)}"))
    # Save the supplemented Fasta_Dict
    dict2fasta(Fasta_Dict, f"{dirname}/pre-Alignment.fasta")
    try:
        # Align the supplemented sequences - called the Raw alignment
        os.system(
            f"{MAFFT_Executable} {dirname}/pre-Alignment.fasta > {dirname}/Raw_Alignment.fasta")
    except BaseException:
        raise RuntimeError("There was an error running MAFFT.")
    # clean up our clutter, a little
    os.remove(f"{dirname}/pre-Alignment.fasta")
    Fasta_Dict_Aligned = fasta2dict(f"{dirname}/Raw_Alignment.fasta")
    if isinstance(sequence, dict):
        # The Raw alignment can then be processed into an alignment that's
        # ready for IQTree
        user_seq_name = [name for name in sequence.keys(
        ) if name in Fasta_Dict_Aligned.keys()][0]
        return_name = Post_MAFFT_processing(
            dirname, Fasta_Dict_Aligned, user_seq_name)
    else:  # If only one
        return_name = Post_MAFFT_processing(dirname, Fasta_Dict_Aligned, False)
    return return_name


def CDHit_Cutoff(identity):
    if (identity <= 0.4 or identity > 1):
        raise ValueError("The CD-Hit identity is invalid")
    elif identity > 0.7:
        return (5)
    elif identity > 0.6:
        return (4)
    elif identity > 0.5:
        return (3)
    elif identity > 0.4:
        return (2)


# This function will run CD-Hit and also ensure that the user sequence(s)
# are not removed
def CDHit(dirname, fin, fout, cutoff, supplement=False):
    try:
        # Run CD-Hit with user parameters
        os.system(
            f'{CDHIT_Executable} -i {dirname}/{fin} -o {dirname}/{fout} -c {cutoff} -n {CDHit_Cutoff(cutoff)}')
    except BaseException:
        print("There was an error running CD-Hit.")
        raise RuntimeError("There was an error running CD-Hit.")
    fout_dict = fasta2dict(f"{dirname}/{fout}")
    realign = False
    if not supplement:
        if isinstance(User_Input_Sequence, dict):
            for key, seq in User_Input_Sequence.items(
            ):  # Ensure all user seqeucnes are in final fasta file
                if key not in fout_dict:
                    realign = True
                    fout_dict.update({key: seq.replace("-", '')})
        else:
            if 'User_Sequence' not in fout_dict.keys():
                realign = True
                fout_dict.update(
                    {'User_Sequence': User_Input_Sequence.replace("-", '')})
    if realign:
        dict2fasta(fout_dict, "Temp.fasta")
        try:
            # Alignment is needed for the Hamming Distance
            os.system(f"{MAFFT_Executable} Temp.fasta > {dirname}/{fout}")
            os.remove("Temp.fasta")
        except BaseException:
            raise RuntimeError("There was an error running MAFFT.")


# This function MODIFIES a fasta file to remove highly-duplicate sequences.
def Remove_Duplicate_Sequences_FASTA(dirname, fpath, supplement=False):
    CDHit(dirname, fpath, f"{fpath[:-6]}_temp.fasta", 0.99, supplement)
    os.remove(f"{dirname}/{fpath}")
    os.remove(f"{dirname}/{fpath[:-6]}_temp.fasta.clstr")
    # Replace the original file with one that has no duplicate sequences.
    os.system(f"mv {dirname}/{fpath[:-6]}_temp.fasta {dirname}/{fpath}")


# returns the Hamming distance for every sequence in an aligned fasta
# dictionary.
def Fasta_Dict_Hamming(Fasta_Dict, sequence):
    Hamming_dict = {}
    for key, value in Fasta_Dict.items():
        Hamming_Diff = 0
        if len(value) != len(sequence):
            print("Hamming distance should be computed with properly aligned sequences.")
            raise ValueError(
                "Hamming distance should be computed with properly aligned sequences.")
        for i, char in enumerate(value):
            if sequence[i] != char:
                Hamming_Diff += 1
        Hamming_dict[key] = Hamming_Diff
    return Hamming_dict


# Find areas of poor coverage on the current tree and get additional
# BlastP search results.
def Supplement_Sequences(dirname, fasta_dict):
    if not (os.path.isdir(f"{dirname}/Sequence_Supplement")
            ):  # Make a directory for this sequence Supplementing process
        os.mkdir(f"{dirname}/Sequence_Supplement")
    # Now we can begin to identify sequences which are similar enough to the
    # rest of the dataset to be retained, but dissimilar enough that they
    # would benefit from supplementing.
    # Be sure there are no gaps in this dictionary, as CD-Hit will reject
    # those.
    dict2fasta(
        fasta_dict,
        f"{dirname}/Sequence_Supplement/Sequences_to_be_supplemented.fasta")
    sequences_list_for_search = []
    CDHit(
        dirname,
        "Sequence_Supplement/Sequences_to_be_supplemented.fasta",
        "Sequence_Supplement/SingleClusters_Supplemented.fasta",
        supplement_Cutoff,
        True)  # Run CD-Hit with supplementation cutoff to identify 'loner' sequences
    # Looking at the CD-Hit clstr file (which holds cluster information)
    with open(f"{dirname}/Sequence_Supplement/SingleClusters_Supplemented.fasta.clstr", 'r') as fin:
        clusters = (fin.read()).split('>Cluster ')
    for cluster in clusters:  # For each cluster
        # Split off the header from each sequence name
        seqs = cluster.split('>')
        if len(seqs) == 2:  # If there is only one sequence in the cluster
            # Seperate and record the name.
            sequences_list_for_search.append((seqs[1])[:-6])
    # Now we can submit a BlastP search for all of the sequences that need to
    # be supplemented, and write all sequences together as one file.
    files_list = [
        f"{dirname}/Sequence_Supplement/Sequences_to_be_supplemented.fasta"]
    print(
        f"Supplementing {len(sequences_list_for_search)} sequences with poor neighboorhood coverage.")
    for sequence_name in sequences_list_for_search:
        if not os.path.exists(
                f"{dirname}/Sequence_Supplement/{sequence_name}_BlastP_Results.fasta"):
            # Now we submit a BlastP search, but override the default expect
            # and hits to get a more narrow set of sequences.
            file = BlastP(
                f"{dirname}/Sequence_Supplement",
                fasta_dict[sequence_name],
                100,
                0.01,
                sequence_name)
            if file:
                # We also record a list of all the output fasta files to
                # concatanate them together later.
                files_list.append(f"{dirname}/Sequence_Supplement/{file}")
    # Now we need to write all of our Supplemented sequence searches together
    # as one fasta file. Just smoosh 'em all together
    with open(f"{dirname}/Sequence_Supplement/Supplemented_BlastP_Sequences.fasta", "a+") as fout:
        for fname in files_list:
            with open(fname, 'r') as supfile:
                fout.write(supfile.read())
    Remove_Duplicate_Sequences_FASTA(
        dirname, "Sequence_Supplement/Supplemented_BlastP_Sequences.fasta", True)
    return (f"Sequence_Supplement/Supplemented_BlastP_Sequences.fasta")


# Remove sequnces that cause insertions in the alighment
def Remove_Insertions(
        fasta_dict,
        User_Sequence,
        deletion_percentage=0.02,
        termini_length=0.05):
    num_pos = len(User_Sequence)
    acceptable_num_gaps = round(
        len(User_Sequence.replace('-', '')) * deletion_percentage)
    if acceptable_num_gaps < 5:
        acceptable_num_gaps = 5
    user_gap_pos = [i for i in range(num_pos) if User_Sequence.startswith(
        '-', i)]  # record all positions with a gap in the user sequence
    # The object 'key_sequence_list' contains a list of tuples with the name
    # of the sequence, the sequence, and a list of all positions with an amino
    # acid, in that order.
    # Sorry for the mess of a line, but it prevents unessecary looping.
    key_sequence_gap_list = [
        (key, [
            i for i in range(num_pos) if not seq.startswith(
                '-', i)]) for key, seq in fasta_dict.items()]
    keys_to_pop = []
    # The ax is already at the root of the trees, and every tree that does not
    # produce good fruit will be cut down and thrown into the fire
    for key, no_gap_list in key_sequence_gap_list:  # Now, for all sequences
        count = 0
        for pos in no_gap_list:  # Evaluate the positions where it has amino acids
            if (pos in user_gap_pos) and (pos > num_pos * termini_length) and (pos < num_pos * (1 - termini_length)
                                                                               ):  # If it has an amino acid where the User_Sequence has a gap, and isn't in the N- or C- terminus
                count += 1  # Tally a strike
            else:
                # If we come to a place of agreement, reset the count.
                count = 0
            # If the number of insertions in a row (count) is higher than what
            # is acceptable
            if count > acceptable_num_gaps:
                if key != 'User_Sequence':
                    keys_to_pop.append(key)  # Record this key as one to remove
                # What further testimony do we need? We have heard it ourselves
                # from his own lips.
                break
        else:
            continue
    for key in keys_to_pop:
        fasta_dict.pop(key)


# Remove sequnces that cause insertions in the alighment
def Remove_Deletions(
        fasta_dict,
        User_Sequence,
        deletion_percentage=0.02,
        termini_length=0.05):
    num_pos = len(User_Sequence)
    acceptable_num_gaps = round(
        len(User_Sequence.replace('-', '')) * deletion_percentage)
    if acceptable_num_gaps < 2:
        acceptable_num_gaps = 2
    user_gap_pos = [
        i for i in range(num_pos) if User_Sequence.startswith(
            '-', i)]  # record all positions with a gap
    # The object 'key_sequence_list' contains a list of tuples with the name
    # of the sequence, the sequence, and a list of all positions with a gap,
    # in that order.
    # Sorry for the mess of a line, but it prevents unessecary looping.
    key_sequence_gap_list = [
        (key, [
            i for i in range(num_pos) if seq.startswith(
                '-', i)]) for key, seq in fasta_dict.items()]
    keys_to_pop = []
    # The ax is already at the root of the trees, and every tree that does not
    # produce good fruit will be cut down and thrown into the fire
    for key, gap_list in key_sequence_gap_list:  # Now, for all sequences
        count = 0
        for pos in gap_list:  # Evaluate the positions where it has gaps
            if (pos not in user_gap_pos) and (pos > num_pos * termini_length) and (pos < num_pos * (1 - termini_length)
                                                                                   ):  # If it has gaps where the User_Sequence does not, and isn't in the N- or C- terminus
                count += 1  # Tally a strike
            else:
                # If we come to a place of agreement, reset the count.
                count = 0
            # If the number of gaps in a row (count) is higher than what is
            # acceptable
            if count > acceptable_num_gaps:
                if key != 'User_Sequence':
                    keys_to_pop.append(key)  # Record this key as one to remove
                # What further testimony do we need? We have heard it ourselves
                # from his own lips.
                break
        else:
            continue
    for key in keys_to_pop:
        fasta_dict.pop(key)


def Clean_all_gaps(fasta_dict):
    # For each position, go thorugh all sequences. If any sequence has a
    # residue, leave the position in at the end.
    sequences_list = [i for i in fasta_dict.values()]
    pos_to_leave = []
    for pos in range(len(sequences_list[0])):
        for seq in sequences_list:
            if (seq[pos] != '-'):
                pos_to_leave.append(pos)
                break
        else:
            continue
    # Remove all positions of all gaps from all sequences in the alignment
    for key, sequence in fasta_dict.items():  # Remove all the gaps
        fasta_dict.update(
            {key: ''.join([char for i, char in enumerate(sequence) if i in pos_to_leave])})


# Modifications after the alignment, mostly having to do with gaps.
def Post_MAFFT_processing(
        dirname,
        fasta_dict,
        multisequence,
        dynamic_sequence_reduction=True):
    # These functions **MODIFY** the fasta_dict by doing what their names say
    # they do.
    old_user_length = 0
    if multisequence:  # Just get one of them
        User_Sequence = fasta_dict.get(multisequence)
        # Each of these functions mutate the alignment, so it's important to
        # repeat until we've finsihed finding all misalignments.
        # HOWEVER - with a multi-sequence input, this tends to be overzealous.
        # We'll only go through the loop twice.
        for i in range(multisequence_iterations):
            Remove_Insertions(fasta_dict, User_Sequence)  # Insertions
            Clean_all_gaps(fasta_dict)  # Clean
            User_Sequence = fasta_dict.get(multisequence)  # Update
            Remove_Deletions(fasta_dict, User_Sequence)  # Deletions
            for key, seq in fasta_dict.items(
            ):  # Now we realign the sequences rather than just removing all the gaps - more relaiable
                fasta_dict.update({key: seq.replace("-", '')})
            dict2fasta(fasta_dict, f"{dirname}/Post_Mafft_Working1.fasta")
            try:
                # Alignment is needed for the Hamming Distance
                os.system(
                    f"{MAFFT_Executable} {dirname}/Post_Mafft_Working1.fasta > {dirname}/Post_Mafft_Working2.fasta")
                fasta_dict = {}
            except BaseException:
                raise RuntimeError("There was an error running MAFFT.")
            fasta2dict(f"{dirname}/Post_Mafft_Working2.fasta", fasta_dict)
            os.system(
                f"rm {dirname}/Post_Mafft_Working1.fasta {dirname}/Post_Mafft_Working2.fasta")
            User_Sequence = fasta_dict.get(multisequence)  # Update
    else:
        User_Sequence = fasta_dict.get("User_Sequence")
        # Each of these functions mutate the alignment, so it's important to
        # repeat until we've finsihed finding all misalignments.
        while len(User_Sequence) != old_user_length:
            old_user_length = len(User_Sequence)  # Store length
            Remove_Insertions(fasta_dict, User_Sequence)  # Insertions
            Clean_all_gaps(fasta_dict)  # Clean
            User_Sequence = fasta_dict.get("User_Sequence")  # Update
            Remove_Deletions(fasta_dict, User_Sequence)  # Deletions
            for key, seq in fasta_dict.items(
            ):  # Now we realign the sequences rather than just removing all the gaps - more relaiable
                fasta_dict.update({key: seq.replace("-", '')})
            dict2fasta(fasta_dict, f"{dirname}/Post_Mafft_Working1.fasta")
            try:
                # Alignment is needed for the Hamming Distance
                os.system(
                    f"{MAFFT_Executable} {dirname}/Post_Mafft_Working1.fasta > {dirname}/Post_Mafft_Working2.fasta")
                fasta_dict = {}
            except BaseException:
                raise RuntimeError("There was an error running MAFFT.")
            fasta2dict(f"{dirname}/Post_Mafft_Working2.fasta", fasta_dict)
            os.system(
                f"rm {dirname}/Post_Mafft_Working1.fasta {dirname}/Post_Mafft_Working2.fasta")
            User_Sequence = fasta_dict.get("User_Sequence")  # Update
    dict2fasta(fasta_dict, f'{dirname}/Post_MAFFT_Cleaned_Penultimate.fasta')
    # If on, this code will reduce the sequence alignment down to less than
    # set number of sequences. It's on by default
    if dynamic_sequence_reduction and len(fasta_dict) > Final_Dataset_Size:
        keep_going = True
        identity = 0.98
        for key, item in fasta_dict.items():
            fasta_dict.update({key: item.replace("-", "")})
        dict2fasta(fasta_dict, f'{dirname}/Post_MAFFT_No_Gaps.fasta')
        # Otherwise, we'll keep lowering the CD-Hit cutoff until we get below
        # the final library size.
        while keep_going:
            try:
                os.system(f"rm {dirname}/Penultimate_Sequences.fasta")
            except BaseException:
                continue
            CDHit(
                dirname,
                "Post_MAFFT_No_Gaps.fasta",
                "Penultimate_Sequences.fasta",
                round(
                    identity,
                    3))  # Run CD-Hit
            identity -= 0.01
            with open(f"{dirname}/Penultimate_Sequences.fasta", "r") as fin:
                lines = fin.readlines()
            if len(lines) <= (Final_Dataset_Size * 2) or identity < 0.80:
                print(f"Number Reached. Lines={len(lines)}")
                keep_going = False
            if len(lines) <= 100:
                keep_going = False
                try:
                    os.system(f"rm {dirname}/Penultimate_Sequences.fasta")
                except BaseException:
                    continue
                CDHit(
                    dirname,
                    "Post_MAFFT_No_Gaps.fasta",
                    "Penultimate_Sequences.fasta",
                    round(
                        (identity + 0.01),
                        3))  # Run CD-Hit
        try:
            os.system(f"rm {dirname}/Post_MAFFT_No_Gaps.fasta")
            # align the sequences passed into the function
            os.system(
                f"{MAFFT_Executable} {dirname}/Penultimate_Sequences.fasta > {dirname}/Final_Sequences.fasta")
        except BaseException:
            print("There was an error creating the sequence alignemnt.")
            raise RuntimeError(
                "There was an error creating the sequence alignemnt.")
        return ("Final_Sequences.fasta")
    else:
        os.system(
            f'cp {dirname}/Post_MAFFT_Cleaned_Penultimate.fasta {dirname}/Final_Sequences.fasta')
        os.system(f'rm {dirname}/Post_MAFFT_Cleaned_Penultimate.fasta')
        return ("Final_Sequences.fasta")


# Construct a phylogenetic tree, and determine a good model.
def IQTree_Phylo(dirname, finname):
    if not os.path.isdir(f"{dirname}/IQTree_Phylo"):  # Make the directory
        os.mkdir(f"{dirname}/IQTree_Phylo")
    try:
        # IQTree section
        os.system(
            f'{IQTREE_Executable} -s {dirname}/{finname} -pre {dirname}/IQTree_Phylo/Phylo -alrt 20000 -bb 20000 -nt AUTO')
    except BaseException:
        print("There was an error building the phylogeny in IQTree.")
        raise RuntimeError(
            "There was an error building the phylogeny in IQTree.")


def IQTree_ASR(dirname, finname):  # Using the tree/model from the phylogney, do the ASR
    model = ''
    try:  # this can fail because for some reason HPCC tries to run it before the IQTree_Phylo has finished running.
        # we need to go to the .iqtree file from the phylogony and find the
        # best model for this data set. This will save a ton of time.
        with open(f"{dirname}/IQTree_Phylo/Phylo.iqtree", "r") as fin:
            iqtree_file_lines = fin.readlines()
        if "Best-fit model" in iqtree_file_lines[42]:
            model = iqtree_file_lines[42].split(':')[1].strip()
        # If we can't find the Best-fit model where we expect it to be, then
        # check the whole file.
        else:
            for line in iqtree_file_lines:
                if "Best-fit model" in line:
                    model = line.split(':')[1].strip()
    except BaseException:
        pass
    if not os.path.isdir(f"{dirname}/IQTree_ASR"):  # Make the directory
        os.mkdir(f"{dirname}/IQTree_ASR")
    # Preventing errors; if we can't find the model, then we'll just use the
    # default behavior.
    if model == '':
        model = 'MFP'
    try:
        # IQTree section
        # Ask what model is best for ASR, use that one here.
        os.system(f'{IQTREE_Executable} -s {dirname}/{finname} -te {dirname}/IQTree_Phylo/Phylo.treefile -pre {dirname}/IQTree_ASR/ASR -asr -m {model} -redo -nt AUTO')
    except BaseException:
        print("There was an error conducting ASR in IQTree.")
        raise RuntimeError("There was an error conducting ASR in IQTree.")
    with open(f"{dirname}/IQTree_Phylo/Phylo.treefile", "r") as fin:
        Phylo_Tree = fin.read()
    with open(f"{dirname}/IQTree_ASR/ASR.treefile", "r") as fin:
        ASR_Tree = fin.read()
    with open(f"{dirname}/Final_Tree.treefile", "w+") as fout:
        fout.write(f"{Phylo_Tree}{ASR_Tree}")
    # Returns a dictionary out of *.state file
    ASR_Statefile_Dict = Statefile_to_Dict(dirname, "IQTree_ASR/ASR.state")
    return (ASR_Statefile_Dict)


# Select ancestral nodes which have a high enough confidence to resurect
# ancestors from them.
def Select_Ancestor_Nodes(dirname):
    Supports = {}
    UFB_Supports = []
    SHALRT_Supports = []
    # This will evaluate the tree, which is in Newick format, and select nodes of high enough confidence to reconstruct an ancestral library.
    # The Phlyo.treefile has confidence values as (SH-aLRT/UFB), and the
    # ASR.treefile has node names.
    with open(f'{dirname}/IQTree_ASR/ASR.treefile', 'r') as ASRtreefin:
        ASRtreefile = ASRtreefin.read()
    with open(f'{dirname}/IQTree_Phylo/Phylo.treefile', 'r') as Phylotreefin:
        Phylotreefile = Phylotreefin.read()
    confident_nodes = []  # these will be the nodes with good values. This isn't quite the same value between *.contree and *.trefile, but close
    # Let's split off just the information for each node, which is stored
    # after every close parenthiesis.
    ASRnodes = [n.split(',')[0].split(':') for n in ASRtreefile.split(')')]
    Phylonodes = [n.split(',')[0].split(':') for n in Phylotreefile.split(')')]
    ASRnodes.pop(0)  # the first split is not actually a node.
    Phylonodes.pop(0)
    nodes_to_skip = []
    # We exclude the last node, because that's the root, which will mess up
    # the below lines.
    for i in range(len(ASRnodes) - 1):
        try:
            UFB_i = float(Phylonodes[i][0].split('/')[1])
            SHALRT_i = float(Phylonodes[i][0].split('/')[0])
            Supports[ASRnodes[i][0]] = ((UFB_i, SHALRT_i))
            UFB_Supports.append(UFB_i)
            SHALRT_Supports.append(SHALRT_i)
            if (SHALRT_i > 80) and (
                    UFB_i > 95):  # If the SH-aLRT >80% and the ultrafast bootstraping is >95%
                # Record the name of the high-confidence nodes.
                confident_nodes.append(ASRnodes[i][0])
        # This is nessecary because sometimes the info on a node is genuinely
        # too bad to write down.
        except BaseException:
            nodes_to_skip.append(ASRnodes[i][0])
            print(f"Failed to read support information on: {ASRnodes[i][0]}")
    SHALRT_mean = sum(SHALRT_Supports) / len(SHALRT_Supports)
    UFB_mean = sum(UFB_Supports) / len(UFB_Supports)
    with open(f"{dirname}/IQTree_Phylo/Supports.txt", 'w+') as fout:
        fout.write(
            f"The phylogenetic tree has a mean ultra-fast bootstrap support of {round(UFB_mean)}%. For this method, 95% is considered the threshold of quality.\n")
        fout.write(
            f"{len([n for n in UFB_Supports if n > 95])} out of {len(UFB_Supports)} nodes have an ultra-fast bootstrap support above 95%.\n")
        fout.write(
            f"The standard deviation of ultra-fast bootstrap support is {round(((1/len(UFB_Supports))*sum([(x-UFB_mean)**2 for x in UFB_Supports]))**0.50,1)}%.\n\n")
        fout.write(
            f"The phylogenetic tree has a mean SH-aLRT support of {round(SHALRT_mean)}%. For this method, 80% is considered the threshold of quality.\n")
        fout.write(
            f"{len([n for n in SHALRT_Supports if n > 80])} out of {len(SHALRT_Supports)} nodes have a SH-aLRT support above 80%.\n")
        fout.write(
            f"The standard deviation of SH-aLRT support is {round(((1/len(SHALRT_Supports))*sum([(x-SHALRT_mean)**2 for x in SHALRT_Supports]))**0.50,1)}%.\n")
    UFB = []
    SHALRT = []
    with open(f"{dirname}/IQTree_Phylo/Supports.csv", 'w+') as fout:
        fout.write("Node,Ultra-Fast Bootstrap,SH-aLRT\n")
        for i in range(len(ASRnodes) - 1):
            # see above- no KeyErrors allowed here
            if (ASRnodes[i][0]) not in nodes_to_skip:
                fout.write(
                    f"{ASRnodes[i][0]},{Supports[ASRnodes[i][0]][0]},{Supports[ASRnodes[i][0]][1]}\n")
                UFB.append(Supports[ASRnodes[i][0]][0])
                SHALRT.append(Supports[ASRnodes[i][0]][1])
    try:
        import matplotlib.pyplot as plt
        # Make the sequence prediction histogram
        plt.rcParams["figure.figsize"] = [7.00, 5.00]
        plt.rcParams["figure.autolayout"] = True
        n_bins = 101
        # Plot the histogram
        plt.hist(UFB, n_bins)
        plt.title('Ultra-Fast Bootstrap Node Confidence Values')
        plt.xlabel("Confidence (%)")
        # Save the histogram
        plt.savefig(f"{dirname}/IQTree_Phylo/UFB_Confidences.png")
        plt.clf()
    except BaseException:
        print("Failed to make histogram of Ultra-Fast Bootstrap node confidence values")
    try:
        # Make the sequence prediction histogram
        plt.rcParams["figure.figsize"] = [7.00, 5.00]
        plt.rcParams["figure.autolayout"] = True
        n_bins = 101
        # Plot the histogram
        plt.hist(SHALRT, n_bins)
        plt.title('SH-aLRT Node Confidence Values')
        plt.xlabel("Confidence (%)")
        # Save the histogram
        plt.savefig(f"{dirname}/IQTree_Phylo/SHaLRT_Confidences.png")
        plt.clf()
    except BaseException:
        print("Failed to make histogram of SH-aLRT node confidence values")
    # return the list of node names whose value is above the cutoff.
    return confident_nodes


# Parse the statefile into a dictionary and record the confidence values.
def Statefile_to_Dict(dirname, fname):
    # This is a dicitonary made out of the statefile - its keys are the node
    # names,
    statefile_dict = {}
    # and its values are a list of tuples with the amino acid and list of
    # amino acid distributions at each position.
    node_distribution = []
    statelines = []
    # This will parse the .state file for the desired nodes to get their AA
    # distribution
    # Read in each line, skipping the header.
    with open(f'{dirname}/{fname}', 'r') as statefin:
        for i, line in enumerate(statefin):
            if i > 8:
                statelines.append(line)
    # Now let's pull the data from each line
    working_node = statelines[0].split()[0]  # prime the working node
    for line in statelines:  # For every line in the state file
        line_list = line.split()  # Break up the line into columns
        # If we're still working on the same node
        if working_node == line_list[0]:
            int_line_list = list(map(float, line_list[3:]))
            # record the probability distribution and the maximum probability
            # of the position
            node_distribution.append(int_line_list)
        else:  # If we've come to the end of a node,
            # Add a key-value pair to the statefile dictionary that is the node
            # name and the node's list
            statefile_dict[working_node] = node_distribution
            working_node = line_list[0]  # update the working_node value
            node_distribution = []  # Clear the working node lists
            int_line_list = list(map(float, line_list[3:]))
            # add to the node list the probabbilty distribution and the maximum
            # probability at that position
            node_distribution.append(int_line_list)
    # Be sure to add the last node into the dictionaries too!
    statefile_dict[working_node] = node_distribution
    # Dictionary format looks like {NodeX:[[probs],[probs],[probs],...]}
    return (statefile_dict)


# Do a binary ASR to determine where gaps in the ancestral sequences
# reside. Returns the Binary_Statefile_Dict
def Binary_Gap_Analysis(dirname, finname):
    if not os.path.isdir(f"{dirname}/IQTree_Binary"):  # Make the directory
        os.mkdir(f"{dirname}/IQTree_Binary")
    binary_dict = fasta2dict(f"{dirname}/{finname}")
    # These lines of code are to fix the bizare bug that has been showoing up with the dynamic sequence reduction that I spent all of spring 2022 trying to fix.
    # The bug is that binary_dict somehow imported all of the non-alighed
    # sequences from post_mafft_cleaned.fasta along with the sequences from
    # Final_Sequences.fasta, but not in a redundant way?
    alignment_length = max([len(seq)] for seq in binary_dict.values())
    to_pop = []
    for key, seq in binary_dict.items():
        if len(seq) != alignment_length[0]:
            to_pop.append(key)
    for key in to_pop:
        binary_dict.pop(key)
    for key, seq in binary_dict.items(
    ):  # Make a binary alignment of the input fasta
        # This line is taken from Ben - how to give proper credit?
        binary_dict.update(
            {key: (''.join(['0' if aa == '-' else '1' for aa in seq]))})
        # Note that a 0 is a gap and a 1 is an AA
    # Write to file for IQTree
    dict2fasta(binary_dict, f"{dirname}/IQTree_Binary/Binary_Alignment.fasta")
    try:
        os.system(f'{IQTREE_Executable} -s {dirname}/IQTree_Binary/Binary_Alignment.fasta -te {dirname}/IQTree_Phylo/Phylo.treefile -pre {dirname}/IQTree_Binary/Binary -blfix -asr -m GTR2+FO -redo -nt AUTO')
    except BaseException:
        print("There was an error determining gaps in the ancestral sequence")
        raise RuntimeError(
            "There was an error determining gaps in the ancestral sequence")
    try:
        Binary_Statefile_Dict = Statefile_to_Dict(
            dirname, "IQTree_Binary/Binary.state")
    except FileNotFoundError:
        print("The binary gap analysis was not successful.")
    ASR_Statefile_Dict = Statefile_to_Dict(dirname, "IQTree_ASR/ASR.state")
    Consensus_Ancestors_with_Gaps = {}
    Pos_with_Gaps = {}  # dictionary of {NodeX:[list of gaps at NodeX]}
    # Find positions that are actually gaps in the ASR
    for node, item in Binary_Statefile_Dict.items():
        gap_pos = []  # list of positions with gaps
        for i, pos in enumerate(item):  # each position
            if float(pos[0]) > 0.5:  # If the posotion has majority gap
                # The reason I've written it this was is that when the chances
                # are  close together (0.501 to 0.499) IQTree puts a gap in the
                # binary gap analysis *facepalm*
                gap_pos.append(i)
        # At each node, record the positions in the ancestral sequence that has
        # majority node.
        Pos_with_Gaps[node] = gap_pos
    # Merge the Sequence ASR with the gap ASR
    for node, cons_list in ASR_Statefile_Dict.items(
    ):  # node is name, cons_list is list of (list of AA confidence values) for each position
        consensus_seq = ''
        for i, pos in enumerate(
                cons_list):  # pos is the position of ancestor at node
            # If this position at this node is likely a gap, add a gap to the
            # consensus sequence
            if i in (Pos_with_Gaps[node]):
                consensus_seq += '-'
            else:
                # Otherwise, add the amino acid from ASR
                consensus_seq += AA_key[pos.index(max(pos))]
        Consensus_Ancestors_with_Gaps[node] = consensus_seq
    # Clean_all_gaps(Consensus_Ancestors_with_Gaps)
    dict2fasta(Consensus_Ancestors_with_Gaps,
               f"{dirname}/Consensus_Ancestors_with_Gaps.fasta")
    Consensus_Ancestors_without_Gaps = {}
    for key, item in Consensus_Ancestors_with_Gaps.items():
        Consensus_Ancestors_without_Gaps.update({key: item.replace("-", "")})
    dict2fasta(Consensus_Ancestors_without_Gaps,
               f"{dirname}/Consensus_Ancestors_without_Gaps.fasta")
    return (Binary_Statefile_Dict)


# This function takes a list of AAs for ONE POSITION and makes a
# degenerate codon for them.
def Degenerate_Nucleotide_Codon(AA_List, source='EColi'):
    # This is an area where the library size could be significantly improved
    # upon - take a more detailed look here - high priority
    if len(AA_List) == 0:
        return ''
    while ('X' in AA_List):
        AA_List.remove('X')
        # Unknown amino acids can occur, but we're going to ignore them.
    # This function takes a list of AAs for ONE POSITION and makes a
    # degenerate codon for them.
    for AA in AA_List:
        if not Is_Valid_AA(AA):
            raise ValueError("That's not an amino acid abreviation....")
    # Codon-based Lookup
    # A list of all the codons that need to be coded for (eg,
    # ['atc','agc','gta'])
    codon_list = []
    if (source == 'EColi'):
        for AA in AA_List:
            codon_list.append(AA_to_Codon_Ecoli[AA])
    elif (source == 'Human'):
        for AA in AA_List:
            codon_list.append(AA_to_Codon_Human[AA])
    else:
        raise NameError(
            "Please Specify EColi or Human as an expression organism.")
    # The three-charecter degenerate codon made from the codons in codon_list
    degenerate_codon = ''
    for pos in range(3):  # For every position in the codon
        # bases_at_pos will be a list of all the nucleotide bases wanted at
        # this position of the codon.
        bases_at_pos = ''
        for codon in codon_list:  # For every codon
            # If that codon's base is not in the bases_at_pos already,
            if (codon[pos] not in bases_at_pos):
                bases_at_pos += codon[pos]  # Add it
        # What degenerate codon codes for the desired bases?
        degenerate_codon += Mixed_Bases_lookup[bases_at_pos]
    return degenerate_codon


def Build_DNA_Sequence(Primer_Request, source='EColi'):
    primer_txt = ''
    if source == 'EColi':
        for pos in Primer_Request:  # pos is the list of AA at the respective position
            if (len(pos) == 1):
                # For a gap in the sequence, we will record nothing.
                if (pos[0] != '-'):
                    primer_txt += (AA_to_Codon_Ecoli[pos[0]])
            elif (len(pos) == 2):  # For a request of two AAs we use a lookup table
                primer_txt += (AA_Pair_lookup_EColi[pos[0] + pos[1]])
            else:  # Otherwise we bruteforce the degenerate DNA
                primer_txt += (Degenerate_Nucleotide_Codon(pos))
    elif source == 'Human':  # Same as above
        for pos in Primer_Request:
            if (len(pos) == 1):
                if (pos[0] != '-'):
                    primer_txt += (AA_to_Codon_Human[pos[0]])
            elif (len(pos) == 2):
                primer_txt += (AA_Pair_lookup_Human[pos[0] + pos[1]])
            else:
                primer_txt += (Degenerate_Nucleotide_Codon(pos))
    else:
        raise NameError(
            "The requested organism is not specified - please specify Human or EColi")
    if not Is_Valid_Codon(primer_txt):
        raise ValueError(
            "The primer generated was not a valid DNA sequence. No clue how that happened. If you're seeing this error, it's proabbly caused by a bug.")
    return (primer_txt)


# Make Uncertianty Libraries for all ancestral sequences with a high
# enough confidence at thier node
def Make_Uncertianty_Libraries(
        dirname,
        ASR_Statefile_Dict,
        Binary_Statefile_Dict,
        Cutoff=0):
    if ((Cutoff > 0.5) or (Cutoff < 0)):  # Error catching
        raise ValueError(
            "DNA Library cutoff value is invalid. It must be a number between 0.5 and 1")
    if not (os.path.isdir(f"{dirname}/DNA_Libraries")):  # Make directory
        os.mkdir(f"{dirname}/DNA_Libraries")
    Good_Ancestor_Nodes = Select_Ancestor_Nodes(dirname)
    if Cutoff == 0:
        Degen_Seqs_Size = {}
        for cutoff_value in [0.25, 0.2, 0.15, 0.125, 0.1, 0.075, 0.05, 0.025]:
            # For every node of sufficently high quality, we're going to make a
            # DNA template with uncertianty cutoffs.
            for node in Good_Ancestor_Nodes:
                # A request is a list of (list of amino acids needed) at each
                # position
                node_request = []
                Positions = ASR_Statefile_Dict[node]
                for i, pos in enumerate(Positions):  # For every position,
                    # If this position isn't a gap as determined by the Binary
                    # ASR
                    if (Binary_Statefile_Dict[node][i][0]) < 0.5:
                        pos_AAs = []
                        for j, prob in enumerate(
                                pos):  # For every probability at that position,
                            if (prob) > cutoff_value:  # If the amino acid is above the threshold
                                pos_AAs.append(
                                    AA_key[j])  # Record that AA at that position
                        if not bool(pos_AAs):  # Empty lists evaluate as false
                            # If none of the amino acids have a high enough
                            # probability to pass the threshold, we'll just
                            # record the most likely one.
                            pos_AAs = AA_key[pos.index(max(pos))]
                        node_request.append(pos_AAs)
                # Make a DNA sequence with degnerate bases
                Degenerate_DNA = Build_DNA_Sequence(node_request)
                if node in Degen_Seqs_Size:  # If there's already an object
                    Degen_Seqs_Size[node].append(
                        Library_Size_Count(Degenerate_DNA))
                else:  # If we need to intialize the list for the node
                    Degen_Seqs_Size[node] = [
                        Library_Size_Count(Degenerate_DNA)]
                with open(f"{dirname}/DNA_Libraries/Library_{cutoff_value*100}%_Cutoff.fasta", 'a+') as fout:
                    fout.write(f">{node}\n{Degenerate_DNA}\n")
        with open(f"{dirname}/DNA_Libraries/Library_Size_Information.csv", 'w+') as fout:
            fout.write(
                f"Confidece Threshold,{','.join(Good_Ancestor_Nodes)},\n")
            for i, value in enumerate(
                    [0.25, 0.2, 0.15, 0.125, 0.1, 0.075, 0.05, 0.025]):
                row_list = [str(Degen_Seqs_Size[node][i])
                            for node in Good_Ancestor_Nodes]
                fout.write(
                    f"{value}% Confidence threshold library,{','.join(row_list)},\n")
            # The below is an older version that wrote column-wise instead of
            # row-wise. I changed so that I can more easily append to the
            # bottom when new libraries are generated.
            '''
            fout.write("Ancestral Node, Sequences in 25% Confidence Library,Sequences in 20% Confidence Library,Sequences in 15% Confidence Library,Sequences in 12.5% Confidence Library,Sequences in 10% Confidence Library,Sequences in 7.5% Confidence Library,Sequences in 5% Confidence Library,Sequences in 2.5% Confidence Library,\n")
            for key,degen_seq_size_list in Degen_Seqs_Size:
                fout.write(f"{key},{','.join(degen_seq_size_list)},\n")
                #fout.write(f"{key},{degen_seq_list[0]},{degen_seq_list[1]},{degen_seq_list[2]},{degen_seq_list[3]},{degen_seq_list[4]},\n")
            '''
    else:
        row_list = []
        # For every node of sufficently high quality, we're going to make a DNA
        # template with uncertianty cutoffs.
        for node in Good_Ancestor_Nodes:
            # A request is a list of (list of amino acids needed) at each
            # position
            node_request = []
            Positions = ASR_Statefile_Dict[node]
            for i, pos in enumerate(Positions):  # For every position,
                # If this position isn't a gap as determined by the Binary ASR
                if (Binary_Statefile_Dict[node][i][0]) < 0.5:
                    pos_AAs = []
                    for j, prob in enumerate(
                            pos):  # For every probability at that position,
                        if (prob) > Cutoff:  # If the amino acid is above the threshold
                            pos_AAs.append(
                                AA_key[j])  # Record that AA at that position
                    if not bool(pos_AAs):  # Empty lists evaluate as false
                        # If none of the amino acids have a high enough
                        # probability to pass the threshold, we'll just record
                        # the most likely one.
                        pos_AAs = AA_key[pos.index(max(pos))]
                    node_request.append(pos_AAs)
            # Make a DNA sequence with degnerate bases
            Degenerate_DNA = Build_DNA_Sequence(node_request)
            row_list.append(str(Library_Size_Count(Degenerate_DNA)))
            with open(f"{dirname}/DNA_Libraries/Library_{Cutoff*100}%_Cutoff.fasta", 'a+') as fout:
                fout.write(f">{node}\n{Degenerate_DNA}\n")
        # This will append a new line to the bottom of the .csv containing the
        # sizes for the library with the new treshold.
        with open(f"{dirname}/DNA_Libraries/Library_Size_Information.csv", 'a+') as fout:
            fout.write(
                f"{Cutoff*100}% Confidence threshold library,{','.join(row_list)},\n")


def Write_Confidences(dirname, ASR_Statefile_Dict, Binary_Statefile_Dict):
    nodes_data = {}
    for node, cons_list in ASR_Statefile_Dict.items(
    ):  # node is name, cons_list is list of (list of AA confidence values) for each position
        node_confidences = []
        for i, pos in enumerate(
                cons_list):  # pos is AA distribution the position i of ancestor at node node
            # If this position at this node has a greater than 50% chance of
            # being a gap
            if (Binary_Statefile_Dict[node][i][0]) > 0.5:
                # append the confidnce that there is a gap.
                node_confidences.append(Binary_Statefile_Dict[node][i][0])
            else:  # if there's an amino acid
                # append the confidnce of that amino acid.
                node_confidences.append(max(pos))
        # confidences is now a list of the confidence of the most likely state
        # for each position at the node.
        # Each node has a tuple which is the average confidence and the number
        # of positions where confidence is below 85%
        nodes_data[node] = ((sum(node_confidences) / len(node_confidences)),
                            (len([i for i in node_confidences if i < 0.85])))
    with open(f"{dirname}/IQTree_ASR/Ancestral_Sequence_Confidences.csv", "w+") as fout:
        fout.write(
            "Node,Average Confidence,Positions below 85 percent confidence\n")
        for node, data in nodes_data.items():
            fout.write(f"{node},{round((data[0]*100),3)},{data[1]}\n")
    with open(f"{dirname}/IQTree_ASR/Ancestral_Sequence_Confidences.txt", "w+") as fout:
        fout.write(f"For the ASR in {dirname} overall:\n\
            {len(ASR_Statefile_Dict)} nodes have an average confidence of {round((sum([n[0] for n in nodes_data.values()])/len(nodes_data)*100),2)}% \n\
            These sequences have an average of {round(sum([n[1] for n in nodes_data.values()]) / len(nodes_data),1)} positions below 85% confidence out of {len(node_confidences)} total positions.\n\
            The topology of these nodes can be found at {dirname}/IQTree_Phylo/Phylo.contree.\n")
    try:
        import matplotlib.pyplot as plt
        # Make the sequence prediction histogram
        plt.rcParams["figure.figsize"] = [7.00, 5.00]
        plt.rcParams["figure.autolayout"] = True
        n_bins = 101
        # Plot the histogram
        data = [round((n[0] * 100), 3) for n in nodes_data.values()]
        plt.hist(data, n_bins)
        plt.title('Average Confidnece of each Ancestral Sequence')
        plt.xlabel("Confidence (%)")
        # Save the histogram
        plt.savefig(f"{dirname}/IQTree_ASR/Confidences.png")
        plt.clf()
    except BaseException:
        print(
            "Failed to make histogram of the average confidnece of each ancestral sequence.")


def seq_heatmap(Seqs, fname, n_col=50):  # Written by Patrick Finneran
    # Seqs is a dictionary where labels are the keys
    # and the values are two dictionaries: seq and score
    # seq is used for the annotation
    # score is used for the color value
    # Seqs = {node_name:{'seqs':sequence,'score':scores}}
    # Sequences must be aligned
    # Sequence must be a string
    # Score must be a list of same length of Sequence

    Sequences = {}
    Scores = {}
    for label, item in Seqs.items():
        seqName = label
        seq_length = len(item['seq'])
        Sequences[label] = []
        Scores[label] = []
        for i in range(seq_length):
            Sequences[label].append(item['seq'][i])
            Scores[label].append(item['score'][i])
    remainder = seq_length % n_col
    for i in range(n_col - remainder):
        for label in Seqs.keys():
            Sequences[label].append(' ')
            Scores[label].append(1)
    # Setup DataFrames
    df_Scores = pd.DataFrame(Scores)
    df_Scores = df_Scores.T
    df_Seqs = pd.DataFrame(Sequences)
    df_Seqs = df_Seqs.T

#    fig, ax = plt.subplots(nrows=int(ceil(seq_length/n_col)+1),figsize=(10,int(ceil(seq_length/n_col)+1)))
    fig, ax = plt.subplots(nrows=int(ceil(seq_length / n_col) + 1),
                           figsize=(10, int(ceil(seq_length / n_col) + 1) / 3))
    for i in range(int(ceil(seq_length / n_col))):
        h1 = sns.heatmap(df_Scores[df_Scores.columns[i * n_col:(i + 1) * n_col]],
                         square=True,
                         xticklabels=False,
                         vmin=0,
                         vmax=1,
                         annot=df_Seqs[df_Seqs.columns[i * n_col:(i + 1) * n_col]],
                         fmt='',
                         ax=ax[i],
                         cbar=False,
                         cmap='hot')
        ax[i].set_yticklabels(ax[i].get_yticklabels(), rotation=0)
    cmap_ticks = [0]
    for i in range(5):
        cmap_ticks.append(cmap_ticks[-1] + 1 / 5)
    mappable = h1.get_children()[0]
    plt.colorbar(mappable,
                 cax=ax[-1],
                 ticks=cmap_ticks,
                 orientation='horizontal')
    title = '{0} Consensus at Each Position'.format(seqName)
    ax[0].set_title(title, fontsize=14)
    plt.savefig(fname, format='pdf', dpi=1200, bbox_inches='tight')
    plt.close()


def Library_Size_Count(Degenerate_Sequence):
    if not Is_Valid_Codon(Degenerate_Sequence):
        raise ValueError("Not a valid DNA sequence")
    Codons = [Degenerate_Sequence[i:i + 3]
              for i in range(0, len(Degenerate_Sequence) - 2, 3)]
    Library_Size = 1
    for Codon in Codons:
        CodedAminoAcids = []
        for base1 in Degenerate_Base_lookup[Codon[0]]:
            for base2 in Degenerate_Base_lookup[Codon[1]]:
                for base3 in Degenerate_Base_lookup[Codon[2]]:
                    CodedAminoAcids.append(Codon_to_AA[base1 + base2 + base3])
        if ' stop ' in CodedAminoAcids:
            Library_Size *= (len(set(CodedAminoAcids)) - 1)
        else:
            Library_Size *= len(set(CodedAminoAcids))
    return (Library_Size)


if __name__ == '__main__':
    parser = argparse.ArgumentParser(
        prog='AP-LASR',
        description='A script for the automation of the generation of protein combinatorial libraries by ancestral sequence reconstruction.',
        epilog='If you really need help, reach out to the author: jjvanantwerp [at] gmail [dot] com')
    # Adding subparsers allows this script to be run in multiple 'modes' for
    # different functions.
    subparsers = parser.add_subparsers(
        title='Modes of operation',
        dest='mode',
        help='Specify a mode to use AP-LASR: {python AP-LASR.py *mode*}.')
    parser_mode1 = subparsers.add_parser(
        'Assist', help='Get assistance in running AP-LASR.')  # Help mode
    # Normal use is ASR mode, which does ASR and builds libraries
    parser_mode2 = subparsers.add_parser(
        'ASR', help='Conduct ASR and generate protein libraries.')
    parser_mode2.add_argument(
        '-n',
        '--name',
        action='store',
        dest='directory',
        default='ASR',
        help='The name of the directory where ASR output will be stored (default \'ASR\').')
    parser_mode2.add_argument(
        '-i',
        '--input',
        action='store',
        dest='input',
        required=True,
        help='REQUIRED: sequence input as raw sequence text or FASTA file name')
    parser_mode2.add_argument(
        '-s',
        '--size',
        action='store',
        dest='FinalDatasetSize',
        default=500,
        type=int,
        help='The desired number of modern sequences in the final dataset. Defaults to 500.')
    parser_mode2.add_argument(
        '--supplement',
        action='store',
        dest='supplementationCutoff',
        type=float,
        help='supplementation is needed when dataset diversity must be improved. We reccomend beginning with a value of 75%.')
    parser_mode2.add_argument(
        '--cdhit',
        action='store',
        dest='cdhitexe',
        default='cd-hit',
        help='The executable for CDHit (default cd-hit)')
    parser_mode2.add_argument(
        '--mafft',
        action='store',
        dest='mafftexe',
        default='mafft',
        help='The executable for MAFFT (default mafft)')
    parser_mode2.add_argument(
        '--iqtree',
        action='store',
        dest='iqtreeexe',
        default='iqtree',
        help='The executable for IQTree (default iqtree)')
    parser_mode2.add_argument(
        '--MSI',
        action='store',
        dest='multisequence_iterations',
        default=2,
        type=int,
        help='The number of iterations in \"Post_MAFFT_Processing\" when using the multi-sequence input feature. Increase for smaller datasets, decrease for larger ones.')
    # This mode makes figures from previously generated ASR data
    parser_mode3 = subparsers.add_parser(
        'MakeFigures', help='Make figures from a previous ASR run\'s data.')
    parser_mode3.add_argument(
        '-n',
        '--name',
        action='store',
        dest='directory',
        default='ASR')
    # This mode generates new protein libraries from a previous ASR run with a
    # new cutoff value
    parser_mode4 = subparsers.add_parser(
        'RemakeLibraries',
        help='Generate new combinatorial libraries for ancestors of a previous ASR run with a specified cutoff value.')
    parser_mode4.add_argument(
        '-n',
        '--name',
        action='store',
        dest='directory',
        default='ASR')
    parser_mode4.add_argument(
        '-t',
        '--threshold',
        action='store',
        dest='threshold',
        required=True,
        type=float)
    # parser.set_defaults(func=lambda args: parser.print_help())
    args = parser.parse_args()
    if (args.mode == 'ASR') or (args.mode == 'MakeFigures') or (
            args.mode == 'RemakeLibraries'):
        directory = args.directory

    if args.mode == 'Assist':
        print("\n")
        print(Help_String)
        print("\n")
        print(Software_Prerequistes)
        print("\n\n")
        print(Directory_Structure)
        print("\n\n")

    if args.mode == 'ASR':
        CDHIT_Executable = args.cdhitexe
        MAFFT_Executable = args.mafftexe
        IQTREE_Executable = args.iqtreeexe
        Final_Dataset_Size = args.FinalDatasetSize
        multisequence_iterations = args.multisequence_iterations
        if not (Final_Dataset_Size > 50):
            print("User provided value for the final dataset size was too small. Preceding with a final dataset size of 60 sequences.")
            Final_Dataset_Size = 60
        if args.supplementationCutoff is not None:
            if (100 > args.supplementationCutoff > 40):
                print(
                    f"Interpreting supplementation cutoff {args.supplementationCutoff} as {args.supplementationCutoff}%.")
                supplement_Cutoff = (args.supplementationCutoff / 100)
            else:
                supplement_Cutoff = args.supplementationCutoff
            if not (0.9 > supplement_Cutoff > 0.4) or (supplement_Cutoff == 0):
                print("User provided value for the supplementation cutoff was outside the valid range. Please specify a value between 90% and 40%.")
                raise ValueError(
                    "User provided value for the supplementation cutoff was outside the valid range. Please specify a value between 90% and 40%.")
        else:
            supplement_Cutoff = 0
        if '.' not in args.input:
            sequence = args.input.upper().replace("-", '')
            if any([True for n in sequence if n not in AA_key]):
                raise ValueError(
                    "Sequence option read as raw sequence - The provided sequence could not be used. Please be sure no amino acids are \'X\'")
            Blastp_out_name = BlastP(directory, sequence)
            Final_Name = Sequence_Processing(
                directory, Blastp_out_name, sequence)
        elif (os.path.exists(f"{directory}/Final_Sequences.fasta")):
            print(
                f"Previously completed BalstP search and dataset curation detected. Proceding with {directory}/Final_Sequences.fasta to IQTree reconstruction.")
            Final_Name = "Final_Sequences.fasta"
        else:
            if not (os.path.exists(args.input)):
                raise ValueError(
                    "The specified input fasta file does not exist.")
            try:
                User_Input_Sequence = {}
                temp_User_Input_Sequence = fasta2dict(args.input)
                for key, item in temp_User_Input_Sequence.items():
                    if '.' in key:
                        key = (key.split("."))[0]
                    User_Input_Sequence.update({key: item})
            except BaseException:
                raise ValueError("The file could not be read as a fasta file.")
            if len(User_Input_Sequence) == 1:
                User_Input_Sequence = list(User_Input_Sequence.values())[0]
            else:
                print(
                    "Detected input fasta contained multiple sequences. This is an acceptable input, but is more prone to errors.")
            Blastp_out_name = BlastP(directory, User_Input_Sequence)
            Final_Name = Sequence_Processing(
                directory, Blastp_out_name, User_Input_Sequence)
        with open(f"{directory}/{Final_Name}") as fin:
            if len(fin.readlines()) < 80:
                print("WARNING: There are very few sequences in the final sequence alignment. This can indicate many highly similar sequences were returned from BlastP or that there are not many sequences available.")

        IQTree_Phylo(directory, Final_Name)
        ASR_Statefile_Dict = IQTree_ASR(directory, Final_Name)
        Binary_Statefile_Dict = Binary_Gap_Analysis(directory, Final_Name)
        Write_Confidences(directory, ASR_Statefile_Dict, Binary_Statefile_Dict)
        Make_Uncertianty_Libraries(
            directory,
            ASR_Statefile_Dict,
            Binary_Statefile_Dict)

    if args.mode == 'MakeFigures':
        if not os.path.exists(f"{directory}/IQTree_Binary"):
            print(
                "Could not make heatmap figures - ASR has not been completed. Please finish ASR first.")
            sys.exit()
        import matplotlib.pyplot as plt
        import pandas as pd
        from math import ceil
        import seaborn as sns
        # Let's make an object of form {NodeX:([sequence],[confidence])}
        Consensus_Ancestors_with_Gaps_and_Confidence = {}
        ASR_Statefile_Dict = Statefile_to_Dict(
            directory, "IQTree_ASR/ASR.state")
        Binary_Statefile_Dict = Statefile_to_Dict(
            directory, "IQTree_Binary/Binary.state")
        # Find positions that are actually gaps in the ASR
        Pos_with_Gaps = {}  # dictionary of {NodeX:[list of gaps at NodeX]}
        for node, item in Binary_Statefile_Dict.items():
            gap_pos = []  # list of positions with gaps
            for i, pos in enumerate(item):  # each position
                if float(pos[0]) > 0.5:  # If the posotion has majority gap
                    # The reason I've written it this was is that when the
                    # chances are  close together (0.501 to 0.499) IQTree puts
                    # a gap in the binary gap analysis *facepalm*
                    gap_pos.append(i)
            # At each node, record the positions in the ancestral sequence that
            # has majority node.
            Pos_with_Gaps[node] = gap_pos
        # Merge the Sequence ASR with the gap ASR
        for node, cons_list in ASR_Statefile_Dict.items(
        ):  # node is name, cons_list is list of (list of AA confidence values) for each position
            consensus_seq = []
            confidence = []
            for i, pos in enumerate(
                    cons_list):  # pos is the position of ancestor at node
                # If this position at this node is likely a gap, add a gap to
                # the consensus sequence
                if i in (Pos_with_Gaps[node]):
                    consensus_seq.append('-')
                    confidence.append(Binary_Statefile_Dict[node][0])
                else:
                    # Otherwise, add the amino acid from ASR
                    consensus_seq.append(AA_key[pos.index(max(pos))])
                    confidence.append(pos)
            Consensus_Ancestors_with_Gaps_and_Confidence[node] = {
                "seqs": consensus_seq, 'score': confidence}
        for node, item in Consensus_Ancestors_with_Gaps_and_Confidence.items():
            seq_heatmap({node: item},
                        f"{directory}/Confidence_Heatmaps/{node}_Confidences.pdf")

    if args.mode == 'RemakeLibraries':
        if (100 > args.threshold > 50):
            print(
                f"Interpreting confidence threshold {args.threshold} as {args.threshold}%.")
            cutoff = (args.threshold / 100)
        else:
            cutoff = args.threshold
        if (cutoff > 0.5) or (cutoff < 0):
            print("Confidence threshold value must be between 0.5 and 0, with a reccomended value between 0.025 and 0.25.")
            raise ValueError(
                "Confidence threshold value must be between 0.5 and 0, with a reccomended value between 0.025 and 0.25.")
        # DO LIBARY MAKING
        try:
            ASR_Statefile_Dict = Statefile_to_Dict(
                directory, "IQTree_ASR/ASR.state")
            Binary_Statefile_Dict = Statefile_to_Dict(
                directory, "IQTree_Binary/Binary.state")
        except BaseException:
            print(
                "Failed to read in data from previous ASR run - be sure the correct directory is specified.")
            raise RuntimeError(
                "Failed to read in data from previous ASR run - be sure the correct directory is specified.")
        try:
            Make_Uncertianty_Libraries(
                directory,
                ASR_Statefile_Dict,
                Binary_Statefile_Dict,
                cutoff)
        except BaseException:
            print("Reconstruction of ASR libraries failed.")
            raise RuntimeError("Reconstruction of ASR libraries failed.")
