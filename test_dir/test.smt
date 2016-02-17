#  Substructural definitions for profiling sets of structures using
#  the BOTHROPS program
#  
#  These definitions are aimed at functional groups which are likely to
#  interact with target protein or nucleic acid in a directional manner
#  and it is intended that STRUCT_ANAL be used with the dt_vmatch (matches
#  only the head atoms of the SMARTS target  
#
#  
#  WARNING:  This file is subject to alteration at ANY time!  Make your own
#            copy NOW to avoid any unpleasant surprises LATER!
#
#___________________________________________________________________________
#
#  bothrops_groups_1:  PWK  11-2-97
#
#  Atom types/ functional groups defined:
#
#      1)  Lipohilic (i.e. not protonated) cations
#      2)  Probable anions at physiological pH
#      3)  Probable protonated cations at physiological pH
#      4)  Neutral hydrogen bond acceptors
#      5)  Neutral hydrogen bond donors
#      6)  Simultaneous donor/acceptor
#      7)  Potential electrophiles
#
#  Notes on these definitions:
#      
#      1)  The struct_anal program uses the dt_vmatch function for matching
#          SMARTS definitions against SMILES which contributes a count of 
#          exactly one for each successful match.
#      2)  Since pyridones contribute a single donor and a single acceptor
#          the problem of their being registered as the wrong tautomer has 
#          not been addressed.
#      3)  Nitro oxygen and halogens are not considered to be hydrogen bond
#          acceptors
#      4)  Amidines and guanidines are treated as protonated unless an 
#          electron withdrawing group is on one of the nitrogens
#      5)  Most acids and bases in the IPC are registered as neutral species
#          and definitions have been included to identify neutral groups
#          which are likely to be charged at physiological pH.  Atoms 
#          which are explicitly protonated or deprotonated are also provided
#          for. Specifically protonated aromatic nitrogen is treated as a
#          neutral hydrogen bond acceptor.
#      6)  Charged groups are counted directly rather than in terms of number
#          of protons removed/added.  Therefore phosphonate does not get
#          counted twice because two protons can be removed and guanidine
#          is not counted three times because charge is delocalised over
#          three nitrogens.  This does require some thought as to how the
#          groups are defined:
#          PHOSPHONATE = P(=O)([O;H,-][O;H,-]  not [O;H,-]P(=O)[O;H,-] which
#          will count the number of acidic protons
#
#--------------------------------------------------------------------------
#  MODIFICATIONS FROM ORIGINAL DEFINITIONS
#  
#_________________________________________________________________________
#  
#   Format:
# SMARTS_NAME  SMARTS_DEFINITION  SYMMETRY_NUMBER  ANALYSIS_FLAG
#
#   Notes on format:
#
#      1) SYMMETRY_NUMBER normally set to 1 but can to be set to integer N 
#         (>0) if substructural count will always lead to a multiple of N.
#         This allows the count to be scaled for symmetrical substructural
#         targets. Symmetry numbers other than 1 must not be used with 
#         these definitions
#      2) ANALYSIS_FLAG tells program whether to perform analysis for this 
#         target.  Can either be set to 1 (perform analysis) or 0 (do not 
#         perform analysis).   
#      3) "#" indicates comment line
#
#________________Some general definitions___________________________ 
# sp3 carbon
SP3    [C;X4]                                                   1   0
# Halogens
HALO   [F,Cl,Br,I]                                              1   0
#   Doubly connected aromatic nitrogen ( hydrogen bond acceptor)
NARM   [n;D2;!H;!+]                                             1   0
#   Triply connected aromatic nitrogen 
N1A    [n;D3;!+][$SP3]                                          1   0
N1B    [n;D2;H;!+]                                              1   0
N1     [$N1A,$N1B]                                              1   0
#   Amino and alkylamino substituents
N1SUB  [N;!+;H2]                                                1   0
N2SUB  [N;!+;H][$SP3]                                           1   0
N3SUB  [N;!+]([$SP3])[$SP3]                                     1   0
NSUB   [$N1SUB,$N2SUB,$N3SUB]                                   1   0
#  For defining delocalised anions
CHAN   [C;H,H2,-]                                               1   0
NHAN   [N;H,-]                                                  1   0
CNH    [$CHAN,$NHAN]                                            1   0
CO     C=O                                                      1   0
OSCO   [$CO,O,S]                                                1   0
#   Aromatic carbons: c-R (R=H, alkyl, amino or alkylamino)
#   These are defined for the identification of heterocycles which
#   are likely to be protonated at physiological pH  
C1A    [c;H]                                                    1   0
C1B    [c;!H][$SP3,$NSUB]                                       1   0
C1     [$C1A,$C1B]                                              1   0
#_____________Electron withdrawing groups____________________________
NITRO1 N(=O)=O                                                  1   0
NITRO2 [N;+]([O;-])=O                                           1   0
NITRO [$NITRO1,$NITRO2]                                         1   0
CYANO  C#N                                                      1   0
GENACY [C,c,S]([C,N,O&!H])=[O,S]                                1   0
IMIN   C=N                                                      1   0
ESINK  [$NITRO,$CYANO,$GENACY,$IMIN]                            1   0
#--------------------------------------------------------------------
#
#_____________ Anions (deprotonated at physiological pH)_____________
#  Acyl, sulfonyl and heterocyclic sulfonamides
ACSFM  [N;H,-]([C,S]=O)S(=O)=O                                  1   0
HTSFM1 [N;H,-](S(=O)(=O))cn                                     1   0
HTSFM2 [N;H,-](S(=O)(=O))c[a;R2][a;R2][$NARM]                   1   0
ARSFM1 [N;H,-](S(=O)(=O))cc[$ESINK]                             1   0
ARSFM2 [N;H,-](S(=O)(=O))c1ccc([$ESINK])cc1                     1   0
SFM    [$ACSFM,$HTSFM1,$HTSFM2,$ARSFM1,$ARSFM2]                 1   0
#  
#  Tetrazoles (consider both tautomers) and define the appropriate
#  nitrogen at the same time for elimination from neutral hydrogen
#  bond definitions
TETRZ1 c1[n;H,-]nnn1                                            1   0
TETRZ2 c1n[n;H,-]nn1                                            1   0
TETRZ  [$TETRZ1,$TETRZ2]                                        1   0
NTETZ1 n[$TETRZ]                                                1   0
NTETZ2 n[$NTETZ1]                                               1   0
#  Conventional oxyanions including explicit oxyanions
OACID  [C,S,P](=O)[O;H,-]                                       1   0
#  Hydroxamic acids
HYDXAM  [O;H,-]NC=O                                             1   0
#  Sulfonylnitromethanes
SUNM   [C;H,H2,-](S(=O)(=O))[$NITRO]                            1   0
#  Potential carbon acids                                        
CESNK3 [C;H,-]([$ESINK])([$ESINK])[$ESINK]                      1   0
#  Combine carbon acid definitions                              1   0
CACID  [$SUNM,$CESNK3]                                          1   0
#  Some delocalised anions found by DT 
#  See J. Med. Chem. v39 (1996) , p 5228
DELOC1 [O,S]=c1[n;H,-]cn[s,o]1                                  1   0
DELOC2 O=S1[N;H]C=NO1                                           1   0
#  See J. Med. Chem. v40 (1997), p 520
DELOC3 [O,S;H,-]c1n[o,s]cc1                                     1   0
#  The folowing were of interest in the insulin resistance 
#  project see also QSAR v 10 (1991) p109
DELOC4 O=C1[$CNH]C(=O)C[$OSCO]1                                  1   0
#  Combine delocalised anions
DELOC  [$DELOC1,$DELOC2,$DELOC3,$DELOC4]                         1   0
#    Combine anion definitions 
ANION  [$OACID,$SFM,$TETRZ,$HYDXAM,$CACID,$DELOC]               1   1
#--------------------------------------------------------------------
#
#_____________ Cations (protonated at physiological pH_______________
# Aliphatic amines including explicitly protonated species
AMIN1  [N;H2;!+][$SP3]                                          1   0
AMIN2  [N;H;!+]([$SP3])[$SP3]                                   1   0
AMIN3  [N;!+]([$SP3])([$SP3])[$SP3]                             1   0
AMIN4  [N;+;H,H2,H3]                                            1   0
AMIN   [$AMIN1,$AMIN2,$AMIN3,$AMIN4]                            1   0
# Amidines and guanidines
# Don't allow electron-withdrawing group on nitrogen and define nitrogen
# joined to amidine carbon so that it can be treated as a neutral acceptor
#
N2EW   N[$NITRO,$CYANO,$GENACY]                                 1   0
NAMDN  [N;!+;v3;!$N2EW;!$NITRO]                                 1   0
AMDIN  C(=[$NAMDN])([$NAMDN])[C,c,$NAMDN]                       1   0
N2AMDN N~[$AMDIN]                                               1   0  
# Aromatic nitrogens which may protonate
# Imidazole and its alkylated derivatives
IMDZ   [$NARM]1[$C1][$C1][$N1][$C1]1                            1   0
# Some six-rings which are expected to protonate 
#    2-Aminopyridine and its alkyl derivatives
PRD1   [$NARM]1c([$NSUB])[$C1][$C1][$C1][$C1]1                  1   0
#    4-Aminopyridine and its alkyl derivatives
PRD2   [$NARM]1[$C1][$C1]c([$NSUB])[$C1][$C1]1                  1   0
#    2,4-Diaminopyrimidine and alkyl derivatives  
PRD3   [$NARM]1c([$NSUB])[$NARM]c([$NSUB])[$C1][$C1]1           1   0
#    Combine cation definitions
PROT   [$AMIN,$IMDZ,$AMDIN,$PRD1,$PRD2,$PRD3]                   1   1
#
#____________ Neutral hydrogen bond acceptors________________________
#  Aliphatic alcohols
OL     [O;H][$SP3]                                              1   0
#   Aliphatic ethers (aromatic ethers are not regarded as acceptors
#                     in this set of definitions)
ETHER  [O;D2;!+]([$SP3])[$SP3]                                  1   0
#   Ketone, amide, ester and sulfur or phosphorus analogues
ONE1   [O,S]=[C,S,P]([C,c,N,n])[C,c,N,n,O&!H&!-,S&!H&!-]        1   0
#   N-oxides (anticipate different registration schemes)
ONE2   O=[n&D3,N&D4]                                            1   0
ONE3   [O;-][n&D3,N&D4;+]                                       1   0
#  Formyl derivatives (this definition will not catch formaldehyde but
#             who other than a mortician is interested in formaldehyde)
ONE4   O=[C;H][C,c,N,n,O&!H&!-,S&!H&!-]                         1   0
#   Aliphatic nitrogens (start with the junk)
NAL1   N[$CYANO,$NITRO,$GENACY]                                 1   0
NAL2   Ncn                                                      1   0
NAL3   Nccn                                                     1   0
NAL4   Ncccn                                                    1   0
NALJ   [$NAL1,$NAL2,$NAL3,$NAL4]                                1   0
NAL    [N;v3;!+;!$PROT;!$NALJ;!$N2AMDN]                         1   0
# Protonated aromatic nitrogen is likely to deprotonate when it gets
# buffered to physiological pH so treat at a neutral; but beware of
# legitimately protonated aromatic nitrogen
NARMH  [n;+;H;!$PRD2]                                           1   0 
#   Combine acceptor definitions 
#   (including aromatic nitrogen from above)
ACC1   [$ETHER,$ONE1,$ONE2,$ONE3,$ONE4,$OL]                     1   0
ACC2   [$NARM,$NARMH,$NAL;!$NTETZ1;!$NTETZ2;!$PRD2]             1   0
ACC    [$ACC1,$ACC2]                                            1   1
#
#--------------------------------------------------------------------
#
#____________ Neutral hydrogen bond donors___________________________
#  
PHOL   [O;H]c                                                   1   0
AMD1   [N;H,H2][S,P,C,c,N]=[O,S]                                1   0
AMD2   [n;H;!+;!$NTETZ1;!$NTETZ2]                               1   0
AMD3   [N;H,H2]c                                                1   0
AMD4   [N;H,H2][$CYANO,$NITRO1,$NITRO2]                         1   0
DON    [$OL,$PHOL,$AMD1,$AMD2,$AMD3,$AMD4]                      1   1
#
#___Simultaneous donor & acceptor____________________________________
#
PHENOL [O;H]c1ccccc1                                            1   0
DOAC   [$OL,$PHENOL]                                            1   1
#
#--------------------------------------------------------------------
#_____________ Cations (but not protonated)__________________________
#
NOTCAT [N,n,P,S,s;+][O,S;-]                                     1   0
NCAT1  [N,n;+;D3;!H]                                            1   0
NCAT2  [N,P;+;D4;!H]                                            1   0
OSCAT  [o,O,s,O;!H;+]                                           1   0
LIPCAT [$NCAT1,$NCAT2,$OSCAT;!$NOTCAT]                          1   1
#--------------------------------------------------------------------
#
#-------------Potential to complex with metals (especially Zinc)-----
#
THIOL   [S;H]                                                   1   0
THIEST  [S;D2;!+]([C,c])C=O                                     1   0
METAL   [$THIOL,$THIEST,$OACID,$NARM,$HYDXAM]                   1   1
#
#--------------------------------------------------------------------
#
#____________ Potential electrophiles________________________________ 
#
# Michael acceptors ( >1 electron withdrawing groups on olefin)_____
MICH1  C([$ESINK])=C[$ESINK]                                    1   0
MICH2  C=C([$ESINK])[$ESINK]                                    1   0
MICH   [$MICH1,$MICH2]                                          1   0
#  Reactive carbonyls and analogues
ANHYD  C(=O)OC=O                                                1   0
ACHAL  C(=O)[Cl,Br,I]                                           1   0
ALDHD  [C;H]=O                                                  1   0
KETF1  C(=O)C(F)(F)F                                            1   0
KETF2  C(O)(O)C(F)(F)F                                          1   0
KETF   [$KETF1,$KETF2]                                          1   0
#  Strained lactones and lactams                 
LACT1  C1(=O)[O,N,S]C1                                          1   0
LACT2  C1(=O)[O,N,S]CC1                                         1   0
LACT   [$LACT1,$LACT2]                                          1   0
#  Alkyl halides
ALHAL1 [C;H2]([C,c])[Cl,Br,I]                                   1   0
ALHAL2 [C;H]([C,c])([C,c])[Cl,Br,I]                             1   0
ALHAL3 C([C,c])([C,c])([C,c])[Cl,Br,I]                          1   0
ALHAL  [$ALHAL1,$ALHAL2,$ALHAL3]                                1   0
#Other halides
OTHHAL [Si,S,P][$HALO]                                          1   0
#  Epoxides and episulfides 
EPI    [O,S]1CC1                                                1   0
#
#  Boron compounds
BOR1   [B;D3;!H;!-]                                             1   0
BOR2   [B;D4;-][O,Cl]                                           1   0
BOR    [$BOR1,$BOR2]                                            1   0
#  Electrophilic aromatics
PYRIM1 [$HALO]c1ncnca1                                          1   0
PYRIM2 [$HALO]c1ncccn1                                          1   0
EAROM1 [$HALO]c1c([$ESINK])cc([$NITRO])cc1                      1   0
#  Combine definitions
EPHIL1 [$MICH,$ANHYD,$ACHAL,$ALHAL,$OTHHAL,$KETF,$ALDHD]        1   0
EPHIL2 [$EPI,$PYRIM1,$PYRIM2,$EAROM1,$BOR,$LACT]                1   0
EPHIL  [$EPHIL1,$EPHIL2]                                        1   1  
#_____________________________________________________________________

# rings
aring5 a1~a~a~a~a~1 1 1
aring6 a1~a~a~a~a~a~1 1 1
cacc *~[$ACC] 1 1
cdon *~[$DON] 1 1
cyclohex C1CCCCC1 1 1
pip N1CCNCC1 1 1
piperd N1CCCCC1 1 1
