

import unittest
from unittest import TestCase
import configparser
from lisa import LISA
import lisa.gene_selection as gene_selection
import numpy as np
import lisa.models as models
import lisa.assays as assays
import lisa.utils as utils

class TestGeneSelection(TestCase):

    @classmethod
    def setUpClass(cls):
        cls.lisa = LISA('hg38', verbose = 0)
        cls.geneset, _ = cls.lisa._load_gene_info()
        cls.query_list = 'PSG5,PSG3,GATA4,CD83,APOBEC2,DFNA5,S100A9,PSG4,SLC26A2,PDGFRA,HAS3,PRSS12,RGS8,S100A12,C8orf49,SLC9A7,TMEM65,PSG9,RGS5,SYTL5,CST1,PSG7,TNFAIP3,GPRC5B,CST4,CHN2,FAM135A,RYR2,RGN,RASGRP1,GCNT3,PSG6,APLF,SORBS2,GSTM4,TNIK,PPARGC1A,TM2D1,PRKD3,CYMP,ABCG2,ZNF114,GALC,SMIM24,SIPA1L2,GDPD1,MMD,PSG2,LINC01300,TMC5,FAM46C,THUMPD3-AS1,POLR1D,CDH11,TRABD2B,MBNL2,CACNA1D,GOLT1A,RPS29,FPR3,SELENOP,COL28A1,PSD3,HKDC1,UNC5CL,NLRP7,LINC00479,PTPRH,CCNYL1,SLC9A2,LINC01468,ELFN2,MECOM,PKP4,IGFBP5,KALRN,SSX2IP,LYZ,RBP2,ZNF704,RP11-248J18.2,TRPM7,PRR9,TMEM163,WSB1,NAV2,TMPRSS2,RP11-356I2.4,IGFL1,GDA,SGK1,C9orf152,CTH,RHOU,CALD1,HGSNAT,ZNF654,RRBP1,PRKCA,FAM134B'
        cls.query_list = [gene.upper().strip() for gene in cls.query_list.split(',')]

        cls.background_list = 'ORC2,CHST4,LSAMP,CCDC36,OR10P1,POLA2,CD1D,TCN1,SLC27A4,STK19,CLINT1,GNB1L,AK9,WDR41,NEUROG2,DAPK1,C19orf52,HRH4,CASS4,H2BFM,OR51E2,NOTCH3,LOXL4,HEATR9,ITGBL1,ITPRIP,GNAQ,NGEF,MUC21,FAF2,B4GALT1,CT45A5,COPG2,DBX2,BRCA2,PIGA,OR4D1,ABCA2,AHSG,BEX2,PRMT2,OTOP2,NDUFB1,BHLHA9,PHB2,NBN,PHYHD1,SCPEP1,THUMPD2,UBL7,ARMCX3,VHL,SYCP2L,TRIM55,C19orf66,TRAT1,CCL3L3,C1orf115,EFNA3,NOS1,C3orf67,RASSF6,MROH7,MIOS,ACBD3,ZDHHC22,LSM1,MS4A6E,HIST1H2AL,PPTC7,PTPN20,MOV10L1,ACADVL,DDX41,CAPNS1,CES1,WDR13,LIN28A,TEX29,UBL4A,ZNF324,IRF8,DDX5,SLC25A12,GNG2,PCGF3,PTCD2,CNTN4,C5AR1,AKR1B15,ATXN3,SERPINB11,AGAP1,TXNDC5,SLC35G5,GPRIN1,LRRC10B,DCANP1,AKR1D1,MYH13,KATNAL2,COL4A4,GREM2,SLC10A7,SEMA4C,MRPL16,SURF1,ZNF33A,ADAMTS12,PLA2G2E,SLC16A3,VSTM2L,POMK,SIN3A,ZNF839,DHX9,TBATA,GPR25,TXNDC12,DCBLD1,ATP6AP2,ZNF670,SMAD5,GOLGA6A,IFI27L2,ZNF625,HOXD11,PARD6G,HOMER3,CEP120,BARHL2,KRT75,SCYL1,OR1F1,MOB3C,ABCB8,C10orf2,PRAC1,RNF114,CADM3,ZNF510,CD37,DYNC1H1,SLITRK4,EXOSC4,NELFCD,UQCC3,PURA,SLC35F4,GSTA5,CALR,LYN,FAM126B,RBM5,MRFAP1,MGAT4A,UBXN2B,PRDM7,NOL8,RBFOX3,LRRTM4,COL6A3,SGTA,KIF2B,CHST2,CDKN2D,VMO1,GBP3,DSG2,FMNL3,LGALS13,KRTAP9-9,CFAP65,PCDHGA9,ACOT12,CFL2,FGF23,OR4F4,PCDHB16,C9orf135,ENTPD2,EPN3,PNMA5,SNCAIP,SYNE3,REM1,DACT3,FAM83D,COMMD8,MAP2,TSPAN9,RHPN1,IER5L,GOLGA8H,TMEM67,RALGDS,VTCN1,PSMB1,ALOX12,SMARCA1,PIGQ,MTRNR2L5,IRGM,LOC389199,CDC73,FGFBP2,DUSP13,RAPH1,UBE2D2,FOXP4,PCDHB6,SCAMP1,ADAMTS5,SPINT2,VPS53,SGK494,MSH6,RFX5,BRAT1,KIAA1033,LRRTM1,TMEM41B,TREH,STXBP5L,MYO15A,ZNF48,INS-IGF2,C10orf111,NMNAT2,DUSP22,RBP5,CARD16,OIT3,TNFRSF12A,CLUH,ZNF641,ZNF846,CLEC9A,ARMC7,SEMA4G,GLRX3,TKFC,TOR3A,PEX7,DIAPH1,MEI4,MAP4K2,GGACT,SLC41A3,ISY1-RAB43,OR5AS1,LRRC75B,AGBL4,OSR1,SOCS7,HTR3D,THEMIS,CCDC60,RAB19,OR5T2,FABP2,SLC22A24,CIB4,PSMD11,PHC1,TEX11,PUS3,CSNK2A3,PSME3,DDRGK1,EEF2K,OR4C46,KLHL6,MYH6,IDH3B,LRRC69,NBPF1,DSTN,CCDC174,CYB5A,CT47A2,IGFBP3,GLI1,ARSD,ZNF579,CAPN2,SLC2A9,MARCH7,GIP,CD99,TIRAP,RNF138,DUOXA2,SPANXD,PDE7A,EXOC3,LAMA4,C8orf89,ADCY5,AVPR1A,SSSCA1,LRRC4,SLC3A1,NDUFC2,ZNF773,METTL12,C2orf74,KLK3,LOC100131107,UBXN4,MRPL53,PAPOLB,CFHR2,FLI1,HIST1H1C,ZBTB4,CCPG1,CNPY4,GADL1,ABHD2,FAM219A,MET,PIEZO2,ARHGEF15,OR51G1,RAX2,SMDT1,QRICH1,ITPR2,CCL1,PBXIP1,NEU2,FNDC3A,PLD5,RFWD3,LY6E,RDH16,MFAP1,SLC39A14,F8,NANOS1,PFDN6,C17orf75,SLFN5,KPNB1,DNA2,OR2T8,FAM65B,PIK3CD,MARCH6,SIK2,SIGLEC5,TLR3,C11orf88,CSNK1G3,TRPV1,TMEM80,FAM183A,SLC16A14,CCT5,RSPRY1,NKAIN2,SLC27A2,PPHLN1,FABP5,CAMTA1,STAMBP,C1orf106,ZNF354A,KRTAP4-9,CNTD2,KIF22,TM7SF2,ARMCX4,ZNF705G,TUBA3C,ADGRL4,NCF2,SERF1A,PPP4R1,SMARCD3,SFTPC,SELE,MYC,KLF8,YEATS2,TRAK1,MAN1A1,F12,FAM107B,EVI2B,CYP2D6,CLEC4G,C18orf54,SYNRG,DMTN,FGF8,FHOD1,SYPL1,BAIAP2L1,SELP,AKR7A2,SLC25A22,RABGGTA,DCUN1D3,AFG3L2,SNAP25,CEP95,ACSL6,OR4D2,OBP2B,ARPC3,CHP2,PEX11G,MROH9,RNF219,NEIL1,RASL12,TRAPPC3,ALDH18A1,RAB25,KDSR,PRSS38,LTA,FUT1,AXIN1,COG2,TNIP2,SLC15A1,CHURC1-FNTB,LYSMD1,REV1,TMEM200C,MYL2,USP26,DLG2,BFSP1,PATL2,FUT5,PRR12,ACSS1,BTF3,PILRB,CSGALNACT2,FOXO3,ZNF506,CADM2,GPR63,PYCRL,RNF10,RPAP3,HSPB8,PAIP2B,VPS4A,CFHR3,CDC123,ANKUB1,CD1B,TVP23A,ZNF587B,NPTXR,APLP2,SPTLC2,FBXL3,ZNF701,GOLGA8J,MOSPD2,NFASC,CNGA1,SLC36A1,LIN54,CRABP2,HEMGN,UGT1A4,OR6B2,CPNE2,CXorf40B,C2orf61,ZADH2,FAM129A,ZWINT,APOH,SHISA9,OGG1,RASSF7,RDH14,ISLR,GPR87,GAD1,ASTN2,WDR54,MRAP2,MT1B,C22orf46,UQCRQ,PHOX2A,MEPE,TRIM72,PEX13'
        cls.background_list = [gene.upper().strip() for gene in cls.background_list.split(',')]

    def validate_sampling(self, label_dict, gene_info, num_genes = 3000):

        query_genes = [k for k, v in label_dict.items() if v]
        background_genes = [k for k,v in label_dict.items() if not v]

        self.assertEqual(len(query_genes), 97)
        self.assertEqual(len(background_genes), num_genes)
        self.assertEqual(len(set(query_genes).intersection(set(background_genes))), 0)
        self.assertEqual(len(set(gene_info['query_symbols']).intersection(set(gene_info['background_symbols']))), 0)

    def test_tad_background_sampling(self):

        label_dict, infodict = gene_selection.select_genes(self.query_list, self.geneset, seed = 0)

        self.validate_sampling(label_dict, infodict)

    def test_random_background_sampling(self):

        label_dict, infodict = gene_selection.select_genes(self.query_list, self.geneset, background_strategy= 'random')

        self.validate_sampling(label_dict, infodict)

    def test_background_provided(self):

        label_dict, infodict = gene_selection.select_genes(self.query_list, self.geneset, background_strategy='provided', 
            background_list= self.background_list)

        self.validate_sampling(label_dict, infodict, 499)     

    def test_gene_matching(self):

        genes_returned = self.geneset.match_user_provided_genes(self.query_list)

        self.assertEqual(len(genes_returned), 97)


def make_random_labels(num_samples, num_pos):
    label = np.zeros(num_samples)
    label[np.random.choice(num_samples, num_pos, replace = False)] = 1
    label = label.astype(np.bool)
    return label


class TestModels(TestCase):

    @classmethod
    def setUpClass(self):

        self.dataset_label = make_random_labels(1000, 10)
        self.gene_label = make_random_labels(3000, 100)
        
        matrix_shape = (self.gene_label.size, self.dataset_label.size)

        design_matrix = self.gene_label[:,np.newaxis] * self.dataset_label[np.newaxis, :]

        self.dataset_directionality = np.random.choice([-1.0,1.0], (1, matrix_shape[1])) # ~B(1/2)
        signal = self.dataset_directionality * (np.random.normal(0,1, matrix_shape) + np.random.lognormal(0, 0.2, matrix_shape)) # ~ N(lognormal(0,0.5), 1)
        background = np.random.normal(0, 1, matrix_shape) # ~N(0,1)

        rp_matrix =  background + design_matrix * signal
        self.rp_matrix = rp_matrix - rp_matrix.min()

        selector = models.LR_BinarySearch_SampleSelectionModel(200, 10)

        self.selections = selector.fit(self.rp_matrix, self.gene_label)
        
    def test_dataset_selection_model(self):

        dataset_label_index = np.where(self.dataset_label)[0]
    
        self.assertEqual(len(set(dataset_label_index).difference(set(self.selections))), 0)

    def test_chromatin_model(self):

        model_matrix = self.rp_matrix[:, self.selections]

        chrom_model = models.LR_ChromatinModel({'C' : list(10.0**np.arange(-2,4.1,0.5))}, penalty = 'l2')

        chrom_model.fit(model_matrix, self.gene_label)

        coefs = chrom_model.model.coef_

        self.assertTrue(((self.dataset_directionality.reshape(-1)[self.selections] > 0) == (coefs.reshape(-1) > 0)).all())

        test_datacube = np.random.rand(self.gene_label.size, self.dataset_label.size, 1000)

        self.assertTrue(chrom_model.get_deltaRP_activation(test_datacube).shape == (self.gene_label.size, 1000))


class TestPeakRP_Assay(TestModels):

    def setUp(self):

        self.model = assays.PeakRP_Assay(None, None, 1, utils.Log(verbose = False))

    def test_peak_rp_p_values(self):

        self.model.rp_matrix = self.rp_matrix

        pvals = self.model.predict(np.ones(len(self.gene_label)).astype(np.bool), self.gene_label)

        lowest_pvals = np.argsort(pvals)

        dataset_label = np.argwhere(np.logical_and(self.dataset_label.reshape(-1), self.dataset_directionality.reshape(-1) > 0))

        self.assertTrue(np.isin(dataset_label, lowest_pvals[:10]).all())


if __name__ == "__main___":
    unittest.main()




