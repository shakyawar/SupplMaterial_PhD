import csv
import string


##### Proteomics data ###################################################
#####################################################################
f_proteomics = csv.reader(open('LM_Allproteins_Pawar.csv', 'rb'), delimiter=',')
a=0
prot_geneID_list = []
for row in f_proteomics:
  a=a+1
  if a>1:
    prot_accNo = row[0]
    prot_geneID = row[1]
    prot_geneID = string.replace(prot_geneID,"LmjF.","LmjF")

    prot_name = "("+row[2]+")"
    prot_length = row[3]
    prot_MolWt = row[7]
    if a>2:
      prot_geneID_list.append(prot_geneID)  
      #print prot_geneID
#print len(prot_geneID_list)
#print prot_geneID_list [1]

##### RNAseq data ######################################################
######################################################################
f_RNAseq = csv.reader(open('LM_RNAseq2013_Rastrojo.csv', 'rb'), delimiter=',')
a=0
RNAseq_geneID_list = []
RNAseq_transcrptID_geneID_FPKM = dict()
for row in f_RNAseq:
  a=a+1
  if a>1:
    RNAseq_TranscriptID = row[0]
    RNAseq_geneID = row [3]
    RNAseq_FPKM = row [8]
    key = RNAseq_TranscriptID +"_"+RNAseq_geneID
    
    RNAseq_geneID_list.append(RNAseq_geneID)
    RNAseq_transcrptID_geneID_FPKM[key] = float(RNAseq_FPKM)
#print len(RNAseq_geneID_list)
#print len(RNAseq_transcrptID_geneID_FPKM.keys())


##### iAC560 model  ######################################################
######################################################################

f_iAC560 = csv.reader(open('iAC560_gene_pathways.csv', 'rb'), delimiter=',')
a=0
iAC560_geneID_list = []
iAC560_pathways_geneID = dict() 
for row in f_iAC560:
  a=a+1
  if a>1:
    f_iAC560_pathway = row[0]
    gene_list = row[1:]
    ## creating gene list ################
    for genes in gene_list:
      if genes !="":
        if genes not in iAC560_geneID_list:
          iAC560_geneID_list.append(genes)
     ###############################     
    while "" in gene_list :
      gene_list.remove("")

    ## creating pathways and associated gene dict #########  
    if len(gene_list) >0:
      for item in gene_list:
        if f_iAC560_pathway in iAC560_pathways_geneID:    
          iAC560_pathways_geneID[f_iAC560_pathway].append(item)
        else:
          iAC560_pathways_geneID[f_iAC560_pathway] =  [item]
    else:
      iAC560_pathways_geneID[f_iAC560_pathway] = ""
    ##############################################
      
#print iAC560_pathways_geneID.keys()[7]
#print iAC560_pathways_geneID[iAC560_pathways_geneID.keys()[7]]
#print "no of pathways in iAC560 =", len(iAC560_pathways_geneID.keys())
#print "no of pathways in iAC560 =", len(iAC560_geneID_list)
    

############### common genes studied  ########################################
print "prot_geneID_list =", len(prot_geneID_list)
print "RNAseq_geneID_list =", len(RNAseq_geneID_list)
print "iAC560_geneID_list =", len(iAC560_geneID_list)

RNSseq_proteomic = list(set(RNAseq_geneID_list).intersection(prot_geneID_list))
print "RNSseq_proteomic = ", len(RNSseq_proteomic)

RNSseq_iAC560 = list(set(RNAseq_geneID_list).intersection(iAC560_geneID_list))
print "RNSseq_iAC560 =", len(RNSseq_iAC560)

proteomics_iAC560 = list(set(prot_geneID_list).intersection(iAC560_geneID_list))
print "proteomics_iAC560 =",len(proteomics_iAC560)

s=[]
s.append(prot_geneID_list)
s.append(RNAseq_geneID_list)
s.append(iAC560_geneID_list)
proteomics_RNAseq_iAC560 = set.intersection(*map(set,s))
print "proteomics_RNAseq_iAC560 =", len(proteomics_RNAseq_iAC560)

#########################################################################




import matplotlib.pyplot as plt
import numpy as np
from scipy.stats import gaussian_kde
data = RNAseq_transcrptID_geneID_FPKM.values()

print len(data)
density = gaussian_kde(data)
xs = np.linspace(0,400,200)
density.covariance_factor = lambda : .25
density._compute_covariance()
plt.plot(xs,density(xs))
#plt.show()



o=open('out.txt','w')
for i in iAC560_geneID_list:
  o.write(i+'\n')
#print 
#print 

