import multiprocessing
from multiprocessing import Pool
from datetime import datetime
from elasticsearch import Elasticsearch
import sys
import numpy as np
import pandas as pd
from sklearn.externals import joblib
from sklearn.preprocessing import StandardScaler
import elasticsearch.helpers
from io import StringIO
import csv


es = Elasticsearch('http://borgdb.cbrc.kaust.edu.sa:9200/', timeout=30, max_retries=10, retry_on_timeout=True)
es.indices.put_settings(index = "pubmed-2021.v2",
   body = {
      "index": {
         "max_result_window": 20000000
      }
   })

dic_hybrid_scores={}
dic_gene={}
dic_rs={}
dic_var={}
dic_uniprot={}
dic_hgvs={}
dic_sift={}
dic_vt={}
dic_polyphen={}
dic_af={}
var2pheno={}
hpo2label={}
dic_ppi={}
vcf={}
#------------------------------
var2pmid={}
hpo2pmid={}
gene2pmid={}
hgnc2uniprot={}

def load_hgnc():
    fname="hgnc2uniprot.txt"
    with open (fname, 'r') as f:
     for line in f:
       lst=line.strip().split("\t")
       if len(lst)==2:
           hgnc2uniprot[lst[0]]=list(lst[1].split("##"))
def load_variants():
    fname="out.rs_pmid.txt"
    with open (fname, 'r') as f:
     for line in f:
       lst=line.strip().split("\t")
       if len(lst)==2:
           var2pmid[lst[0]]=set(lst[1].split(","))
       else:
           var2pmid[lst[0]]=set(str(0))        

  
def load_phenotypes():
    fname="merged.hpo2pmid.txt"
    with open (fname, 'r') as f:
     for line in f:
       lst=line.strip().split("\t")
       if len(lst)==3:
           hpo2pmid[lst[0]]=set(lst[2].split(","))
       else:
           hpo2pmid[lst[0]]=set(str(0))

def load_genes():
  fname="out.uniprot_pmid.txt"
  with open (fname, 'r') as f:
       for line in f:
         lst=line.strip().split("\t")
         if len(lst)==2:
           gene2pmid[lst[0]]=set(lst[1].split(","))
         else:
           gene2pmid[lst[0]]=set(str(0))
   
  fname="out.gene.symbol_pmid.txt"
  with open (fname, 'r') as f:
       for line in f:
         lst=line.strip().split("\t")
         if len(lst)==2:
           gene2pmid[lst[0]]=set(lst[1].split(","))
         else:
           gene2pmid[lst[0]]=set(str(0))


def search(index,qry):
   r=[]
   qry=qry.strip()
   qry='"'+qry+'"'

   res = es.search(index=index,
                   track_total_hits= True,
                   _source=["_id"],
                    size=20000000,
                    body={
                        "query": {
                        "query_string": {
                            "query":qry }}},
                   request_timeout=(10*60))
   for hit in res['hits']['hits']:
      pid = hit['_id']
      r.append(pid)
   qry=qry.replace('"','')
   hpo2pmid[qry]=set(r)


def readfromfile(fname):
    vt_weights={}
###Uploaded_variation     Location        Allele  Gene    Feature Feature_type    Consequence     cDNA_position   CDS_position    Protein_position        Amino_acids     Codons  Existing_variation      IMPACT  DISTANCE        STRAND  FLAGS   SYMBOL  SYMBOL_SOURCE   HGNC_ID BIOTYPE CANONICAL       TSL     SOURCE  SIFT    PolyPhen        EXON    INTRON  HGVSc   HGVSp   HGVS_OFFSET     AF      gnomAD_AF       gnomAD_AFR_AF   gnomAD_AMR_AF   gnomAD_ASJ_AF   gnomAD_EAS_AF   gnomAD_FIN_AF   gnomAD_NFE_AF   gnomAD_OTH_AF   gnomAD_SAS_AF   CLIN_SIG        SOMATIC PHENO   PAVS    GO_CLASSES      PHENOTYPE       PPI
    reader=csv.reader(open("var_type_weights.txt"), delimiter="\t")
    for row in reader:
      vt_weights[row[0]]=float(row[1])

    with open(fname, 'r') as f:
       cnt=0
       for line in f:
         lst=line.strip().split("	")
         #gene_list=[]
         #uniprot_list=[]
         clist=[]  
         ulist=[]
         hgvslist=[]
         aflist=[]
         #consider canonical variants only
         if ((len(lst)>5) and (lst[0] not in "#Uploaded_variation") and (lst[21] == "YES") and ("protein_coding" in line)):
           cnt=cnt+1
           vcf[cnt]=line.strip()
           dic_var[cnt]=lst[1]+"_"+lst[2]
           tmp=lst[12].split(",")
           dic_rs[cnt]=tmp[0]
           dic_gene[cnt]=lst[17]
           tmp=lst[12].split(",")
           for t in tmp:
             if "rs" in t:
               dic_rs[cnt]=t  
           if lst[19] not in "-" and lst[19] in hgnc2uniprot:
             ulist=hgnc2uniprot[lst[19]]
             if len(ulist)>0:
               dic_uniprot[cnt]=ulist
           if "-" not in lst[6]:
             tmp=lst[6].split(",")
             dic_vt[cnt]=float(vt_weights[tmp[0]])
           if lst[28] not in "-":
             tmp=lst[28].split(":")
             if len(tmp)>1:
               hgvslist.append(tmp[1])
             else:
               hgvslist.append(tmp[0])
           if lst[29] not in "-":
             tmp=lst[29].split(":")
             if len(tmp)>1:
               hgvslist.append(tmp[1])
             else:
               hgvslist.append(tmp[0])

           dic_hgvs[cnt]=hgvslist
      
           if lst[24] not in "-":
              tmp=lst[24].split("(")
              tmp=tmp[1].split(")")
              dic_sift[cnt]=1.0-float(tmp[0])
           else:
              dic_sift[cnt]=float(0) 
 
           if lst[25] not in "-":
              tmp=lst[25].split("(")
              tmp=tmp[1].split(")")
              dic_polyphen[cnt]=float(tmp[0])
           else:
              dic_polyphen[cnt]=float(0)
           for i in range (32,40):
             if lst[i] not in "-":
               aflist.append(float(lst[i]))
             else:
               aflist.append(float(0))
           dic_af[cnt]=max(aflist)
           #PPI DATA
           if (len(lst)>46):
             if len(lst[47])>0:          
               tmp=lst[47].split("--")
               if len(tmp)>1:
                 dic_ppi[cnt]=tmp[1]        

def calc_jaccard_ppi (cnt, hpo_list):
 l=[]
 tmp_var=set()
 if cnt in dic_ppi:
   l=dic_ppi[cnt].split("##")
     
 if len(l)>0:
     #tmp_var=set()
     for qry in l:
        if qry in gene2pmid:
          tmp_var=tmp_var.union(gene2pmid[qry])

 if len(tmp_var)>0:
  res_lst=[]
  for pheno_label in hpo_list:
    if pheno_label in hpo2pmid and len(hpo2pmid[pheno_label])>0:

       res1=len(tmp_var.intersection(hpo2pmid[pheno_label]))
       res2=len(tmp_var.union(hpo2pmid[pheno_label]))

       if res2>0:
         res_lst.append(res1/res2)
       else:
         res_lst.append(0)
    else:
      res_lst.append(0)
 else: 
   jaccard=0
   return jaccard

 jaccard=sum(res_lst)/len(hpo_list)
 return jaccard

def calc_jaccard_prot (cnt,hpo_list):
  qry_var=[]
  gene=dic_gene[cnt]
  if len(gene)>2:  ## if length of gene symbol >2
     qry_var.append(gene)

  if cnt in dic_uniprot and len(dic_uniprot[cnt])>0:
     qry_var=qry_var+dic_uniprot[cnt]

  if len(qry_var)>0: 
     tmp_var=set()
     for qry in qry_var:
       if qry in gene2pmid:
         tmp_var=tmp_var.union(gene2pmid[qry])
  else:
    jaccard=0
    return jaccard

  res_lst=[]
  for pheno_label in hpo_list:
    if pheno_label in hpo2pmid and len(hpo2pmid[pheno_label])>0:

      res1=len(tmp_var.intersection(hpo2pmid[pheno_label]))
      res2=len(tmp_var.union(hpo2pmid[pheno_label]))

      if res2>0:
       res_lst.append(res1/res2)
      else:
       res_lst.append(0)
    else:
      res_lst.append(0)

  jaccard=sum(res_lst)/len(hpo_list)
  return jaccard


def calc_jaccard_rs (cnt, hpo_list):
 qry_var=[]
 rsid = dic_rs[cnt]
 if rsid not in "-" and rsid in var2pmid:
    qry_var.append(rsid)

 if cnt in dic_hgvs and len(dic_hgvs[cnt])>0:
   for item in dic_hgvs[cnt]:
      if item in var2pmid:
          qry_var.append(item)
 if len(qry_var)>0:
   tmp_var=set()
   for qry in qry_var:
     tmp_var=tmp_var.union(var2pmid[qry])
 else:
   jaccard=0
   return jaccard
 
 if len(hpo_list)<1:
   jaccard=0
   return jaccard

 res_lst=[]
 for pheno_label in hpo_list:
   if pheno_label in hpo2pmid and len(hpo2pmid[pheno_label])>0:
       res1=len(tmp_var.intersection(hpo2pmid[pheno_label]))       
       res2=len(tmp_var.union(hpo2pmid[pheno_label]))
       if res2>0:
         res_lst.append(res1/res2)
       else:
         res_lst.append(0)
   else:
       res_lst.append(0)
 jaccard=sum(res_lst)/len(hpo_list)
 return jaccard   

def get_vt(cnt,hpo_list):
  return dic_vt[cnt]

def get_poly(cnt,hpo_list):
  return dic_sift[cnt]

def get_sift(cnt,hpo_list):
  return dic_polyphen[cnt]
##-------------------------------------------------------------------------------
print ("vcf file="+sys.argv[1])
print ("preferred genomic evidence based model (SIFT/PloyPhen-2)="+sys.argv[2])
print ("user provided symptoms:"+sys.argv[3])
if sys.argv[1] =="":
   print ("provide an input file annotated by VEP in VCF format. USAGE: python STARVar.py inputfile option symptoms")
   sys.exit(0)
if (sys.argv[2] != "s" and sys.argv[2] !="p" and sys.argv[2] !="v"):
   print ("provide an option for genomic evidence based model: s for SIFT; p for PolyPhen-2; v for variant type. USAGE: python STARVAR.py inputfile option symptoms")
   sys.exit(0)

if sys.argv[3] =="":
   print ("provide a set of patient symptoms, use semicolon(;) to seperate them. USAGE: python STARVar.py inputfile option symptoms")
   sys.exit(0)
if "HP_" in sys.argv[3]:
   print ("HPO IDs should be in the form of HP:XXXXXX")
   sys.exit(0)

#hpo_labels=["HP:0001250","HP:0001249","HP:0002445","HP:0001263","HP:0000253"]
sys.argv[3]=sys.argv[3].replace('"','')
hpo_labels=list((sys.argv[3]).split(";"))

load_hgnc()

load_phenotypes()

for label in hpo_labels:
   if label not in hpo2pmid:
      search("pubmed-2021.v2",label)

load_variants()
load_genes()

infile=sys.argv[1]
readfromfile(infile)

f_rs=[]
f_prot=[]
f_ppi=[]
lst=[]
f_vt=[]
f_sift=[]
f_poly=[]

for cnt in range (1,(len(dic_var)+1)):
   lst.append((cnt,hpo_labels))

size=multiprocessing.cpu_count()-1
#print ("Nof cpus="+str(size))

process_pool = multiprocessing.Pool(size)

f_rs = process_pool.starmap(calc_jaccard_rs, lst)
f_prot=process_pool.starmap(calc_jaccard_prot, lst)
f_ppi=process_pool.starmap(calc_jaccard_ppi, lst)

if sys.argv[2] is "v":
  f_vt=process_pool.starmap(get_vt, lst)
if sys.argv[2] is "p":
  f_poly=process_pool.starmap(get_poly, lst)
if sys.argv[2] is "s":
  f_sift=process_pool.starmap(get_sift, lst)

if sys.argv[2] is "v": 
  testdata=pd.DataFrame({'Jaccard_rs': f_rs,
     'Jaccard_prot': f_prot,
     'Jaccard_ppi': f_ppi,
     'Variant_Cons':f_vt
    })
  testdata['Variant_Cons'] = testdata['Variant_Cons'].astype(float)
if sys.argv[2] is "p":
  testdata=pd.DataFrame({'Jaccard_rs': f_rs,
     'Jaccard_prot': f_prot,
     'Jaccard_ppi': f_ppi,
     'POLYPHEN2':f_poly
    })
  testdata['POLYPHEN2'] = testdata['POLYPHEN2'].astype(float)
if sys.argv[2] is "s":
  testdata=pd.DataFrame({'Jaccard_rs': f_rs,
     'Jaccard_prot': f_prot,
     'Jaccard_ppi': f_ppi,
     'SIFT':f_sift
    })
  testdata['SIFT'] = testdata['SIFT'].astype(float)

testdata['Jaccard_rs'] = testdata['Jaccard_rs'].astype(float)
testdata['Jaccard_prot'] = testdata['Jaccard_prot'].astype(float)
testdata['Jaccard_ppi'] = testdata['Jaccard_ppi'].astype(float) 
loaded_model = joblib.load('Lit.model.sav')

if sys.argv[2] is "v":
  ynew = loaded_model.predict_proba(testdata.drop(['Variant_Cons'], axis=1))
  testdata['LITCOMP']=ynew[:,1]
  ensemble_model = joblib.load('Lit.VC.model.sav')

if sys.argv[2] is "s":
  ynew = loaded_model.predict_proba(testdata.drop(['SIFT'], axis=1))
  testdata['LITCOMP']=ynew[:,1]
  ensemble_model = joblib.load('Lit.Sift.model.sav')

if sys.argv[2] is "p":
  ynew = loaded_model.predict_proba(testdata.drop(['POLYPHEN2'], axis=1))
  testdata['LITCOMP']=ynew[:,1]
  ensemble_model = joblib.load('Lit.Poly.model.sav')


testdata.drop(['Jaccard_rs', 'Jaccard_prot','Jaccard_ppi'], axis=1, inplace=True)

yens = ensemble_model.predict_proba(testdata)

for i in range(0,len(yens)):
   
  dic_hybrid_scores[vcf[(i+1)]]=float(yens[i][1])

newdic={k: v for k, v in sorted(dic_hybrid_scores.items(), key=lambda item: item[1], reverse=True)}
rank=0

for key in newdic:
  rank=rank+1
  print (key+"\t"+str(dic_hybrid_scores[key])+"\t"+str(newdic[key])+"\t"+str(rank))

now = datetime.now()
current_time = now.strftime("%H:%M:%S")

