import timeit
import mysql.connector
import xml.etree.ElementTree as ET
import requests
import re
import math
import nltk 
import pandas as pd
import traceback

start = timeit.default_timer()
print("start_time:", start)

mydb = mysql.connector.connect(host="",user="",password="",database="") #Database connection
mycursor = mydb.cursor() 

mycursor.execute("TRUNCATE TABLE ") # Truncate old BIONDA 
#mycursor.execute("DROP TABLE sen_v1")
df_disease = pd.read_sql_query("select ID,Name from ", mydb) # Select disease from DB


# Select Marker: Gene/Proteins, mirna and lncRNA from MySql Dictionaries
df = pd.read_sql_query("select * from ", mydb)

df_mirna =pd.read_sql_query("Select * from ;", mydb)

df_lncrrna = pd.read_sql_query("SELECT * FROM;", mydb)

#normalize Genenames
df11=df.assign(CA=df['CA'].str.split(' ')).explode('CA')

mask=(df11['CA'].str.len()>2)

df_biomarker=df11.loc[mask].reset_index()

#Url for EuropePMC
print("*********************************************************************************************************************************************\n")
urla='https://www.ebi.ac.uk/europepmc/webservices/rest/search?query=(%22biomarker%22%20OR%20%22biomarkers%22%20OR%20%22biological%20marker%22%20OR%20%22biological%20markers%22)%20AND%20(SRC%3A%22AGR%22%20OR%20SRC%3A%22CBA%22%20OR%20SRC%3A%22CTX%22%20OR%20SRC%3A%22ETH%22%20OR%20SRC%3A%22HIR%22%20OR%20SRC%3A%22MED%22%20SRC%3A%22NBK%22%20OR%20SRC%3A%22PAT%22%20OR%20SRC%3A%22PMC%22)%20AND%20(LANG%3A%22eng%22%20OR%20LANG%3A%22en%22%20OR%20LANG%3A%22us%22)%20AND%20(HAS_ABSTRACT%3Ay)&resultType=idlist&pageSize=1000&format=xml'
urlb='https://www.ebi.ac.uk/europepmc/webservices/rest/search?query=(%22biomarker%22%20OR%20%22biomarkers%22%20OR%20%22biological%20marker%22%20OR%20%22biological%20markers%22)%20AND%20(SRC%3A%22AGR%22%20OR%20SRC%3A%22CBA%22%20OR%20SRC%3A%22CTX%22%20OR%20SRC%3A%22ETH%22%20OR%20SRC%3A%22HIR%22%20OR%20SRC%3A%22MED%22%20SRC%3A%22NBK%22%20OR%20SRC%3A%22PAT%22%20OR%20SRC%3A%22PMC%22)%20AND%20(LANG%3A%22eng%22%20OR%20LANG%3A%22en%22%20OR%20LANG%3A%22us%22)%20AND%20(HAS_ABSTRACT%3Ay)&resultType=idlist&pageSize=1000&format=xml'

#download Paper IDs
array1 = []
re1=requests.get(urla)
root = ET.fromstring(re1.content)
for hitCount in root.iter('hitCount'):
    hit_count=int(hitCount.text)
result_value=hit_count-1

print("total papers:", result_value)

hit_count1=hit_count/1000
hit_count=math.ceil(hit_count1)
print("No. of API calls for saving paper_id's:", hit_count)

print("Downloading paper id's......")

for x in range(hit_count):
    re1=requests.get(urla)
    root1 = ET.fromstring(re1.content)
    print("api request:", x)
    for pid in root1.iter('id'):
        array1.append(pid.text)
    for nextCursorMark in root1.iter('nextCursorMark'):
        urla=urlb
        urla =urla+"&cursorMark="+nextCursorMark.text

#Start to download all Abstracts         
print("Parsing started......")
          
for i in range(result_value):
    
    paper_id=(array1[i]) 
    #download single abstract
    make_url= 'https://www.ebi.ac.uk/europepmc/webservices/rest/search?query=ext_id:'+paper_id+'&resultType=core&format=xml'
    re2=requests.get(make_url)
    try:
        root2 = ET.fromstring(re2.content)
    
        Jornal='NA'
         #extract Metadata 
        if root2.find('.//journal'):
            #print(root2.find('.//journal'))
            for journal in root2.find('.//journal'):
                if(journal.tag=='title'):
                    Jornal=journal.text
        for cit_value in root2.iter('citedByCount'):
            citation=int(cit_value.text)
            
        if paper_id.startswith("PPR"):
            paper_type="Abstract / Preprint"
        else:
            paper_type="Abstract"
     
        
        for authorString in root2.iter('authorString'):
            Autor=authorString.text
        for firstPublicationDate in root2.iter('firstPublicationDate'):
            PublicationDate=firstPublicationDate.text
        for abstractText in root2.iter('abstractText'): #NLK tokenize to sentence
            abstract_text_without_tags= re.sub(r"<[^>]*>"," ",abstractText.text )#remove existing <> tags from text
            nltk_tokens = nltk.sent_tokenize(abstract_text_without_tags)#tokenize the text
            #print(nltk_tokens)
            #normalisiierung für sen 
            for sen in nltk_tokens: #normalisition of abstracts
                sen_N = sen.replace("[", "")
                sen_N = sen_N.replace("]", "")
                sen_N = sen_N.replace("(", "")
                sen_N = sen_N.replace(")", "")
                sen_N = sen_N.replace(",", " ")
                sen_N = sen_N.replace("\'", "")
                sen_N = sen_N.replace("-", "")
                sen_N = sen_N.replace("α", "alpha")
                sen_N = sen_N.replace("β", "beta")     
                #print(sen_N)
                newlist=[[df_disease.at[dis_index, 'Name'],df_disease.at[dis_index, 'ID'],paper_id,citation,Jornal, Autor, PublicationDate,sen,sen_N] for dis_index in range(len(df_disease))  if df_disease.at[dis_index,'Name'].lower() in sen.lower() if(re.search(rf"\b{re.escape(df_disease.at[dis_index,'Name'])}\b", sen_N, re.IGNORECASE))] 
                            
    
                if newlist:
                    #check for genes/proteins from gene_df
                    for gene_index in range(len(df_biomarker)):
                        biomarker=df_biomarker.at[gene_index, 'CA']
                        proteinName = df_biomarker.at[gene_index, 'ProteinName']
                        protein = df_biomarker.at[gene_index, 'Protein']
                        if(biomarker in sen_N):
                            if(re.search(rf"\b{re.escape(biomarker)}\b",sen_N)):
                                query_string = #Insert in your final Table
                                for zzz in range(len(newlist)):
                                    mycursor.execute(query_string, newlist[zzz]+[biomarker]+[df_biomarker.at[gene_index, 'PID']]+[df_biomarker.at[gene_index, 'ProteinName']]+[df_biomarker.at[gene_index, 'Protein']]+[paper_type])
                            if(re.search(rf"\b{re.escape(proteinName)}\b",sen_N)):
                                	query_string = #Insert in your final Table                                for zzz in range(len(newlist)):
                                    mycursor.execute(query_string, newlist[zzz]+[biomarker]+[df_biomarker.at[gene_index, 'PID']]+[df_biomarker.at[gene_index, 'ProteinName']]+[df_biomarker.at[gene_index, 'Protein']]+[paper_type])
                            if(re.search(rf"\b{re.escape(protein)}\b",sen_N)):
                                    query_string = #Insert in your final Table
                                for zzz in range(len(newlist)):
                                    mycursor.execute(query_string, newlist[zzz]+[biomarker]+[df_biomarker.at[gene_index, 'PID']]+[df_biomarker.at[gene_index, 'ProteinName']]+[df_biomarker.at[gene_index, 'Protein']]+[paper_type])
                   
                           
                    for mirna_index in range(len(df_mirna)):
                        mirna_name=df_mirna.at[mirna_index, 'miRNAname']
                        if(mirna_name.lower() in sen_N.lower()):
                            if(re.search(rf"\b{re.escape(mirna_name)}\b",sen_N,re.IGNORECASE)):
                                query_string = #Insert in your final Table
                                for ttt in range(len(newlist)):
                                    mycursor.execute(query_string, newlist[ttt]+[mirna_name]+[df_mirna.at[mirna_index, 'idmiRNA']]+[df_mirna.at[mirna_index, 'ProteinName']]+[df_mirna.at[mirna_index, 'Protein_Entry']]+[paper_type])

      
                    for lncRNAindex in range(len(df_lncrrna)):
                        lcnrna_name=df_lncrrna.at[lncRNAindex, 'gene_ID']
                        lncrna_trans=df_lncrrna.at[lncRNAindex, 'trans_ID']
                        if(lcnrna_name.lower() in sen_N.lower()):
                            if(re.search(rf"\b{re.escape(lcnrna_name)}\b",sen_N,re.IGNORECASE)):
                            	query_string = #Insert in your final Table
                                for ttt in range(len(newlist)):
                                    mycursor.execute(query_string, newlist[ttt]+[lcnrna_name]+[df_lncrrna.at[lncRNAindex, 'trans_ID']]+[df_lncrrna.at[lncRNAindex, 'ProteinName']]+[df_lncrrna.at[lncRNAindex, 'Protein_Entry']]+[paper_type])
                        if(lncrna_trans.lower() in sen_N.lower()):
                            if(re.search(rf"\b{re.escape(lncrna_trans)}\b",sen_N,re.IGNORECASE)):
                                query_string = #Insert in your final Table
                            	for ttt in range(len(newlist)):
                                	mycursor.execute(query_string, newlist[ttt]+[lcnrna_name]+[df_lncrrna.at[lncRNAindex, 'trans_ID']]+[df_lncrrna.at[lncRNAindex, 'ProteinName']]+[df_lncrrna.at[lncRNAindex, 'Protein_Entry']]+[paper_type])
                     
    except Exception:
        #print("Some Problem occurred for paper no",i,"with paper_id:", paper_id)
        traceback.print_exc()
        pass

#new script for double


#Delete Duplicates
delete_create= "CREATE TABLE biomarker.temp ..;"
mycursor.execute(delete_create,mydb)
delete_inset ="INSERT INTO temp  SELECT DISTINCT "
mycursor.execute(delete_inset,mydb)
delete_drop ="DROP TABLE "
mycursor.execute(delete_drop,mydb)
delete_rename ="RENAME TABLE temp"
mycursor.execute(delete_rename,mydb)

#Delete cancer/carcrinoma in case of double tag
doubled =""
mycursor.execute(doubled,mydb)
double_delete = "
mycursor.execute(double_delete,mydb)
mydb.commit()
mydb.close()   
   
stop = timeit.default_timer()  
print('Time: ', stop - start)            
  
    
    
