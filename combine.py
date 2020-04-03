def usage():
 print("Available command-line options:\n  -i PDBfile\n  -c Chain of mutation\n  -r Resid of the mutation\n  -w One letter of wild amino acid\n  -m One letter of mutatant amino acid\n  -f Input file\n  -o Output file\n  -d Model type\n  -h Display this command-line summary")
def net_volume(wild,mutation):
 vol = {'A':'88.6','R':'173.4','N':'114.1','D':'111.1','C':'108.5','E':'138.4','Q':'143.8','G':'60.1','H':'153.2','I':'166.7','L':'166.7','K':'168.6','M':'162.9','F':'189.9','P':'112.7','S':'89.0','T':'116.1','W':'227.8','Y':'193.6','V':'140.0'}
 return ('{:.1f}'.format(float(vol.get(mutation,'0'))-float(vol.get(wild,'0'))))
def net_hydrophobicity(wild,mutation):
 hyd = {'A':'0','R':'3.71','N':'3.47','D':'2.95','C':'0.49','E':'1.64','Q':'3.01','G':'1.72','H':'4.76','I':'-1.56','L':'-1.81','K':'5.39','M':'-0.76','F':'-2.2','P':'-1.52','S':'1.83','T':'1.78','W':'-0.38','Y':'-1.09','V':'-0.78'}
 return ('{:.1f}'.format(float(hyd.get(mutation,'0'))-float(hyd.get(wild,'0'))))
def flexibility(wild,mutation):
 flex = {'A':'1','R':'81','N':'36','D':'18','C':'3','E':'54','Q':'108','G':'1','H':'36','I':'9','L':'9','K':'81','M':'27','F':'18','P':'2','S':'3','T':'3','W':'36','Y':'18','V':'3'}
 return (int(flex.get(mutation,'0'))-int(flex.get(wild,'0')))
def mutation_hydrophobicity(wild,mutation):
 if wild in ('I','V','L','F','C','M','A','W'):
  if mutation in ('G','T','S','Y','P','H'):
   return 0
  if mutation in ('N','D','Q','E','K','R'):
   return 1
  if mutation in ('I','V','L','F','C','M','A','W'):
   return 2
 if wild in ('G','T','S','Y','P','H'):
  if mutation in ('G','T','S','Y','P','H'):
   return 3
  if mutation in ('N','D','Q','E','K','R'):
   return 4
  if mutation in ('I','V','L','F','C','M','A','W'):
   return 5
 if wild in ('N','D','Q','E','K','R'):
  if mutation in ('G','T','S','Y','P','H'):
   return 6
  if mutation in ('N','D','Q','E','K','R'):
   return 7
  if mutation in ('I','V','L','F','C','M','A','W'):
   return 8
def mutation_polarity(wild,mutation):
 if wild in ('R','H','K'):
  if mutation in ('A','C','G','I','L','M','F','P','W','V'):
   return 0
  if mutation in ('R','H','K'):
   return 1
  if mutation in ('N','Q','S','T','Y'):
   return 2
  if mutation in ('D','E'):
   return 3
 if wild in ('A','C','G','I','L','M','F','P','W','V'):
  if mutation in ('A','C','G','I','L','M','F','P','W','V'):
   return 4
  if mutation in ('R','H','K'):
   return 5
  if mutation in ('N','Q','S','T','Y'):
   return 6
  if mutation in ('D','E'):
   return 7
 if wild in ('N','Q','S','T','Y'):
  if mutation in ('A','C','G','I','L','M','F','P','W','V'):
   return 8
  if mutation in ('R','H','K'):
   return 9
  if mutation in ('N','Q','S','T','Y'):
   return 10
  if mutation in ('D','E'):
   return 11
 if wild in ('D','E'):
  if mutation in ('A','C','G','I','L','M','F','P','W','V'):
   return 12
  if mutation in ('R','H','K'):
   return 13
  if mutation in ('N','Q','S','T','Y'):
   return 14
  if mutation in ('D','E'):
   return 15
def mutation_size(wild,mutation):
 if wild in ('G','A','S'):
  if mutation in ('C','D','P','N','T'):
   return 0
  if mutation in ('E','V','Q','H'):
   return 1
  if mutation in ('M','I','L','K','R'):
   return 2
  if mutation in ('F','Y','W'):
   return 3
  if mutation in ('G','A','S'):
   return 4
 if wild in ('C','D','P','N','T'):
  if mutation in ('C','D','P','N','T'):
   return 5
  if mutation in ('E','V','Q','H'):
   return 6
  if mutation in ('M','I','L','K','R'):
   return 7
  if mutation in ('F','Y','W'):
   return 8
  if mutation in ('G','A','S'):
   return 9
 if wild in ('E','V','Q','H'):
  if mutation in ('C','D','P','N','T'):
   return 10
  if mutation in ('E','V','Q','H'):
   return 11
  if mutation in ('M','I','L','K','R'):
   return 12
  if mutation in ('F','Y','W'):
   return 13
  if mutation in ('G','A','S'):
   return 14
 if wild in ('M','I','L','K','R'):
  if mutation in ('C','D','P','N','T'):
   return 15
  if mutation in ('E','V','Q','H'):
   return 16
  if mutation in ('M','I','L','K','R'):
   return 17
  if mutation in ('F','Y','W'):
   return 18
  if mutation in ('G','A','S'):
   return 19
 if wild in ('F','Y','W'):
  if mutation in ('C','D','P','N','T'):
   return 20
  if mutation in ('E','V','Q','H'):
   return 21
  if mutation in ('M','I','L','K','R'):
   return 22
  if mutation in ('F','Y','W'):
   return 23
  if mutation in ('G','A','S'):
   return 24
def mutation_hbonds(wild,mutation):
 if wild in ('R','W','K'):
  if mutation in ('A','C','G','I','L','M','F','P','V'):
   return 0 
  if mutation in ('R','W','K'):
   return 1   
  if mutation in ('N','Q','S','T','H','Y'):
   return 2   
  if mutation in ('D','E'):
   return 3   
 if wild in ('A','C','G','I','L','M','F','P','V'):
  if mutation in ('A','C','G','I','L','M','F','P','V'):
   return 4  
  if mutation in ('R','W','K'):
   return 5   
  if mutation in ('N','Q','S','T','Y','H'):
   return 6   
  if mutation in ('D','E'):
   return 7   
 if wild in ('N','Q','S','T','Y','H'):
  if mutation in ('A','C','G','I','L','M','F','P','V'):
   return 8   
  if mutation in ('R','W','K'):
   return 9   
  if mutation in ('N','Q','S','T','Y','H'):
   return 10   
  if mutation in ('D','E'):
   return 11  
 if wild in ('D','E'):
  if mutation in ('A','C','G','I','L','M','F','P','V'):
   return 12   
  if mutation in ('R','W','K'):
   return 13   
  if mutation in ('N','Q','S','T','Y','H'):
   return 14   
  if mutation in ('D','E'):
   return 15
def mutation_chemical(wild,mutation):
 if wild in ('A','G','I','L','P','V'):
  if mutation in ('C','R','H','K'):
   return 0   
  if mutation in ('N','Q'):
   return 1   
  if mutation in ('D','E'):
   return 2   
  if mutation in ('C','M'):
   return 3   
  if mutation in ('S','T'):
   return 4   
  if mutation in ('F','W','Y'):
   return 5   
  if mutation in ('A','G','I','L','P','V'):
   return 6
 if wild in ('C','R','H','K'):
  if mutation in ('C','R','H','K'):
   return 7   
  if mutation in ('N','Q'):
   return 8   
  if mutation in ('D','E'):
   return 9  
  if mutation in ('C','M'):
   return 10   
  if mutation in ('S','T'):
   return 11   
  if mutation in ('F','W','Y'):
   return 12   
  if mutation in ('A','G','I','L','P','V'):
   return 13
 if wild in ('N','Q'):
  if mutation in ('C','R','H','K'):
   return 14   
  if mutation in ('N','Q'):
   return 15   
  if mutation in ('D','E'):
   return 16   
  if mutation in ('C','M'):
   return 17   
  if mutation in ('S','T'):
   return 18   
  if mutation in ('F','W','Y'):
   return 19   
  if mutation in ('A','G','I','L','P','V'):
   return 20
 if wild in ('D','E'):
  if mutation in ('C','R','H','K'):
   return 21   
  if mutation in ('N','Q'):
   return 22   
  if mutation in ('D','E'):
   return 23   
  if mutation in ('C','M'):
   return 24   
  if mutation in ('S','T'):
   return 25   
  if mutation in ('F','W','Y'):
   return 26   
  if mutation in ('A','G','I','L','P','V'):
   return 27   
 if wild in ('C','M'):
  if mutation in ('C','R','H','K'):
   return 28   
  if mutation in ('N','Q'):
   return 29   
  if mutation in ('D','E'):
   return 30   
  if mutation in ('C','M'):
   return 31   
  if mutation in ('S','T'):
   return 32   
  if mutation in ('F','W','Y'):
   return 33   
  if mutation in ('A','G','I','L','P','V'):
   return 34      
 if wild in ('S','T'):
  if mutation in ('C','R','H','K'):
   return 35   
  if mutation in ('N','Q'):
   return 36   
  if mutation in ('D','E'):
   return 37   
  if mutation in ('C','M'):
   return 38   
  if mutation in ('S','T'):
   return 39   
  if mutation in ('F','W','Y'):
   return 40   
  if mutation in ('A','G','I','L','P','V'):
   return 41   
 if wild in ('F','W','Y'):
  if mutation in ('C','R','H','K'):
   return 42   
  if mutation in ('N','Q'):
   return 43   
  if mutation in ('D','E'):
   return 44   
  if mutation in ('C','M'):
   return 45  
  if mutation in ('S','T'):
   return 46   
  if mutation in ('F','W','Y'):
   return 47   
  if mutation in ('A','G','I','L','P','V'):
   return 48
def mutation_ala(wild,mutation):
 if wild=='A' or mutation=='A':
  return 1
 else:
  return 0        
def translate_aa(three_letter):
 trans = {'ALA':'A','ARG':'R','ASN':'N','ASP':'D','CYS':'C','GLU':'E','GLN':'Q','GLY':'G','HIS':'H','ILE':'I','LEU':'L','LYS':'K','MET':'M','PHE':'F','PRO':'P','SER':'S','THR':'T','TRP':'W','TYR':'Y','VAL':'V'}
 return (trans[three_letter])
def mutation_aa_label(three_letter):
 aa = {'ALA':'1','ARG':'2','ASN':'3','ASP':'4','CYS':'5','GLU':'6','GLN':'7','GLY':'8','HIS':'9','ILE':'10','LEU':'11','LYS':'12','MET':'13','PHE':'14','PRO':'15','SER':'16','THR':'17','TRP':'18','TYR':'19','VAL':'20'}
 return (aa.get(three_letter,0))
def mutation_type(wild,mutation):
 wild_lists = ['A','F','C','D','N','E','Q','G','H','L','I','K','M','P','R','S','T','V','W','Y']
 mutation_lists = ['A','F','C','D','N','E','Q','G','H','L','I','K','M','P','R','S','T','V','W','Y']
 label_1=0
 label_2=0
 for i in wild_lists:
  for j in mutation_lists:
   if i != j:
    label_1 += 1
    if wild == i:
     if mutation == j and i != j:
      label_2 = label_1
      break
 return label_2
#def mutation_pdb(pdbid):
# from Bio.PDB import PDBList
# pdbl = PDBList()
# pdbl.retrieve_pdb_file(pdbid,file_format='pdb',pdir='.',overwrite=True)
#def mutation_pdb(pdbid):
# null_value=os.system('wget \'https://files.rcsb.org/view/'+pdbid+'.pdb\' -O pdb'+pdbid+'.ent -q')
def mutation_pdb(pdbid):
 import urllib
 url="https://files.rcsb.org/view/"+pdbid+".pdb"
 urllib.urlretrieve(url, pdbid)
def mutation_sequence(pdbid,resid,chain):
 label_index=1
 resid_label=[]
 resid_label_aa=[]
 mutation_coordinate=[]
 resid_label_aa=[0 for _ in range(11)] #how many sequence to use
 for i in range(len(resid_label_aa)):
  resid_int=filter(lambda ch: ch in '-0123456789', resid.strip())
  resid_label.append(int(resid_int)-(len(resid_label_aa)-1)/2+i)
  for line in open(pdbid):
   pdbstr=line.strip()
   if pdbstr[0:4]=="ATOM":
    if pdbstr[21:22]==chain:
     if pdbstr[22:27].strip()==str(resid_label[i]) or (pdbstr[22:27].strip()==str(resid) and label_index==1):
      if pdbstr[13:15]=="CA":
       if pdbstr[22:27].strip()==str(resid):
        mutation_coordinate=[float(pdbstr[29:38].strip()),float(pdbstr[38:46].strip()),float(pdbstr[46:55].strip())]
        label_index=0
       resid_label_aa[i]=pdbstr[17:20]
       break
 if len(mutation_coordinate) !=3:
  error_index=1
 else:
  error_index=0
 return resid_label_aa,mutation_coordinate,error_index
def mutation_distance(pdbid,chain,mutation_coordinate):
 resid_label_dis_aa=[]
 resid_label_aa=[]
 resid_label_distance=[]
 for line in open(pdbid):
   pdbstr=line.strip()
   if pdbstr[0:4]=="ATOM":
    if pdbstr[13:15]=="CA":
     if pdbstr[21:22]!=chain:
      mutation_coordinate1=[float(pdbstr[29:38].strip()),float(pdbstr[38:46].strip()),float(pdbstr[46:55].strip())]
      resid_label_aa.append(pdbstr[17:20])
      resid_label_distance.append(np.sqrt(np.square(mutation_coordinate[0]-mutation_coordinate1[0])+np.square(mutation_coordinate[1]-mutation_coordinate1[1])+np.square(mutation_coordinate[2]-mutation_coordinate1[2])))      
 b=zip(resid_label_distance,range(len(resid_label_distance)))
 b.sort(key = lambda x : x[0])
 c = [x[1] for x in b]
 sequ_num=10 #10 sequence in total,use last 7sequence
 if len(c)>=sequ_num:
  for j in range(0,sequ_num,1):
   if b[j][0]<=10: 
    resid_label_dis_aa.append(resid_label_aa[c[j]])
   else:
    resid_label_dis_aa.append(0)
 else:
  for j in range(0,len(c),1): 
   if b[j][0]<=10: 
    resid_label_dis_aa.append(resid_label_aa[c[j]])
   else:
    resid_label_dis_aa.append(0)
  while len(resid_label_dis_aa) < sequ_num:
   resid_label_dis_aa.append(0)
 return  resid_label_dis_aa
def mutation_pdb_information(pdbid):
 reso=0
 r_value=0
 temp=0
 ph=0
 for line in open(pdbid):
  pdbstr=line.strip()
  if pdbstr[0:22]=="REMARK   2 RESOLUTION.":
   try:
    reso=float(pdbstr[26:30].strip())
   except ValueError:
    reso=0
  if pdbstr[0:45]=="REMARK   3   R VALUE            (WORKING SET)":
   try:
    r_value=float(pdbstr[49:54].strip())
   except ValueError:
    r_value=0
  if pdbstr[0:23]=="REMARK 200  TEMPERATURE":
   try:
    temp=float(pdbstr[45:48].strip())
   except ValueError:
    temp=0
  if pdbstr[0:14]=="REMARK 200  PH":
   ph=filter(lambda ch: ch in '0123456789.', pdbstr[45:48].strip())
   try:
    ph=float(ph)
   except ValueError:
    ph=0
   break
 return reso,r_value,temp,ph
def file_loop(pdb_id,mutation_chain,mutation_resid,wild_aa,mutation_aa,model_type):
 delete_index=0
 if len(wild_aa) == 3:
  wild_aa=translate_aa(wild_aa)
 if len(mutation_aa) == 3:
  mutation_aa=translate_aa(mutation_aa)
 label = []
 label.append(net_volume(wild_aa,mutation_aa))
 label.append(net_hydrophobicity(wild_aa,mutation_aa))
 label.append(flexibility(wild_aa,mutation_aa))
 label.append(mutation_hydrophobicity(wild_aa,mutation_aa))
 label.append(mutation_polarity(wild_aa,mutation_aa))
 label.append(mutation_type(wild_aa,mutation_aa))
 label.append(mutation_size(wild_aa,mutation_aa))
 label.append(mutation_hbonds(wild_aa,mutation_aa))
 label.append(mutation_chemical(wild_aa,mutation_aa))
 #label.append(mutation_ala(wild_aa,mutation_aa))
 if os.path.isfile(pdb_id) == 0: 
  #mutation_pdb(pdb_id)
  #pdb_id="pdb"+pdb_id+".ent"
  return('Please check the PDB file')
 (resid_label_aa,mutation_coordinate,error_index)=mutation_sequence(pdb_id,str(mutation_resid),mutation_chain)
 if error_index == 1:
  return('Please check the PDB file, there is no coordinate at the mutation site')
 label_aa_distance=mutation_distance(pdb_id,mutation_chain,mutation_coordinate)
 for i in resid_label_aa:
  label.append(mutation_aa_label(i))
 for i in label_aa_distance:
  label.append(mutation_aa_label(i))
 (reso,r_value,temp,ph)=mutation_pdb_information(pdb_id)
 label.append(reso)
 label.append(temp)
 label.append(ph)
 #if delete_index ==1:
  #os.remove(pdb_id)
 return(pred_feature(label,model_type))
def pred_feature(label,model_type):
 if str(model_type) == '1':
  model=xgb.Booster(model_file='regression.model')
 else:
  model=xgb.Booster(model_file='classification.model')
 x=np.array(label)
 #x=pd.DataFrame(x)
 #x=x.iloc[:,[-1]].values
 x=x.reshape((1,len(label)))
 x=xgb.DMatrix(x)
 y_pred = model.predict(x)
 if str(model_type) == '1':
  if y_pred[0]>0:
   return ("%.2f"%y_pred[0]+" Destabilizing")
  else:
   return ("%.2f"%y_pred[0]+" Stabilizing")
 else:
  if y_pred[0]>0.5:
   return('Disruptive')
  else:
   return('Nondisruptive') 
 #return (y_pred[0])
def output_format(model):
 if str(model)=="1":
  return("Structure_file Chain Position Wild Mutant ddG(kcal/mol) Type")
 else:
  return("Structure_file Chain Position Wild Mutant Type")
import sys, getopt
import os
import numpy as np
import xgboost as xgb
#import pandas as pd
opts, args = getopt.getopt(sys.argv[1:], "hi:c:r:w:m:f:d:o:", ["help=", "pdbfile=","chain=","resid=","wild=","mutation_list","model","output="])
for op, value in opts:
 if op == "-i":
  pdb_id = value
  #pdb_id = str.lower(pdb_id)
 elif op == "-c":
  mutation_chain = value
 elif op == "-r":
  mutation_resid = value
 elif op == "-w":
  wild_aa= value
 elif op == "-m":
  mutation_aa = value
 elif op == "-d":
  model_type = value
 elif op == "-o":
  output_name = value
 elif op == "-f":
  file_name = value
 elif op == "-h":
  usage()
  sys.exit()
try:
 file_name
except NameError:
 try:
  output_name
 except NameError:
  pred=file_loop(pdb_id,mutation_chain,mutation_resid,wild_aa,mutation_aa,model_type)
  print pred
 else:
  f = open(output_name, 'w')
  pred=file_loop(pdb_id,mutation_chain,mutation_resid,wild_aa,mutation_aa,model_type)
  print >>f,pred
  f.close()
else:
 try:
  output_name
 except NameError:
  for line in open(file_name):
   info=line.strip().split(' ')
   pred=file_loop(pdb_id,info[0],info[1],str(info[2]),info[3],model_type)
   print pred
  sys.exit()
 else:
  f = open(output_name, 'w')
  print >>f,output_format(model_type)
  for line in open(file_name):
   info=line.strip().split(' ')
   pred=file_loop(pdb_id,info[0],info[1],str(info[2]),info[3],model_type)
   print >>f,pdb_id,info[0],info[1],str(info[2]),info[3],pred
  f.close()
#pred=file_loop(pdb_id,mutation_chain,mutation_resid,wild_aa,mutation_aa)
