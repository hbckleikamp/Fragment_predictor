# -*- coding: utf-8 -*-
"""
Created on Mon Aug  8 11:29:37 2022

@author: ZR48SA
"""

#%% modules
import numpy as np
import pandas as pd
import itertools
#%% Define masses

mass_H=1.00783
mass_C=12
mass_N=14.00307
mass_O=15.99491
mass_P=30.97376
mass_S=31.97207
mass_Se=79.91652
ele_mass=[mass_H,mass_C,mass_O,mass_N,mass_P,mass_S,mass_Se]
elements=["Hydrogen","Carbon","Oxygen","Nitrogen","Phosphor","Sulfur","Selenium"]

aa_molecular_formulas=pd.DataFrame([
#AA  H    C O N P S Se 
#["", 0,   0,0,0,0,0,0],
["A",5,   3,1,1,0,0,0],
["C",5,   3,1,1,0,1,0],
["D",5,   4,3,1,0,0,0],
["E",7,   5,3,1,0,0,0],
["F",9,   9,1,1,0,0,0],
["G",3,   2,1,1,0,0,0],
["H",7,   6,1,3,0,0,0],
["I",11,  6,1,1,0,0,0],
["J",11,  6,1,1,0,0,0],
["K",12,  6,1,2,0,0,0],
["L",11,  6,1,1,0,0,0],
["M",9,   5,1,1,0,1,0],
["N",6,   4,2,2,0,0,0],
["P",7,   5,1,1,0,0,0],
["Q",8,   5,2,2,0,0,0],
["R",12,  6,1,4,0,0,0],
["S",5,   3,2,1,0,0,0],
["T",7,   4,2,1,0,0,0],
["V",9,   5,1,1,0,0,0],
["W",10, 11,1,2,0,0,0],
["Y",9,   9,2,1,0,0,0],
# ["U",5,   3,1,1,0,0,1],
# ["O",19, 12,2,3,0,0,0], #no side chain neutral losses are described for U and O in expert system
#modified residues
["Mox", 9,5,2,1,0,1,0], #oxidized methionine
["Ccam",8,5,2,2,0,1,0], #carbamidomethylated cysteine "H(3) C(2) N O"
["Sph", 6,3,5,1,1,0,0], #phosphorylation "H O(3) P"
["Tph", 8,4,5,1,1,0,0],
["Yph",10,9,5,1,1,0,0]],
    columns=["AA"]+elements).set_index("AA")
#add acetylations #no side chain neutral losses are described for acetylated amino acids
# ac=aa_molecular_formulas+[2,2,1,0,0,0,0] # "H(2) C(2) O"
# ac.index=aa_molecular_formulas.index+"ac"
# aa_molecular_formulas=pd.concat([aa_molecular_formulas,ac])
aa_molecular_formulas=aa_molecular_formulas.T.to_dict(orient='list')


ion_molecular_formulas=pd.DataFrame([
#ion  H  C  O  N P S Se 
["a",-2,-1,-2, 0,0,0,0],
["b",-2, 0, 0, 0,0,0,0],
["c", 1, 0,-1, 1,0,0,0],
["x",-2, 1, 1, 0,0,0,0],
["y", 0, 0, 0, 0,0,0,0],
["z",-3, 0, 0,-1,0,0,0],  
["im",-2,-1, 0,-1,0,0,0]], #immonium
    columns=["ion"]+elements).set_index("ion")

neutral_loss_formulas=pd.DataFrame([
#loss      H    C   O   N  P  S Se 
["",       0,   0,  0,  0, 0, 0,0],
["H2O",   -2,   0, -1,  0, 0, 0,0], 
["H3PO4", -3,  -0, -4, -0,-1,-0,0], 
["HPO3",  -1,  -0, -3, -0,-1,-0,0],
["CH2S",  -2,  -1, -0, -0,-0,-1,0],
["H2O2",  -2,  -0, -2, -0,-0,-0,0],
["CO2",   -0,  -1, -2, -0,-0,-0,0],
["C2H4O2",-4,  -2, -2, -0,-0,-0,0],
["C2H4",  -4,  -2, -0, -0,-0,-0,0],
["C2H5N", -5,  -2, -0, -1,-0,-0,0], 
["C4H9N", -9,  -4, -0, -1,-0,-0,0], 
["C4H11N",-11, -4, -0, -1,-0,-0,0],
["C3H9N", -9,  -3, -0, -1,-0,-0,0],
["C3H6",  -6,  -3, -0, -0,-0,-0,0],
["C4H8",  -8,  -4, -0, -0,-0,-0,0],
["C2H4S", -4,  -2, -0, -0,-0,-1,0],
["C3H6S", -6,  -3, -0, -0,-0,-1,0],
["CH4SO", -4,  -1, -1, -0,-0,-1,0],
["C3H8SO",-8,  -3, -1, -0,-0,-1,0],
["C3H6SO",-6,  -3, -1, -0,-0,-1,0],
["NH3",   -3,  -0, -0, -1,-0,-0,0],
["CH3NO", -3,  -1, -1, -1,-0,-0,0],
["C2H5NO",-5,  -2, -1, -1,-0,-0,0],
["C3H5NO",-5,  -3, -1, -1,-0,-0,0],
["CH2N2", -2,  -1, -0, -2,-0,-0,0],
["C3H9N3",-9,  -3, -0, -3,-0,-0,0], 
["CH4O",  -4,  -1, -1, -0,-0,-0,0],
["C2H4O", -4,  -2, -1, -0,-0,-0,0],
["C3H6",  -6,  -3, -0, -0,-0,-0,0],
["C8H7N", -7,  -8, -0, -1,-0,-0,0],
["C9H9N", -9,  -9, -0, -1,-0,-0,0],
["CO",    -0,  -1, -1, -0,-0,-0,0],

#IONS
["a",     -2,  -1, -2,  0, 0, 0,0],
["b",     -2,   0,  0,  0, 0, 0,0],
["c",      1,   0, -1,  1, 0, 0,0],
["x",     -2,   1,  1,  0, 0, 0,0],
["y",      0,   0,  0,  0, 0, 0,0],
#["z",     -3,   0,  0, -1, 0, 0,0],  #z ions do not occur in Expert systems outline
["im",    -2,  -1,  0, -1, 0, 0,0]
],
columns=["loss"]+elements).set_index("loss")
#should I write these as dataframes or as dictionaries?
neutral_loss_formulas=neutral_loss_formulas.T.to_dict(orient='list')

#loss dict for side chain neutral losses
loss_dict={
"Sph":["H3PO4","HPO3"],
"Tph":["H3PO4","HPO3"],
"Yph":["H3PO4","HPO3"],
"C":["CH2S"],
"D":["H2O","CO2","C2H4O2"],
"E":["H2O","C2H4O2"],
"I":["C2H4"],
"K":["C2H5N","C4H9N","C4H11N","C3H9N"],
"L":["C3H6","C4H8"],
"M":["C2H4S","C3H6S"],
"Mox":["CH4SO","C3H8SO","C3H6SO"],
"N":["NH3","CH3NO","C2H5NO"],
"Q":["NH3", "CH3NO", "C2H5NO", "C3H5NO"],
"R":["NH3", "CH2N2", "C3H9N3"],
"S":["H2O","CH4O"],
"T":["H2O","C2H4O"],
"V":["C3H6"],
"W":["C8H7N","C9H9N"]}
# for k,v in loss_dict.items():
#     loss_dict.update({k:v+["a","b","c","y","x"]})


water=[2,0,1,0,0,0,0]


fragment_list=[]



#%%

#precompute possible combinations of amino acids, ions and neutral loss

compositions=[]
for k,v in aa_molecular_formulas.items():
    l=[[str([k]+v)]]
    
    #add side chain
    side_chain=loss_dict.get(k)
    if side_chain==None:
        side_chain=[""]
    else:
        side_chain=[""]+side_chain
    l+=[[str([sc]+neutral_loss_formulas.get(sc)) for sc in side_chain]]  

    #add water-loss
    l+=[[str(["",       0,   0,  0,  0, 0, 0,0]),
        str(["H2O",   -2,   0, -1,  0, 0, 0,0])]]
    #add ammonium-loss
    l+=[[str(["",       0,   0,  0,  0, 0, 0,0]),
        str(["NH3",   -3,  -0, -0, -1,-0,-0,0])]]
    #add CO-loss
    l+=[[str(["",       0,   0,  0,  0, 0, 0,0]),
        str(["CO",    -0,  -1, -1, -0,-0,-0,0])]]
    #add ions
    l+=[[str([ion]+neutral_loss_formulas.get(ion)) for ion in ["","a","b","c","y","x"]]]
    
    

    #compute all combiniations
    combinations=list(itertools.product(*l))
    carr=[np.array([eval(c) for c in cc]) for cc in combinations]
    for c in carr:
        title=[i for i in c[1:,0] if i !=""]
        title.sort()
        title="_".join([c[0,0]]+title)
        compositions.append([title]+list(c[:,1:].astype(int).sum(axis=0)))  
    
    # if k=="E":
    #     break
 
#eval back as lists
#compute chemical compositions, filter impossible compositions
df=pd.DataFrame(compositions,columns=["ion"]+elements).set_index("ion")
predf=df[~(df<0).any(axis=1)].reset_index().drop_duplicates().set_index("ion")


#%%

charge=2
#peptide='GMQCFSVTEIWR' #example
peptide="SLENETphLNK" #key: all modifications cam,ox,ac,ph are written lowercase

#%%

Orbitrap=True
fixed_Ccam=True


if type(peptide)!=type([]): 
    uppers=[ix for ix,i in enumerate(peptide) if i.isupper()]+[len(peptide)]
    peptide=[peptide[uppers[ix]:uppers[ix+1]]for ix,i in enumerate(uppers[:-1])]
    #split at uppercase
    
if fixed_Ccam:
    for aai,aa in enumerate(peptide):
        if aa=="C":
            peptide[aai]="Ccam"
        
#compute all combinations of side chain neutral losses of parent
nested_peptide=[[i] for i in peptide]
for ix,aa in enumerate(peptide):
    ld=loss_dict.get(aa)
    if ld==ld:
        res = [aa+"("+nl+")" for nl in ld]
        nested_peptide[ix].extend(res)
combinations=list(itertools.product(*nested_peptide))
parents=[i for i in combinations if "".join(i).count("(")<=3] #filter for less than 4 neutral losses

# for each combination add regular ions a,b,y (x if phosphorylated)
abc=[[[str(i)]+list(c[0:i]) for i in range(1,len(peptide))] for c in parents]
xyz=[[[str(i)]+list(c[-i:  ]) for i in range(1,len(peptide))] for c in parents]

# internal ions (aby)
res = list(set([tuple([str(i)+":"+str(j-1)]+list(parent[i: j])) for parent in parents for i in range(len(parent))
          for j in range(i + 1, len(parent) + 1)])) 
res=[list(i) for i in res]


#special neutral losses
sN=[]
# 90 Neutral loss at C-Term-H2O if rule49 is true and (annotation is C-terminal or parent ion) 
#and there is no modification at the C-term and neutral loss count equals 0 then add annotation with neutral loss of H2O from C-Term
if len(peptide[-1])==1: #unmodified cterm, and no neutral losses
    #internal xyz Cterm loss
    iC=[i for i in res if ":"+str(len(peptide)-1) in i[0]]
    iC=[["int"+i[0]]+list(i[1:]) for i in iC if "(" not in "".join(i)]
    for ix,i in enumerate(iC):
        iC[ix][-1]=iC[ix][-1]+"(H2O)"
    #parent Cterm loss
    p=["p"]+list(parents[0])
    p[-1]+="(H2O)"
    iC.append(p)
    #regular xyz Cterm loss
    xyzC=[["xyz"+i[0]]+i[1:] for i in xyz[0]]
    for ix,i in enumerate(xyzC):
        xyzC[ix][-1]=xyzC[ix][-1]+"(H2O)"
    iC.extend(xyzC)
    sN.extend(iC)

# 91 Neutral loss at N-Term-NH3 if rule49 is true and fragment sequence starts not with P 
#and annotation is N-terminal and there is no modification at the N-term and neutral loss count equals 0 then add annotation with neutral loss of NH3 from N-Term
if len(peptide[0])==1: #unmodified Nterm, and no neutral losses
    #internal abc Nterm loss
    iN=[i for i in res if "0:" in i[0]]
    iN=[["int"+i[0]]+list(i[1:]) for i in iN if "(" not in "".join(i)] #no modifications
    for ix,i in enumerate(iN):
        iN[ix][1]=iN[ix][1]+"(NH3)"
    #parent Nterm loss
    p=["p"]+list(parents[0])
    p[1]=p[1]+"(NH3)"
    iN.append(p)
    #regular abc Nterm loss
    abcN=[["abc"+i[0]]+i[1:] for i in abc[0]]
    for ix,i in enumerate(abcN):
        abcN[ix][1]=abcN[ix][1]+"(NH3)"
    iN.extend(abcN)
    sN.extend(iN)

# 92 Neutral loss at Parent-CO if rule49 is true and annotation is parent ion then add annotation with neutral loss of CO from parent ion
sN.extend([["p(CO)"] +list(i) for i in parents])

# 93 Neutral loss at internal fragment-CO if rule49 is true and annotation is internal fragment then add annotation with neutral loss of CO from internal fragment
resCO=[list(i) for i in res]
for ix,i in enumerate(resCO):
    resCO[ix][0]="int"+resCO[ix][0]+"(CO)"
sN.extend(resCO)

sN=[i for i in sN if "".join(i).count("(")<=3] #filter for less than 4 neutral losses

#extend regular ions
abc=[list(i) for i in set([tuple(i) for i in sum(abc,[])])]
xyz=[list(i) for i in set([tuple(i) for i in sum(xyz,[])])]

#abc
abc=[[i[0]+x]    +i[1:-1]+[i[-1]+x] for i in abc                                  for x in ["(a)","(b)","(c)"]]
abc+=[[i[0][3]+x]+i[1:-1]+[i[-1]+x] for i in [i for i in sN if i[0][0:3]=="abc"]  for x in ["(a)","(b)","(c)"]] #add sN abc

#xyz and pad
ions=["(y)"]
if "ph" in "".join(peptide): ions+=["(x)"]
    
xyz=[[i[0]+x    ]+[""]*(len(peptide)-int(i[0]))   +[i[1]+x]+i[2:] for i in xyz for x in ions] 
xyz+=[[i[0][3]+x]+[""]*(len(peptide)-int(i[0][3]))+[i[1]+x]+i[2:] for i in [i for i in sN if i[0][0:3]=="xyz"]  for x in ions] 

#res aby and pad


         #ion_type  pad_left           add_ion  add_rest    #pad_right                                                                                  #not a1b1/y1
res_aby= [[i[0]+x] +[""]*int(i[0][0]) +i[1:-1] +[i[-1]+x] +[""]*(len(peptide)-int(i[0][2])-1) for i in res for x in ["(a)","(b)"]                      if i[0][0]!="0"] 
res_aby+=[[i[0]+x] +[""]*int(i[0][3]) +i[1:-1] +[i[-1]+x] +[""]*(len(peptide)-int(i[0][5])-1) for i in sN  for x in ["(a)","(b)"]  if i[0][0:3]=="int" if i[0][3]!="0"]
res_aby+=[[i[0]+x] +[""]*int(i[0][0]) +[i[1]+x]+i[2:]     +[""]*(len(peptide)-int(i[0][2])-1) for i in res for x in ["(y)"]                            if i[0][0]!=str(len(peptide)-1)] 
res_aby+=[[i[0]+x] +[""]*int(i[0][3]) +[i[1]+x]+i[2:]     +[""]*(len(peptide)-int(i[0][5])-1)               for i in sN  for x in ["(y)"]        if i[0][0:3]=="int" if i[0][3]!=str(len(peptide)-1)]

aa_labels=["aa"+str(i) for i in range(len(peptide))]
parents=[["p"]+list(i) for i in parents]
all_fragments=pd.DataFrame(parents+res_aby+abc+xyz,columns=["ion_type"]+aa_labels).fillna("")

#%% aggregate neutral losses

uni_aa=np.unique(all_fragments[aa_labels].values)[1:]
titles=[]
for i in uni_aa:
    title=[i.replace(")","") for i in i.split("(")]
    title[1:].sort()
    titles.append([i,"_".join(title)])
titles=[i for i in titles if i[1] in predf.index]
m=predf.loc[[i[1] for i in titles]]
m.index=[i[0] for i in titles]
all_fragments=all_fragments[~all_fragments[aa_labels].isin(set(uni_aa[1:])-set(m.index)).any(axis=1)].reset_index(drop=True)

z=pd.DataFrame([0]*len(elements)).T
z.columns=m.columns
z.index=[""]
m=pd.concat([z,m],axis=0)

fd=pd.concat([m.loc[all_fragments[i]].reset_index(drop=True) for i in aa_labels])
s=fd.groupby(fd.index).sum().values
all_fragments[elements]=s
all_fragments["mass"]=np.multiply(s,ele_mass).sum(axis=1)


#tomorrow: test with some examples
#fast ppm difference search, put that stuff on github!
#x=np.multiply(aa_molecular_formulas.loc[list(peptide)][::-1]+water+ion_molecular_formulas.loc["x",:],ele_mass).cumsum()


# mas=[]
# for u in uni_aa:
#     ms=[i.replace(")","") for i in u.split("(")]
#     m=pd.concat([aa_molecular_formulas.loc[ms[0:1]],neutral_loss_formulas.loc[ms[1:]]]).sum()
#     m.name=u
#     mas.append(m)
# mas=pd.concat(mas,axis=1).T
#this can be precomputed.

#create index and neutral loss
#retrieve
#sum per index

#%% compute chemical compositions

#compute chemical compositions of each unique amino acid/neutral loss combination
#remove rows that contain negative compositions (give them -10 selenium)

#for each index, retrieve chemical compositions for column in columns
#sum by index

#remove CO if CO index

#%% add immonium



#immonium
# fragment_list=[]
# 2 Ready for immonium ions if mass analyzer was Orbitrap then immonium ions are possible else immonium ions are not possible
# rules 4-28 immonium AA if rule2 is true and peptide sequence contains AA and m/z value in spectrum then add immonium ion of AA to candidate list

 #immonium
# if Orbitrap:

#     im=aa_molecular_formulas.loc[set(peptide)]+ion_molecular_formulas.loc["im",:]
#     im.index="Im_"+im.index
#     fragment_list.append(im)

# 89 Neutral loss at IM K-NH3 if annotation is immonium ion and fragment sequence contains K and neutral loss count is equals 0 then add annotation with neutral loss of NH3 from immonium ion K


# side chain

# 3 Ready for side chain ions if mass analyzer was Orbitrap then side chain loss fragments are possible else side chain loss fragments are not possible

# 29 side chain by R if rule3 is true and peptide sequence contains R and m/z value in spectrum then add side chain loss of C4H6N2 to candidate list
# 30 side chain by H if rule3 is true and peptide sequence contains H and m/z value in spectrum then add side chain loss of C4H6N2 to candidate list
# 31 side chain by K if rule3 is true and peptide sequence contains K and m/z value in spectrum then add side chain loss of C5H11N to candidate list
# 32 side chain by W if rule3 is true and peptide sequence contains T and m/z value in spectrum then add side chain losses of C10H9N, C9H9N and C9H7N to candidate list

#These side chain fragments are highly ambiguous between the two papers with:
    #Michalski, Annette, et al. "A systematic investigation into the nature of tryptic HCD spectra." Journal of proteome research 11.11 (2012): 5479-5491.
    # here they are descibed as c5h9n3+h+ for r, c5h9n+h+ for k and c9h9n+h+, c9h7n+h+ for W

# from pyteomics import mass
# if Orbitrap:
#     for aa in peptide:
#         if aa=="R":
#             fragment_list.append(["sidechain_C4H6N2_R",mass.calculate_mass(formula='C5H9N3')+water_mass])
#         if aa=="H":
#             fragment_list.append(["sidechain_C4H6N2_H",mass.calculate_mass(formula='C5H9N3')+water_mass])

#         if aa=="K":
            
#         if aa=="W":

#%% add side chain losses

#%% remove impossible compositions



# # remove unstable regular ions (c1,b1,a1,a2)

# #41 unstable a1 if annotation is a1 and has no modification a1 is unstable go to rule40
# abc=[i for i in abc if i[0][2]!="a" or i[0][0]!="1" or len(i[1].split("(")[0])>1]
# # 42 unstable a2 if annotation is a2 and contains Q or N or has a Carbamidometylated C a2 is unstable go to rule40
# abc=[i for i in abc if i[0][2]!="a" or i[0][0]!="2" if sum([x in "".join(i) for x in ["Q","N"]])]
# # 43 unstable b1 if annotation is b1 and peptide sequence is not starting with R, H, K or Acetylated M,S, A b1 is unstable go to rule40



#%% compute masses

#%% compute mzs

#%% add scores

#%% groupby mass, pick best score

#%%


# #standard ions
# # 33 Parent ion if fragmentation is HCD or CID add parent ion to queue
# parent=pd.DataFrame((aa_molecular_formulas.loc[list(peptide)]+water).sum()).T
# parent.index=["parent_"+"".join(peptide)]
# fragment_list.append(parent)


# # 34 a-ion series if fragmentation is HCD or CID add a ion series to queue
# a=(aa_molecular_formulas.loc[list(peptide)]+water+ion_molecular_formulas.loc["a",:]).cumsum()
# a.index=["a"+str(ix+1)+"_"+"".join(a.index[0:ix+1]) for ix,_ in enumerate(a.index)]
# # 41 unstable a1 if annotation is a1 and has no modification a1 is unstable go to rule40
# if len(peptide[0])>1: #len>1 means modified
#     fragment_list.append(pd.DataFrame(a.iloc[0]).T)
# # 42 unstable a2 if annotation is a2 and contains Q or N or has a Carbamidometylated C a2 is unstable go to rule40
# if sum([i in peptide[1] for i in ["Q","N","Ccam"]]):
#     fragment_list.append(pd.DataFrame(a.iloc[1]).T)
# fragment_list.append(a.iloc[2:])

# # 35 b-ion series if fragmentation is HCD or CID add b ion series to queue
# b=(aa_molecular_formulas.loc[list(peptide)]+water+ion_molecular_formulas.loc["b",:]).cumsum()
# b.index=["b"+str(ix+1)+"_"+"".join(b.index[0:ix+1]) for ix,_ in enumerate(b.index)]
# # 43 unstable b1 if annotation is b1 and peptide sequence is not starting with R, H, K or Acetylated M,S, A b1 is unstable go to rule40
# if peptide[1] not in ["R", "H", "K","Mac","Sac","Aac"]:
#     fragment_list.append(pd.DataFrame(b.iloc[0]).T)
# fragment_list.append(b.iloc[1:])

# # 36 y-ion series if fragmentation is HCD or CID add y ion series to queue
# y=(aa_molecular_formulas.loc[list(peptide)][::-1]+water+ion_molecular_formulas.loc["y",:]).cumsum()
# y.index=["y"+str(ix+1)+"_"+"".join(y.index[0:ix+1]) for ix,_ in enumerate(y.index)]
# fragment_list.append(y)

# # 37 stable c1 if second amino acid is a Gln or Asn add c1 ion to queue 
# if sum([i in peptide[1] for i in ["Q","N"]]): #not sure if contains acetylation, wording of 37 is diff3erent from 42
#     c=(aa_molecular_formulas.loc[list(peptide)]+water+ion_molecular_formulas.loc["c",:]).cumsum()
#     c.index=["c"+str(ix+1)+"_"+"".join(c.index[0:ix+1]) for ix,_ in enumerate(c.index)]
#     fragment_list.append(c.iloc[0]) #c1 only

# # 38 stable x by Phosphorylation if peptide is phosphorylated add x ions to queue
# if "phe" in "".join(peptide):
#     x=(aa_molecular_formulas.loc[list(peptide)][::-1]+water+ion_molecular_formulas.loc["x",:]).cumsum()
#     x.index=["x"+str(ix+1)+"_"+"".join(x.index[0:ix+1]) for ix,_ in enumerate(x.index)]
#     fragment_list.append(x)

# # Neutral
# fragment_list=pd.concat(fragment_list)
# fragment_list["neutral_loss"]=""

# # 49 Neutral loss at phS if fragment sequence contains S and has a Phosphorylation on S and neutral loss count is less than 3 then add annotation with neutral losses of H3PO4 and HPO3 from Phosphorylated S to the queue
# # 50 Neutral loss at phT if fragment sequence contains T and has a Phosphorylation on T and neutral loss count is less than 3 then add annotation with neutral losses of H3PO4 and HPO3 from Phosphorylated T to the queue
# # 51 Neutral loss at phY if fragment sequence contains Y and has a Phosphorylation on Y and neutral loss count is less than 3 then add annotation with neutral losses of H3PO4 and HPO3 from Phosphorylated Y to the queue
# has_p=fragment_list[fragment_list["Phosphor"]>0]
# has_p[elements]-neutral_loss_formulas.loc["H3PO4"]
# has_p["neutral_loss"]=has_p["neutral_loss"]+"-"+"H3PO4"
# fragment_list=pd.concat([fragment_list,has_p])

# # 52 Ready for Neutral loss if annotation is in spectrum and neutral loss count is less than 3 then annotation can have neutral losses else annotation cannot have neutral losses

# # 53 Ready for Internal Fragmentation if annotation is from a-,b- or y-series then annotation can have internal fragments else annotation cannot have internal fragments
# res = [[str(i)+":"+str(j),"".join(peptide[i: j])] for i in range(len(peptide))
#           for j in range(i + 1, len(peptide) + 1)]

# internal_ions=[]
# for r in res:
#     internal_ions.extend([["a"+r[0]+"_"+r[1]]+(aa_molecular_formulas.loc[list(r[1])]+water+ion_molecular_formulas.loc["a",:]).cumsum().values.tolist()[0],
#                           ["b"+r[0]+"_"+r[1]]+(aa_molecular_formulas.loc[list(r[1])]+water+ion_molecular_formulas.loc["b",:]).cumsum().values.tolist()[0],
#                           ["y"+r[0]+"_"+r[1]]+(aa_molecular_formulas.loc[list(r[1])]+water+ion_molecular_formulas.loc["y",:]).cumsum().values.tolist()[0]])

# internal_ions=pd.DataFrame(internal_ions,columns=["ion"]+elements).set_index("ion")
# internal_ions["neutral_loss"]=""
# fragment_list=pd.concat([fragment_list,internal_ions])




#remove any impossible compositions (-1)
#calculate mass
#x=np.multiply(aa_molecular_formulas.loc[list(peptide)][::-1]+water+ion_molecular_formulas.loc["x",:],ele_mass).cumsum()

#%% charge
# 44 Charge 1+ if annotation charge is 1+ and peptide charge is greater than or equals 1+ and m/z value in spectrum add annotation to candidate list
# 45 Charge 2+ if annotation charge is 2+ and peptide charge is greater than or equals 2+ and m/z value in spectrum increase charge to 2+ and add annotation to candidate list
# 46 Charge 3+ if annotation charge is 3+ and peptide charge is greater than or equals 3+ and m/z value in spectrum increase charge to 3+ and add annotation to candidate list
# 47 Charge 4+ if annotation charge is 4+ and peptide charge is greater than or equals 4+ and m/z value in spectrum increase charge to 4+ and add annotation to candidate list
# 48 Charge 5+ if annotation charge is 5+ and peptide charge is greater than or equals 5+ and m/z value in spectrum increase charge to 5+ and add annotation to candidate list


#%% 

#add scores

# 94 General Enumeration end if there are more elements in the queue then go to rule39
# 95 Priority Enumeration start if there is an element in the candidate list get next elements in the candidate list and initialize score of the annotation to 0
# 96 Filter for null values if annotation is null discard annotation and remove from candidate list and go to rule95
# 97 Priority B Rule if annotation ion type is b then increase score by 100
# 98 Priority Y Rule if annotation ion type is y then increase score by 100
# 99 Priority Parent Rule if annotation ion type is parent ion then increase score by 99
# 100 Priority A Rule if annotation ion type is a then increase score by 98
# 101 Priority C Rule if annotation ion type is c then increase score by 98
# 102 Priority IM Rule if annotation is immonium ion then increase score by 97
# 103 Priority Internal with P Rule if annotation is internal fragment and fragment sequence starts with P then increase score by 80
# 104 Priority Internal without P Rule if annotation is internal fragment and fragment sequence starts without P then increase score by 70
# 105 Priority Neutral loss Rule if annotation has neutral losses then decrease score by 5 times number of H2O, NH3 and CO and 30 times number of all other losses
# 106 Priority Enumeration end if there are more elements in the candidate list then add scored annotation to result list and go to rule39


#%% Priority MAE

# Only one fragment ion classification is made for each Average Mass Ion even when several are possible
# (for example, analysis of the second, fourth, and fifth peptide bond in DTA files from Sample 5 indicated that 8% of doubly
#  charged ions are isomeric with a singly charged ion). Rules commonly utilized in manual analysis are used to assign the most likely annotation, 
# but all other possibilities are listed in the summary output. Fragment ions are classified in the following order: 
#     1) parent ion and ions generated by neutral losses of water/ammonia, the guanidinium side chain (only allowed when Arg is present), or H2CO3; 
#     2) singly charged “canonical” bn and yn ions produced by cleavage at one peptide bond; 
#     3) singly charged canonical an ions; 
#     4) singly charged dehydrated/deammoniated canonical an, bn, and yn ions and C-terminal rearrangements (bn− 1 + 18); 
#     5) singly charged multiple dehydrated/deammoniated canonical bn ions; (6) 
#         doubly charged canonical bn and yn ions; (7) 
#         doubly charged canonical an ions; (8) 
#         doubly charged dehydrated/deammoniated canonical an, bn, and yn ions and C-terminal rearrangements (bn− 1 + 18); (9) 
#         doubly charged multiple dehydrated/deammoniated canonical bn ions; (10) 
#         triply charged canonical bn and yn ions; (11) 
#         triply charged canonical an ions; (12) 
#         triply charged dehydrated/deammoniated canonical an, bn, and yn ions and C-terminal rearrangements (bn− 1 + 18); (13)
#         triply charged multiple dehydrated/deammoniated canonical bn ions; (14) 
#         b ions generated by internal fragmentation; (15) 
#         a ions from internal fragmentation; and (16) 
#         dehydrated/deammoniated a and b ions from internal fragmentation.

# Rules for Fragment Ions Generated by Secondary Cleavages—
# To decide between alternative ion types and to minimize chance assignments due to the large number of combinatorial possibilities, heuristic rules for fragment ion annotations based on simple chemical rules for multistep cleavages (indicated by →) or parallel reactions that are expected to be independent of each other (indicated by ‖) are applied as follows.
# 1.
# Do not assign a dehydrated (▵) or deammoniated ion (▿) if the corresponding unmodified form is absent 
# (using the rule that “unmodified → ▵/▿”). We use the standard ▵ symbol for dehydration, but in some cases, the ion shows mass
# closer to a deammoniated ion. Stochastic variation in intensity of different isotope peaks produces ambiguity between dehydrated and deammoniated ions
# (see Fig. 2, expanded view panels that show individual ions, keeping in mind that dehydration and deammoniated ions are 1 Da apart and have overlapping 
# isotope peaks). Therefore, we utilize a ▿ symbol for potentially deammoniated ions, instead of the more common −NH3 nomenclature, to emphasize that 
# the evidence for the observed deammoniated ions is inadequate to distinguish from dehydration.
# 2.
# Do not assign an a ion if the corresponding b ion is absent except for a ions generated by internal fragmentation. 
# This rule can be expressed as: “b → a.”
# 3.
# Do not assign a dehydrated/deammoniated (▵/▿) a ion if the ratio of intensities between the ▵/▿ a ion and 
# its corresponding a ion (Ra) is significantly different from the ratio of intensities between the ▵/▿ b ion and 
# its corresponding b ion (Rb) generated by cleavage of the same peptide bond. 
# This assumes that the two reactions are independent: “(a → ▵/▿ a) ‖ (b → ▵/▿ b)”.
# This rule is utilized when |Ra − Rb| < 0.15(Ra + Rb) to avoid misclassifying cases where small changes could be due to stochastic differences in ion counting.
# 4.
# Do not assign a product of C-terminal cleavage to Pro unless no other alternative is possible. 
#This rule is based on the known chemistry of Pro (12, 19) and forces the program to change the assignment order for fragment ions when this situation occurs.
# 5.
# Do not assign internal fragment ions unless there is a moderately high likelihood for chemical cleavage at both ends. 
# This rule assumes that the two cleavages are independent: “(cleavage at site 1) ‖ (cleavage at site 2)”
#  where site 1 and site 2 represent the peptide bonds defining the ends of the internal fragment ion.
#  We assume that a fragmentation event has similar chemistry whether it represents the first or second cleavage. 
#  Thus, we estimate the likelihood of cleavage at each peptide bond from the intensities of the canonical a, b, or y ion generated 
#  by fragmentation at each site that are directly related to the reactivity of that site. 
#  The probability of obtaining an internal fragmentation between the sth and tth peptide bonds is estimated by Equation 6,
# Probabilityst= ∑ Is∑ I×∑ It∑ I 
# Eq. 6
# where Is and It are the intensities of each observed ion identified as sequence-specific a, b, or y 
#ions representing cleavage at the s and t peptide bonds and I is the intensity of each observed ion of any type except the internal fragment ions.
# 6.
# For internal fragment ions, the order of assignment prioritizes b over a fragment ions 
#and unmodified fragment ions over ▵/▿-modified a or b ions (for example Δan). This rule assumes that the loss of H2O/NH3 and CO are independent: “(b → a) ‖ (unmodified b → ▵/▿b).”
# 7.
# When two or three ▵/▿ derivatives of a b ion are present (e.g. bn, ▵bn, ▵▵bn, and ▵▵▵bn), 
#the related ions must follow intensity patterns that assume sequential reactions: “unmodified → ▵ → ▵▵ → ▵▵▵.”
#For example, a set of three related ions should show intensity patterns bn ≥ ▵bn ≥ ▵▵bn, bn ≤ ▵bn ≤ ▵▵bn or bn ≤ ▵bn ≥ ▵▵bn and exclude the pattern bn ≥ ▵bn ≤ ▵▵bn. If the set fails this test, the ▵▵ form should not be assigned. Similar patterns should be assessed for sets of four ions; if the set fails the test, the ▵▵▵ ion should not be assigned, and the first three ions should be reconsidered. This annotation test considers alternative assignments and low intensity ions in the sDTA files when testing for the unmodified element in this series and requires that the peptide sequence has sufficient Ser/Thr to account for the multiple dehydration events (number of Ser/Thr ≥ number of ▵ allowed).