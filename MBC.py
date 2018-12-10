
# coding: utf-8

# # MBC for Molecular Biology Companion

# In[1]:

def Rev_com(string):
    DNA = list(string)
    #print (DNA)
    for j in reversed(DNA):
        #print (j)
        if j == 'a':
            Rev_Com.append('t')
        elif j == 't':
            Rev_Com.append('a')
        elif j == 'g':
            Rev_Com.append('c')
        elif j == 'c':
            Rev_Com.append('g')
        elif j == 'A':
            Rev_Com.append('T')
        elif j == 'T':
            Rev_Com.append('A')
        elif j == 'G':
            Rev_Com.append('C')
        elif j == 'C':
            Rev_Com.append('G')
        else :
            Rev_Com.append(j)
    Rev_Com_str = ''.join(Rev_Com)
    #print (Rev_Com_str)
    return Rev_Com_str


# In[2]:

def exact_search(reference_string,target_string):
    reference = []
    target = []
    for i in reference_string:
        reference.append(i)
    for j in target_string:
        target.append(j)
    length = len(target)
    list_pos_match = []

    for i in range(0, len(reference)) :
        if len(reference) >= i + length:
            list_pos_match.append(reference[i : i + length])

    if target in list_pos_match:
        print ("it is a match")
        return True
    else:
        print ('it is not a match')
        return False


# In[3]:

def bidirec_exact_search(reference_string,target_string):
    rev_com_str = Rev_com(reference_string)
    ans1 = exact_search(reference_string,target_string)
    ans2 = exact_search(rev_com_str,target_string)
    
    if ans1 == True:
        print('the target matches in the same direction')
    if ans2 == True:
        print('the target matches in reverse complement direction')
    else:
        print('no match in either direction')
    return


# In[4]:

def translation(DNA_string):
    codon_table = {'atg':'Met' ,'tgg':'Trp','gct':'Ala', 'gcc':'Ala', 'gca':'Ala', 'gcg':'Ala', 'cgt':'Arg', 'cgc':'Arg', 'cga':'Arg','cgg':'Arg','aga':'Arg','agg':'Arg',
               'aat':'Asn','aac':'Asn','gat':'Asp' , 'gac':'Asp' , 'tgt':'Cys' ,'tgc':'Cys', 'caa':'Gln' , 'cag':'Gln' , 'gaa':'Glu' , 'gag':'Glu', 'ggt':'Gly', 'ggc':'Gly' , 'gga':'Gly' , 'ggg':'Gly',
               'cat':'His' , 'cac':'His'  , 'att':'Ile'  , 'atc':'Ile' , 'ata':'Ile' , 'tta':'Leu' , 'ttg':'Leu' , 'ctt':'Leu' , 'ctc':'Leu' , 'cta':'Leu' , 'ctg':'Leu'  , 'aaa':'Lys' , 'aag':'Lys' , 'ttt':'Phe', 'ttc':'Phe',
               'cct':'Pro' , 'ccc':'Pro' , 'cca':'Pro' , 'ccg':'Pro' , 'tct':'Ser' , 'tcc':'Ser' , 'tca':'Ser' , 'tcg':'Ser' , 'agt':'Ser' , 'agc':'Ser' ,  'act':'Thr' , 'acc':'Thr' , 'aca':'Thr' , 'acg':'Thr' ,
               'tat':'Tyr', 'tac':'Tyr' , 'gtt':'Val' , 'gtc':'Val' , 'gta':'Val' , 'gtg':'Val' }
    length = len(DNA_string)
    if length < 3:
        return
    temp_list = []
    temp_list_del0 = []
    temp_list_del1 = []
    DNA_list = []
    for i in DNA_string:
        DNA_list.append(i)
    for i in range(0, len(DNA_list)):
        temp_list.append(DNA_list[i])
        if len(temp_list_del0) + 1 == len(DNA_list):
            pass
        else:
            temp_list_del0.append(DNA_list[i+1])
            
        if len(temp_list_del1) + 2 == len(DNA_list):
            pass
        else:
            temp_list_del1.append(DNA_list[i+2])


    def translate(list,codon_table):
        codon_list = []
        for i in range(0, len(list)):
            if len(list) >= 3:
                a = list[0]
                b = list[1]
                c = list[2] 
                d = a + b +c
                del(list[0:3])
                codon_list.append(d)
            else:
                break
        protein_list = []
        while len(codon_list) > 0 :
            i = codon_list[0]
            if i != 'atg' and len(protein_list) == 0 :
                codon_list.pop(0)
            elif i == 'taa' or i == 'tag' or i == 'tga':
                break
            else:
                protein_list.append(codon_table[i])
                codon_list.pop(0)
        
        protein_string = '-'.join(protein_list)
        
        return protein_string
    ans1 = translate(frame1_list,codon_table)
    ans2 = translate(frame2_list,codon_table)
    ans3 = translate(frame3_list,codon_table)
    #print("frame 1 reads " , ans1)
    #print("frame 2 reads " , ans2)
    #print("frame 3 reads " , ans3)
    if len(ans1) > len(ans2) and len(ans1) > len(ans3):
        return ans1
    elif len(ans2) > len(ans1) and len(ans2) > len(ans3):
        return ans2
    
    else:
        return ans3


# In[5]:

def letter3to1(threeletteraminoacid_string, sep) :
    amino_3to1letter = { 'Arg':'R' , 'His':'H' , 'Lys':'K' , 'Asp':'D' , 'Glu':'E' , 'Ser':'S' , 'Thr': 'T' , 'Asn':'N' , 'Gln':'Q' , 'Cys':'C' , 'Gly':'G' , 'Pro':'P' ,
                    'Ala':'A' , 'Val' : 'V' , 'Ile' : 'I' , 'Leu':'L' , 'Met':'M' , 'Phe':'F' , 'Tyr':'Y' , 'Trp':'W' }
    temp_list = []
    
    for i in threeletteraminoacid_string:
        temp_list.append(i)
    while '-' in temp_list:
        temp_list.remove('-')
    temp2_list = []   
    while len(temp_list) > 0 :
        a = temp_list[0]
        b = temp_list[1]
        c = temp_list[2]
        d = a + b + c
        temp2_list.append(d)
        del(temp_list[0:3])
        
    oneletter_list = []
    for i in temp2_list :
        x = amino_3to1letter[i]
        oneletter_list.append(x)
    oneletter_string = sep.join(oneletter_list)
    
    return oneletter_string
    
    


# In[6]:

def mol_weight(amino_string):
    weight_dict = { 'A':89, 'R':174,  'N':132 , 'D':133 , 'C':121 , 'Q':146 , 'E':147 , 'G':75 , 'H':155, 'I':131 , 'L':131 , 'K':146 ,'M' :149 ,
               'F':165 , 'P':115 , 'S':105 , 'T':119 , 'W' :204 , 'Y' :181 , 'V' :117}
    temp_list = []
    if amino_string.find('a') ==1 or amino_string.find('e') ==1 or amino_string.find('l') ==1 or amino_string.find('n') ==1 or amino_string.find('r') ==1 or amino_string.find('s') ==1 or amino_string.find('y') ==1   : 
        one_letter_string = letter3to1(amino_string,'')
    else:
        for i in amino_string:
            temp_list.append(i)
        while '-' in temp_list:
            temp_list.remove('-')
        while ' ' in temp_list:
            temp_list.remove(' ')
        one_letter_string = ''.join(temp_list)
    one_letter_list = []
    for i in one_letter_string:
        one_letter_list.append(i)
    mol_weight = 0
    for i in one_letter_list:
        mol_weight = mol_weight + weight_dict[i]
    
    return mol_weight

def aa_composition(amino_string, optional_argument = 'percent'):
	'''only take in one-letter coded amino acid strings.'''
	length = len(amino_string)
	temp_dict = {}
	for i in amino_string:
		if i not in temp_dict.keys():
			temp_dict[i] = 1
		if i in temp_dict.keys():
			temp_dict[i] += 1

    if optional_argument == 'percent':
        percent_dict = temp_dict
        for i in percent_dict:
		  percent_dict[i] = (percent_dict[i]/length)*100
        return percent_dict

    if optional_argument == 'count':
        return temp_dict

    if optional_argument == 'fraction':
        fraction_dict = temp_dict
        for i in percent_dict:
          fraction_dict[i] = (fraction_dict[i]/length)
        return fraction_dict










