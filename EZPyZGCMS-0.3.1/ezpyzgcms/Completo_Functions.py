import numpy as np
import itertools
import math
from collections import Counter
import copy
from functools import reduce
import re


'Constants'
#Different Isotopes and their probabilities
H1=0.999885
H2=.000115
N14=.99632
N15=.00368
Si28=.922297
Si29=.046832
Si30=.030872
S32=.9493
S33=.0076
S34=.0429
S36=.0002
O16=.99758
O17=.00038
O18=.00204
C12=.989184
C13=.010816
Cl35=0.7578
Cl36=0
Cl37=0.2422
P31=1.0


def comboFinder(atom, numberofatoms):
    "This function takes an atom type (ie O) and a molecule size (ie 6) and returns every possible combination with the atom's different isotypes as well as their probabilities "
    isotopes = {'C':[12,13],'H':[1,2],'O':[16,17,18], 'S':[32,33,34,36],'Si':[28,29,30], 'N':[14,15], 'Cl':[35,36,37], 'P':[31]}

    isoProbs = {12: C12, 13: C13, 1: H1, 2: H2,16: O16, 17: O17, 18: O18, 14: N14, 15: N15, 28: Si28,
                29: Si29, 30: Si30, 32: S32, 33: S33, 34: S34, 36:S36, 35:Cl35, 37:Cl37, 31:P31}

    ions = [isotopes[atom][0]*numberofatoms,isotopes[atom][-1]*numberofatoms]
    ions = list(range(ions[0],ions[1]+1))

    combinations = list(itertools.combinations_with_replacement(isotopes[atom], numberofatoms))

    if len(isotopes[atom])==3:
        combinations+list(itertools.combinations_with_replacement((isotopes[atom][0],isotopes[atom][1]), numberofatoms))
        combinations+list(itertools.combinations_with_replacement((isotopes[atom][1],isotopes[atom][2]), numberofatoms))
    final={x:[] for x in ions}
    for x in combinations:
        if int(sum(x)) in ions:
            final[int(sum(x))].append(x)
    finalProbs = {x:[[isoProbs[q] for q in y] for y in final[x]] for x in final}
    return finalProbs

def chemicalFormulator(chemicalFormula):
    "input is a string such as 'C1H1N1' denoting a chemical formula; this function takes in a chemical formula and returns a dictionary"
    #Change all letters to uppercase
    inputup=chemicalFormula.upper()
    fixed=[]
    #change I to lowercase (this accounts for Si)
    if any(x in 'I' for x in inputup):
        inputup=inputup.replace('I','i')

    #Start to iterate through the string and identify composition
    #Adds hyphens and values of 1 when necessary to reformat string
    for i in range(len(inputup)-1):

        if inputup[i].isupper():
            if inputup[i+1].isupper():
                fixed.append(inputup[i]+'-'+'1')
            if inputup[i+1].islower():
                    fixed.append(inputup[i])
            elif inputup[i+1].isnumeric():
                fixed.append(inputup[i]+'-')

        if inputup[i].isnumeric():

            if inputup[i+1].isupper():
                fixed.extend(inputup[i])
            else:
                fixed.extend(inputup[i])

        elif inputup[i].islower():
            if inputup[i+1].isupper():
                fixed.append(inputup[i]+'-')
                fixed.append('1')
            else:
                fixed.append(inputup[i]+'-')

    #This portion takes into account the last character separately
    #This avoids index out of range error
    x=len(inputup)

    if inputup[x-1].isnumeric():
        fixed.append(inputup[x-1])
    else:
        fixed.append(inputup[x-1]+'-'+'1')

    compiled=','.join(fixed)
###Converting to a more dictionary friendly version

    if (i in (',') for i in compiled):
            compiled=compiled.replace(',','')

    if (i in ('-') for i in compiled):
            compiled=compiled.replace('-',':')

    grouped=re.findall('([CHNOS][i]?:[0-9][0-9]?)',compiled)

    d={}
    for i in grouped:
        j=i.split(':')
        d[j[0]]=int(j[1])
    elements1=d
    return elements1



def Pvalues(combolist, numax, clabeled):
    """
    x,y,z are fractional abundances of isotopes possible in nature e.g. O16 input calls .99757 fractional abundance
    x1,y1,z1 refers to the number of that particular isotope
    n referes to the total number of that element regardless of isotope
    p is the probability of an elemental isotopomer of p1 meaning only one elemental isotopomer

    """
    CM = np.zeros((len(combolist),len(combolist)))
    combolistdetails = [[Counter(combolist[x][y]) for y in range(0,len(combolist[x]))] for x in combolist]
    lowerisodict = {0.989184:'C12', 0.999885:'H1',0.99632:'N14', 0.99758:'O16'}
    combolistdog = copy.deepcopy(combolistdetails)

    def copier(combolistykno):
            ykno = copy.deepcopy(combolistykno)
            for x in ykno:

                for y in x:

                    for q in y:

                        if q in lowerisodict:
                            y[q]=y[q]-1

            return ykno
    labelcmds = {}
    if clabeled==True:
        for x in range(1,len(combolistdog)):
            labelcmds[x] = copier(combolistdog)
            combolistdog = copier(copy.deepcopy(combolistdog))


        Otcms = {x:np.zeros((len(combolist),len(combolist))) for x in range(1,len(combolistdog))}

    def p3(x,x1,y,y1,z,z1, numax):
        f=math.factorial
        if x1+y1+z1==numax:
            try:
                pvalue3=(f(numax))*(((x**x1)/f(x1))*((y**y1)/f(y1))*((z**z1)/f(z1)))
            except ValueError:
                pvalue3 = 0
            finally:
                return (pvalue3)
        else:
            pvalue3 = 0
            return pvalue3

    def processor(combolistdetailykno, CM, numax, column):
        for r,q in enumerate(combolistdetailykno):
            for j in q:
                if len(j)==0:
                    print('empty')
                elif len(j)==1:
                    CM[r][column]=CM[r][column]+p3(list(j.keys())[0], list(j.values())[0], 0,0,0,0, numax)
                elif len(j)==2:
                    CM[r][column]=CM[r][column]+p3(list(j.keys())[0], list(j.values())[0], list(j.keys())[1], list(j.values())[1],0,0, numax)
                elif len(j)==3:
                    CM[r][column]=CM[r][column]+p3(list(j.keys())[0], list(j.values())[0], list(j.keys())[1], list(j.values())[1], list(j.keys())[2], list(j.values())[2], numax)
        return CM

    CM = processor(combolistdetails, CM, numax, 0)
    if clabeled == True:
        OtherCMs = {x:processor(labelcmds[x], Otcms[x], numax-x, x)[:-x] for x in labelcmds}
        for x in OtherCMs:
            for y,z in enumerate(OtherCMs[x]):
                CM[y+x] = CM[y+x]+z
    else:
        for x,y in enumerate(CM):
            np.fill_diagonal(CM[x:],CM[x][0])

    return CM



def completoCM(metabolite, listoflabeledatoms):
    """
    completoCM takes a metabolite of class Item() and a list of atoms that are labeled and produces an inverted correction matrix to correct the
    raw abundances for that metabolite.

    metabolite.chemicalFormula = {'C':3, 'H': 2, 'O':1}
    listoflabeledatoms = ['C','H','O'] if all labeled

    """

    metaboChem = chemicalFormulator(metabolite.chemicalFormula)
    combolists = {x:comboFinder(x,metaboChem[x]) for x in metaboChem}
    PvaluesDictUnlab = {x: Pvalues(combolists[x], metaboChem[x], clabeled = False) for x in combolists if x not in listoflabeledatoms}
    PvaluesDictlab = {x: Pvalues(combolists[x], metaboChem[x], clabeled = True) for x in combolists if x in listoflabeledatoms}
    final  = {}
    final.update(PvaluesDictUnlab)
    final.update(PvaluesDictlab)

    try:
        shaper  = (np.diff(metabolite.Ions)[0]+1,np.diff(metabolite.Ions)[0]+1)
    except:
        shaper = (1,1)

    emptyCMs1 = {x: np.lib.pad(final[x], ((0,shaper[0]-np.shape(final[x])[0]),(0,shaper[0]-np.shape(final[x])[0])), 'constant', constant_values=(0)) if np.shape(final[x])[0]<shaper[0] else final[x][:shaper[0]][:,0:shaper[0]] for x in final}
    listofCM=copy.deepcopy(emptyCMs1)

    for z in emptyCMs1:
        try:
            invCM = np.linalg.inv(emptyCMs1[z])
            listofCM[z]=emptyCMs1[z]
        except:
            nCM = np.zeros(shaper)
            for y,x in enumerate(emptyCMs1[z][:,0]):
                np.fill_diagonal(nCM[y:][:,0:shaper[0]-y],x)

            listofCM[z]=nCM

    CM = reduce(lambda x,y: np.dot(x,y), listofCM.values())

    invCM = np.linalg.inv(CM)

    return invCM


def completoData(invCM, datalistmetabo):

    """
    completoData takes in an inverse correction matrix produced by completoCM and raw data for the relevant metabolite. It then takes the dot product
    of the two matrices and returns the corrected abundances

    """
    finaldat = np.dot(invCM, datalistmetabo)
    finalenrich = [x/sum(finaldat) for x in finaldat]
    return finaldat, finalenrich


##Rudimentary functions for reference
def p1(x,x1,n):
    f=math.factorial
    try:
        pvalue1=(f(n))*((x**x1/f(x1)))
    except ValueError:
        pvalue1=0
    finally:
        return (pvalue1)

def p2(x,x1,y,y1,n):
    f=math.factorial

    pvalue2=(f(n))*((x**x1/f(x1))*(y**y1/f(y1)))
    return (pvalue2)

def p3OG(x,x1,y,y1,z,z1,n):
    f=math.factorial

    pvalue3=(f(n))*((x**x1/f(x1))*(y**y1/f(y1))*(z**z1/f(z1)))
    return (pvalue3)
