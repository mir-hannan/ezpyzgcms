from pyms.GCMS.IO.ANDI import ANDI_reader
from pyms.Noise.SavitzkyGolay import savitzky_golay
from pyms.IntensityMatrix import build_intensity_matrix_i
from pyms.TopHat import tophat
import pandas as pd
from Completo_Functions import completoCM, completoData
import openpyxl as xl
from openpyxl.utils.dataframe import dataframe_to_rows
from os import path
from glob import glob
import matplotlib.pyplot as plt


def find_ext(dr, ext):
    return glob(path.join(dr, "*.{}".format(ext)))


def grabber(df, ref, MetaboliteName): 
    # grabs certain metabolites(x) off of a reference
    q = df.iloc[df.index.get_loc(df.loc[MetaboliteName].name):
                df.index.get_loc(df.loc[MetaboliteName].name)+ref[MetaboliteName]]
    return q


class Item:
    """
    This class is the primary information holder for our analysis. It takes parameters from our library (i.e. RT, formula, Ions, and Name) and then
    systematically parses the CDF for abundance information on the relevant metabolite. It then corrects the isotope abundances for natural abundance 
    to get accurate tracing data.

    """
    instances = []

    def __init__(self, Name, RT, Ions, chemicalFormula, mainC, labeledspecies, im):
        self.Name = Name
        self.Derivatization = 'TBDMS'
        self.RT = RT
        self.RTis = RT * 60  # Retention time in seconds
        self.Ions = Ions  # Ions that you're interested in
        self.IonL = list(range(Ions[0], Ions[1]+1))

        for i in self.IonL:
            ic = im.get_ic_at_mass(i)
            ic1 = savitzky_golay(ic)
            ic_smooth = savitzky_golay(ic1)
            ic_base = tophat(ic_smooth, struct="1.5m")
            im.set_ic_at_index(im.get_index_of_mass(i), ic_base)


        self.ics = {x: im.get_ic_at_mass(x) for x in range(self.Ions[0], self.Ions[1] + 1)}
        self.RTIndex = self.ics[self.Ions[0]].get_index_at_time(self.RTis)
        self.RTadj = [self.ics[self.Ions[0]].get_index_at_time(self.RTis-15), self.ics[self.Ions[0]].get_index_at_time(self.RTis+15)]
        self.icsval = {x: im.get_ic_at_mass(x).intensity_array for x in range(self.Ions[0], self.Ions[1]+1)}

        self.IntensityDF = pd.DataFrame(self.icsval)
        self.IntensityDF['RT'] = [x/60 for x in im.get_ic_at_mass(self.Ions[0]). time_list]
        self.IntensityDFF = self.IntensityDF[self.RTadj[0]:self.RTadj[1]].nlargest(1, list(range(self.Ions[0], self.Ions[1]+1)))
        self.FinalRT = self.IntensityDFF[['RT']].values[0][0]

        self.RawAbun = self.IntensityDFF.values[0][:-1]
        self.chemicalFormula = chemicalFormula
        self.mainC = mainC
        self.labeledspecies = labeledspecies
        self.CorrRaw, self.enrichment = completoData(completoCM(self, labeledspecies), self.RawAbun)
        self.enrichmentFixed = [0 if x < 0 else x for x in self.enrichment]
        self.TotalIonCount = {self.Name: sum(self.RawAbun)}
        self.CorrFixed = [0 if x < 0 else x for x in self.CorrRaw]

        # CorrF is the total Ion abundance for the metabolite in question
        self.CorrF = dict([(self.Name[0:-3] + str(x) + ' (M' + str(y) + ')', self.CorrFixed[y]) for y, x in enumerate(self.IonL)])
        self.enrichmentF = dict([(self.Name[0:-3] + str(x)+' (M' + str(y)+')', self.enrichmentFixed[y]) for y, x in enumerate(self.IonL)])

        # self.enrichmentNormalized is the fragment information of the metabolite corrected so that they add up to 1.
        self.enrichmentNormalized = {self.Name[0:-3] + str(x)+' (M' + str(y)+')': self.enrichmentFixed[y]/sum(self.enrichmentFixed) for y, x in enumerate(self.IonL)}

        Item.instances.append(self)

    def __str__(self):
        return self.Name


listofCDFs = find_ext('.', 'CDF')

listofGCMSobs = [ANDI_reader(x) for x in listofCDFs]


def imSurveyor(GCMSob):
    im = build_intensity_matrix_i(GCMSob)

    #  Loading Control
    Nor260 = Item('Nor_260', 16.93, [260, 260], 'C12H30O1N1Si2', '', [], im)
    Nor288 = Item('Nor_288', 16.93, [288, 288], 'C13H30O2N1Si2', '', [], im)

    PYR174 = Item('PYR_174', 7.78, [174, 177], 'C6H12O3N1Si1', 'C3', ['C'], im)
    LAC261 = Item('LAC_261', 12.13, [261, 264], 'C10H25O2Si2', 'C3', ['C'], im)

    # #Amino Acids
    ALA260 = Item('ALA_260', 13.34, [260, 263], 'C11H26O2N1Si2', 'C3', ['C'], im)
    GLY246 = Item('GLY_246', 14.57, [246, 248], 'C10H24O2N1Si2', 'C2', ['C'], im)
    VAL288 = Item('VAL_288', 16.52, [288, 293], 'C13H30O2N1Si2', 'C5', ['C'], im)
    LEU344 = Item('LEU_344', 17.52, [344, 350], 'C17H38N1O2Si2', 'C6', ['C'], im)
    ILE344 = Item('ILE_344', 18.62, [344, 350], 'C17H38N1O2Si2', 'C6', ['C'], im)
    SER432 = Item('SER_432', 25.56, [432, 435], 'C20H46N1O3Si3', 'C3', ['C'], im)
    MET320 = Item('MET_320', 27.18, [320, 325], 'C13H30N1O2Si2S1', 'C5', ['C'], im)
    THR404 = Item('THR_404', 26.14, [404, 408], 'C18H42O3N1Si3', 'C4', ['C'], im)
    THR446 = Item('THR_446', 26.14, [446, 450], 'C21H48N1O3Si3', 'C4', ['C'], im)
    PRO328 = Item('PRO_328', 20.90, [328, 332], 'C16H36O2N1Si2', 'C4', ['C'], im)
    PHE336 = Item('PHE_336', 30.58, [336, 345], 'C17H30O2N1Si2', 'C9', ['C'], im)
    ASP418 = Item('ASP_418', 30.97, [418, 422], 'C18H40O4N1Si3', 'C4', ['C'], im)
    CYS406 = Item('CYS_406', 32.45, [406, 409], 'C17H40N1O2S1Si3', 'C3', ['C'], im)
    GLU432 = Item('GLU_432', 34.1, [432, 437], 'C19H42O4N1Si3', 'C5', ['C'], im)
    ASN417 = Item('ASN_417', 35.67, [417, 421], 'C18H41N2O3Si3', 'C4', ['C'], im)
    LYS488 = Item('LYS_488', 35.91, [488, 494], 'C24H56N2O2Si3', 'C6', ['C'], im)
    GLN431 = Item('GLN_431', 38.67, [431, 436], 'C19H43O3N2Si3', 'C5', ['C'], im)
    GLN473 = Item('GLN_473', 38.67, [473, 478], 'C22H49N2O3Si3', 'C5', ['C'], im)
    HIS440 = Item('HIS_440', 43.51, [440, 446], 'C20H42N3O2Si3', 'C6', ['C'], im)
    TYR466 = Item('TYR_466', 43.02, [466, 475], 'C23H44N1O3Si3', 'C9', ['C'], im)

    # #TCA
    SUC331 = Item('SUC_331', 21.50, [331, 336], 'C15H31O4Si2', 'C5', ['C'], im)
    SUC289 = Item('SUC_289', 21.50, [289, 294], 'C12H25O4Si2', 'C5', ['C'], im)

    FUM287 = Item('FUM_287', 21.93, [287, 291], 'C12H23O4Si2', 'C4', ['C'], im)
    AKG346 = Item('AKG_346', 29.35, [346, 351], 'C14H28O5N1Si2', 'C5', ['C'], im)
    MAL419 = Item('MAL_419', 29.86, [419, 423], 'C19H39O5Si3', 'C4', ['C'], im)
    CIT591 = Item('CIT_591', 41.52, [591, 597], 'C26H55O7Si4', 'C6', ['C'], im)


#  -----> Add them to this list

    Metabolites = Item.instances

    enrichmentDict = {}

    TIC_dict = {}
    ICs_dict = {}
    RT_dict = {}
    Ion_dict = {}

    for x in Metabolites:
        for k, v in x.enrichmentNormalized.items():
            enrichmentDict[k] = v
        for k, v in x.TotalIonCount.items():
            TIC_dict[k] = v
        for k, v in x.CorrF.items():
            ICs_dict[k] = v
        Ion_dict[x.Name+' (M0)'] = len(x.IonL) 
        RT_dict[x.Name] = x.FinalRT

    return enrichmentDict, TIC_dict, Ion_dict, ICs_dict, RT_dict, Metabolites


listofCDFsClean = [x.strip('./CDF') for x in listofCDFs]

enrichmentDict = {}
TIC_dict = {}
ICs_dict = {}
RT_dict = {}
Ion_dict = {}

for GCMSobs, CDFFile in zip(listofGCMSobs, listofCDFsClean):
    (enrichmentDict[CDFFile], TIC_dict[CDFFile],
        Ion_dict[CDFFile], ICs_dict[CDFFile], RT_dict[CDFFile], Metabolites) = imSurveyor(GCMSobs)

enrichment_DF = pd.DataFrame(enrichmentDict)
TIC_DF = pd.DataFrame(TIC_dict)
ICs_DF = pd.DataFrame(ICs_dict)
RT_DF = pd.DataFrame(RT_dict)
Ion_DF = pd.DataFrame(Ion_dict)

listOfDFs = (enrichment_DF, TIC_DF, ICs_DF, RT_DF, Ion_DF)
listOfDFNames = ('Ion Enrichment', 'Total Ion Counts', 'Ion Counts', 'Retention Times', 'Number of Ions')


# Graphing
def plotter(df, title, Metabolite, stacked):
    """
    This function can plot a DF with a title.

    plotter(TICCDF, 'Total Ion Count Pyruvate', Pyr_174, stacked=True)
    """
    if df is enrichment_DF:
        AAchart = grabber(df, Ion_dict[[x.strip('./CDF') for x in listofCDFs][0]], Metabolite)
        AAchart.T.plot(kind='bar', stacked=stacked)  # , yerr = AAcharterror)
        plt.title(Metabolite+' '+str(title)+' Enrichment')
        plt.show()
    elif df is TIC_DF:
        AAchart = pd.DataFrame(df.loc[Metabolite])
        AAchart.plot(kind='bar', stacked=stacked)  # , yerr = AAcharterror)
        plt.title(Metabolite+' '+str(title)+' TIC')
        plt.show()
    else:
        print('error')


i = 0
allMetabolites = [x.Name for x in Metabolites]
print(allMetabolites)

while i == 0:

    whichmetabos = input('Which metabolites would you like to plot? \
                            "i.e. PYR_174 or PYR_174, ALA_260"\n or exit: \n\n')
    try:
        if whichmetabos == 'exit':
            i += 1
        elif len(whichmetabos.split(',')) > 1:
            for x in [q.strip(' ') for q in whichmetabos.split(',')]:
                plotter(enrichment_DF, 'Metabolite Ion Abundance', x+' (M0)', True)
                plotter(TIC_DF, 'Metabolite Ion Abundance', x, True)
        elif whichmetabos == 'all':
            for x in allMetabolites:
                plotter(enrichment_DF, 'Metabolite Ion Abundance', x+' (M0)', True)
                plotter(TIC_DF, 'Metabolite Ion Abundance', x, True)
        else:
            plotter(enrichment_DF, 'Metabolite Ion Abundance', whichmetabos+' (M0)', True)
            plotter(TIC_DF, 'Metabolite Ion Abundance', whichmetabos, True)

    except KeyError:
        print('Invalid entry. Please try again.\n')

wb = xl.Workbook()

for DF, DFName in zip(listOfDFs, listOfDFNames):
    wb.create_sheet(DFName)
    ws = wb[DFName]
    for r in dataframe_to_rows(DF, index=True, header=True):
        ws.append(r)

p = 0
while p == 0:
    save = input("Would you like to save your data as a xlsx? (Y/N)\n")

    match save:
        case "Y":
            filename = input('What would you like to save your file as?\n')
            wb.save(filename+'.xlsx')
            print(filename+".xlsx has been saved in the local directory\n")
            p += 1
        case "N":
            print("Goodbye")
            p += 1
        case _:
            print('Invalid entry. Please try again')
