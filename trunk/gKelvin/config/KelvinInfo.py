"""Contains values and constants used by other parts of the program"""

# a simple struct used as a container
class container(object):
    pass

# the actual program command used to run kelvin
kelvinRunCommand = 'kelvin'

# the different possible filename tags, description, and default value
# format is <tag> <description> <default value> <True if required, or Falset>
fileNameEntry = [('PD', 'Pedigree file (post-makeped)', 'pedfile.dat', True), 
                 ('DF', 'Data file', 'datafile.dat', True), 
                 ('MK', 'Marker file', 'markers.dat', True),
                 ('MP', 'Map file', 'mapfile.dat', True), 
                 ('HE', 'Output: Bayes Ratio', 'br.out', False), 
                 ('PF', 'Output: PPL', 'ppl.out', False)]

def getFileNameTags():
    return [ x[0] for x in fileNameEntry ]

def getFileNameTagDescriptions():
    return [ x[1] for x in fileNameEntry ]

def getFileNameDefaults():
    return [ x[2] for x in fileNameEntry ]

def getFileNameInfo():
    """Returns information about the filename tags, descriptions, and
    default value
    
    """

    entries = fileNameEntry
    result = []

    for eachEntry in entries:
        temp = container()
        temp.tag = eachEntry[0]
        temp.description = eachEntry[1]
        temp.default = eachEntry[2]
        temp.isRequired = eachEntry[3]
        result.append(temp)

    return result

fileNameInfo = getFileNameInfo()

# the possible analysis types, 
# format for each entry is <tag> <description> <# of arguments after tag>
analysisTypeEntry = [('TP', 'two point', 0), 
                     ('SS', 'multi-point sex specific', 1), 
                     ('SA', 'multi-point sex average', 1)]

def getAnalysisTypeTags():
    return [ x[0] for x in analysisTypeEntry ]

def getAnalysisTypeDescriptions():
    return [ x[1] for x in analysisTypeEntry ]

def getAnalysisTypeNumbers():
    return [ x[2] for x in analysisTypeEntry ]

def getAnalysisTypeInfo():
    """Returns information about the analysis type tags, descriptions, and
    number of items after it in a struct"""

    entries = analysisTypeEntry
    result = []

    for eachEntry in entries:
        temp = container()
        temp.tag = eachEntry[0]
        temp.description = eachEntry[1]
        temp.numItems = eachEntry[2]
        result.append(temp)

    return result

analysisTypeTags = getAnalysisTypeTags()

# format is <tag> <number of items after tag> <default value>
unknownPersonId = ('UP', 1, 0)
diseaseAlleles = ('DA', 1, 2)

genotypes = ['DD', 'Dd', 'dd']

genotypeTagToDescription = dict(zip(genotypes, genotypes))
genotypeDescriptionToTag = dict(zip(genotypes, genotypes))
# this is needed when constraint choice boxes are blanked out
genotypeTagToDescription[''] = ''
genotypeDescriptionToTag[''] = ''

# format is <tag> <description> <# of numbers after it>
# the ? is there because it doesn't fit the pattern
traitTypesEntry = [('DT', 'Dichotomous', 0), ('QT', 'Quantitative', '?'),
                   ('QT', 'QT Threshold', '?')]

def getTraitTypeTags():
    return [ x[0] for x in traitTypesEntry ]

def getTraitTypeDescriptions():
    return [ x[1] for x in traitTypesEntry ]

traitTypeTags = getTraitTypeTags()
traitTypeDescriptions = getTraitTypeDescriptions()
traitTypeTagToDescription = dict(zip(traitTypeTags, traitTypeDescriptions))
traitTypeDescriptionToTag = dict(zip(traitTypeDescriptions, traitTypeTags))

# format is <tag> <# of numbers after tag>
affectionStatusEntry = ('AS', 3)

# defaults for affection status
# format is <unknown> <unaffected> <affected>
defaultAffectionDT = (0, 1, 2)
defaultAffectionQT = (-99.99, -88.88, 88.88)

# distributions used for QT
# format is <tag> <description> <# of params>
#distributionEntry = [('T', 'T-Distribution', 3), ('normal', 'Normal', 2),
#                     ('chisq', 'Chi Squared', 2)]
distributionEntry = [('T', 'Normal', 3), ('chisq', 'Chi Squared', 2)]

def getDistributionTags():
    return [ x[0] for x in distributionEntry ]

def getDistributionDescriptions():
    return [ x[1] for x in distributionEntry ]

distributionTags = getDistributionTags()
distributionDescriptions = getDistributionDescriptions()
distributionTagToDescription = dict(zip(distributionTags, 
                                        distributionDescriptions))
distributionDescriptionToTag = dict(zip(distributionDescriptions, 
                                        distributionTags))

# default parameters for different distributions
defaultDistributionNormal = ('0', '1')
defaultDistributionT = ('0', '1', '30')
defaultDistributionChiSquared = ('0', '1')

# format is <tag> <# of numbers after tag>
liabilityClassEntry = ('LC', 1)
defaultLiabilityClass = '1'

# the possible operators in constraints
# format is <operator> <description>
operators = [('==', 'is equal to'), ('!=', 'is not equal to'), 
             ('>', 'is greater than'), ('>=', 'is greater than or equal to')]

def getOperatorTags():
    return [ x[0] for x in operators ]

def getOperatorDescriptions():
    return [ x[1] for x in operators ]

operatorTags = getOperatorTags()
operatorDescriptions = getOperatorDescriptions()
operatorTagToDescription = dict(operators)
operatorDescriptionToTag = dict(zip(operatorDescriptions, operatorTags))
# this is needed when constraint choice boxes are blanked out
operatorTagToDescription[''] = ''
operatorDescriptionToTag[''] = ''

# sex-specific thetas, format is <tag> <description>
thetaConstraintParameters = [('Tm', 'theta-male'), ('Tf', 'theta-female')]

def getThetaTags():
    return [ x[0] for x in thetaConstraintParameters ]

def getThetaDescriptions():
    return [ x[1] for x in thetaConstraintParameters ]

thetaTags = getThetaTags()
thetaDescriptions = getThetaDescriptions()
thetaTagToDescription = dict(thetaConstraintParameters)
thetaDescriptionToTag = dict(zip(thetaDescriptions, thetaTags))
# this is needed when constraint choice boxes are blanked out
thetaTagToDescription[''] = ''
thetaDescriptionToTag[''] = ''

# standard deviation parameter, format is <tag> <description>
stdParameter = ('P1', 'standard deviation')

# this is a list of all the variable number tags used, which
# also contains information 
# format is <tag> <description> <depends on analysis type> 
# <available for two-point> <available for multi-point> 
# <depends on trait type> <available for dt> <available for qt>
variableNumberEntries = \
    [('AF', 'marker allele frequency', False, True, True, False, True, True),
     ('GF', 'disease allele frequency ', False, True, True, False, True, True),
     ('Th', 'theta', True, True, False, False, True, True),
     ('Tm', 'theta-male', True, True, False, False, True, True),
     ('Tf', 'theta-female', True, True, False, False, True, True),
     ('TL', 'trait location', True, False, True, False, True, True),
     ('DD', 'DD', False, True, True, True, True, True),
     ('Dd', 'Dd', False, True, True, True, True, True),
     ('dd', 'dd', False, True, True, True, True, True),
     ('P1', 'standard deviation', False, True, True, True, False, True),
     ('TT', 'trait threshold', False, True, True, True, False, True),
     ('AL', 'heterogeneity (alpha)', False, True, True, False, True, True),
     ('LD', "linkage disequilibrium D'", True, True, False, False, True, True),
    ]

def getVariableNumberInfo():
    """Gets all the information variable number tags in a struct construct"""

    entries = variableNumberEntries
    result = []

    for eachEntry in entries:
        temp = container()
        temp.tag = eachEntry[0]
        temp.description = eachEntry[1]
        temp.analysisType = eachEntry[2]
        temp.tp = eachEntry[3]
        temp.mp = eachEntry[4]
        temp.traitType = eachEntry[5]
        temp.dt = eachEntry[6]
        temp.qt = eachEntry[7]
        result.append(temp)

    return result

variableNumberInfo = getVariableNumberInfo()
variableNumberTagToDescription = dict((x.tag, x.description) for x in 
                                                        variableNumberInfo)
variableNumberDescriptionToTag = dict((x.description, x.tag) for x in 
                                                        variableNumberInfo)

# depending on the trait type, the genotype description should be
# preceded by a word describing it, for DT it should say 'penetrance dd', etc.,
# and for QT it should say 'mean dd', etc.
genotypePrefix = ['penetrance ', 'mean ']
genPrefix = [(pre + gen, gen) for pre in genotypePrefix for gen in genotypes]

# update the dictionary to include the new entries
variableNumberDescriptionToTag.update(genPrefix)

# format is <tag> <# of numbers after tag>
polynomialEvaluationEntry = ('PE', 0)

# format is <tag> <# of numbers after tag>
TMEntry = ('TM', 0)
MMEntry = ('MM', 0)
AMEntry = ('AM', 0)
chromosomeEntry = ('XC', 0)
MINEntry = ('MIN', 1)
MAXEntry = ('MAX', 1)
TMINEntry = ('T_MIN', 1)
TMAXEntry = ('T_MAX', 1)

# list of types of new kelvin config files
# format is <description> <filename>
newConfigFile = [('Two Point', 'twopoint_dt.conf'),
                 ('Multi-point', 'multipoint_dt.conf')]

# some tags don't really fit into any template, so i'll just
# put them here and use them when needed. Create the line
# as a list of parsed strings.
defaultTT = ['TT', '1', ';', '2', ';', '3']
defaultLD = ['LD', '-1', '1', '0.1']

# a priority list of all the tags, used to determine where to put new file
# lines. will be searched to see if the first token matches these tags.
# format is <tag> <comment> <whether this line is a constraint>
priorityTagList = \
    [('PD', 'pedigree file (post-makeped)', False), 
     ('DF', 'data file (list of loci)', False), 
     ('MK', 'marker file (marker allele frequencies)', False), 
     ('MP', 'map file', False), 
     ('OF', 'full likelihood set output file', False), 
     ('HE', 'bayes ratio output file', False), 
     ('PF', 'ppl output file', False), 
     ('TP', 'two-point analysis', False), 
     ('SS', 'multi-point sex specific analysis', False), 
     ('SA', 'multi-point sex average analysis', False), 
     ('LD', 'linkage disequilibrium', False), 
     ('XC', 'X chromosome Analysis', False), 
     ('LC', 'number of liability classes', False), 
     ('DT', 'dichotomous trait', False), 
     ('QT', 'quantitative trait', False), 
     ('AS', 'phenotype codes (affection status)', False), 
     ('MIN', 'min for distribution', False), 
     ('MAX', 'max for distribution', False), 
     ('T_MIN', 'min for distribution', False), 
     ('T_MAX', 'max for distribution', False), 
     ('MM', 'marker to marker analysis (all pairwise)', False), 
     ('AM', 'marker to marker analysis (only adjacent)', False), 
     ('UP', 'unknown person id in pedigree file', False), 
     ('PE', 'using polynomial evaluations to speed up calculations', False), 
     ('TM', 'Calculation will be done at marker locations also', False), 
     ('DD', 'penetrance or mean value', False), 
     ('Dd', 'penetrance or mean value', False), 
     ('dd', 'penetrance or mean value', False), 
     ('DD', 'penetrance or mean constraint', True), 
     ('Dd', 'penetrance or mean constraint', True), 
     ('dd', 'penetrance or mean constraint', True), 
     ('Tm', 'theta constraint', True), 
     ('Tf', 'theta constraint', True), 
     ('P1', 'standard deviation constraint', True), 
     ('GF', 'gene frequency', False), 
     ('AF', 'marker allele frequency', False), 
     ('Th', 'theta', False), 
     ('Tm', 'theta-male', False), 
     ('Tf', 'theta-female', False), 
     ('TL', 'trait location positions', False), 
     ('AL', 'alpha values', False), 
     ('P1', 'standard deviation values', False), 
     ('TT', 'trait threshold', False)
    ]
