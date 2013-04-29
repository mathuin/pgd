import unittest
import datetime
from selenium import webdriver
from pgd_search.models import *
from pgd_core.models import *
#from pgd_splicer.SegmentBuilder import SegmentBuilderTask
from pgd_constants import AA_CHOICES, SS_CHOICES
from math import ceil
from search.SearchForm import SearchSyntaxField

PRO_MIN = -1
PRO_MAX = 3
FIELDS = ['a1','a2','a3','a4','a5','a6','a7','L1','L2','L3','L4','L5','phi','psi','ome','chi1','chi2','chi3','chi4','chi5','bm','bs','bg','h_bond_energy','zeta']
FIELDS_DICT = {}
for i in range(1, len(FIELDS)+1):
    #shift values into decimels
    FIELDS_DICT[FIELDS[i-1]] = i*.01

class SearchParserValidation(unittest.TestCase):

    def calculateAA(self, chainIndex):
        #aa_choice = chainIndex-1 if chainIndex-1 < len(AA_CHOICES) else chainIndex-1-len(AA_CHOICES)
        return AA_CHOICES[(chainIndex-1)%len(AA_CHOICES)][0]

    def calculateSS(self, chainIndex):
        #ss_choice = chainIndex-1 if chainIndex-1 < len(SS_CHOICES) else chainIndex-1-len(SS_CHOICES)
        return SS_CHOICES[(chainIndex-1)%len(SS_CHOICES)][0]

    #calculates a value for a field based on other properties
    #   1) whole number indicates protein id
    #   2) first 2 digits after decimal indicate field in FIELDS lookup list/dict
    #   3) last 2 digits indicate chainIndex
    #
    #   example:  1.0305  =>  1.  03   05  =>  protein=1  field=03   chainindex=05
    def calculateResidueField(self, proteinID, chainIndex, field):
        return proteinID + FIELDS_DICT[field] + (chainIndex*.0001)

    def create_protein(self, i, **kwargs):
        protein = Protein()
        protein.code            = '%+i' %i
        protein.threshold       = i
        protein.resolution      = i + .01
        protein.rfactor         = i + .02
        protein.rfree           = i + .03
        protein.pdb_date        = datetime.date(2001,1,1)
        protein.__dict__.update(kwargs)
        protein.save()
        return protein

    def create_chain(self, protein):
        chain = Chain()
        chain.id = '%s%s' % (protein.code, 'A')
        chain.protein = protein
        chain.code = 'A'
        chain.save()
        return chain

    def create_residue(self, i, z, protein, chain, **kwargs):
        residue = Residue()
        residue.protein = protein
        residue.chain = chain
        residue.chainID = chain.code
        residue.chainIndex = z
        #set choice fields.
        residue.aa = self.calculateAA(i)
        residue.ss = self.calculateSS(i)
        for field in FIELDS:
            residue.__dict__[field] = self.calculateResidueField(i, z, field)
        residue.__dict__.update(kwargs)
        residue.save()
        return residue

    def create_bulk_residues(self, i, numResidue, **kwargs):
        protein = self.create_protein(i, **kwargs)
        chain = self.create_chain(protein)
        chainList = [self.create_residue(i, 1, protein, chain) for z in range(numResidue)]
        return chainList


    #creates a set objects with predictable values so we can predict search results
    def setUp(self):
        self.tearDown()

    def tearDown(self):
        Protein.objects.all().delete()

    def testSearchQueryStrings(self):

        # create protein
        protein = self.create_protein(1)
        #create chain
        chain = self.create_chain(protein)
        #create residues
        chainList = [self.create_residue(1, z, protein, chain, a1=z) for z in range(6)]

        # create Search
        search = Search(segmentLength=5)

        # create associated Search_residues
        data = {}
        #data['index'] = 5
        data['residues'] = 1
        data['a1_i_0'] = 1
        data['a1_0'] = "1-5"

        #Set search data equal to search residue parameters
        search.data = data
        search.save()

        self.assertEqual(
            # See that the intended query is executed by parse_search; upper bound of the expected range is NON-INCLUSIVE
            set(chainList[1:6]),
            set(Search.parse_search(search)),
            "Query strings search test failed on range '%s'.  Expected value was %s    Actual values were %s"%(data['a1_0'], set(chainList[1:5]), set(Search.parse_search(search)))
        )

        data['a1_0'] = "1-3"
        search.data = data
        search.save()
        self.assertEqual(
            # See that the intended query is executed by parse_search; upper bound of the expected range is NON-INCLUSIVE
            set(chainList[1:4]),
            set(Search.parse_search(search)),
            "Query strings search test failed on range '%s'"%data['a1_0']
        )


        data['a1_0'] = '<6'
        search.data = data
        self.assertEqual(
            # See that the intended query is executed by parse_search
            set(chainList[0:6]),
            set(Search.parse_search(search)),
            "Query strings search test failed on range '%s'.  Expected value was %s    Actual values were %s"%(data['a1_0'], set(chainList[0:6]), set(Search.parse_search(search)))
        )

        data['a1_0'] = '<=6'
        search.data = data
        self.assertEqual(
            # See that the intended query is executed by parse_search
            set(chainList[0:7]),
            set(Search.parse_search(search)),
            "Query strings search test failed on range '%s'.  Expected value was %s    Actual values were %s"%(data['a1_0'], set(chainList[0:7]), set(Search.parse_search(search)))
        )

        data['a1_0'] = '>5'
        search.data = data
        self.assertEqual(
            # See that the intended query is executed by parse_search
            set(chainList[6:7]),
            set(Search.parse_search(search)),
            "Query strings search test failed on range '%s'.  Expected value was %s    Actual values were %s"%(data['a1_0'], set(chainList[6:7]), set(Search.parse_search(search)))
        )

        data['a1_0'] = '>=5'
        search.data = data
        self.assertEqual(
            # See that the intended query is executed by parse_search
            set(chainList[5:7]),
            set(Search.parse_search(search)),
            "Query strings search test failed on range '%s'.  Expected value was %s    Actual values were %s"%(data['a1_0'], set(chainList[5:7]), set(Search.parse_search(search)))
        )


    def testSearchResolution(self):

        # create protein
        protein1 = self.create_protein(1, resolution=1.1)
        #create chain
        chain = self.create_chain(protein1)
        #create residues
        chainList = [self.create_residue(1, z, protein1, chain) for z in range(3)]

        # create protein
        protein2 = self.create_protein(2, resolution=0.2)
        #create chain
        chain2 = self.create_chain(protein2)
        #create residues
        chainList2 = [self.create_residue(2, i, protein2, chain2) for i in range(3)]

        #create Search and search constraints
        search = Search(segmentLength=0)
        data = {}
        data['residues'] = 0
        data['resolutionMin'] = .1
        data['resolutionMax'] = 1.2

        #Set search data equal to search residue parameters
        search.data = data
        #search.save()
        self.assertEqual(
            # See that both proteins are returned
            set(chainList+chainList2),
            set(Search.parse_search(search)),
            "Resolution search failed on %f-%f   Expected results were %s   Returned results were %s"%(data['resolutionMin'],data['resolutionMax'],list(x.protein.resolution for x in chainList+chainList2),list(y.protein.resolution for y in Search.parse_search(search)))
        )

        data['resolutionMin'] = 0.1
        data['resolutionMax'] = 0.3
        #Set search data equal to search residue parameters
        search.data = data
        #search.save()
        self.assertEqual(
            # See that one protein is returned
            set(chainList2),
            set(Search.parse_search(search)),
            "Resolution search failed on %f-%f   Expected results were %s   Returned results were %s"%(data['resolutionMin'],data['resolutionMax'],list(x.protein.resolution for x in chainList2),list(y.protein.resolution for y in Search.parse_search(search)))
        )

        data['resolutionMin'] = 1
        data['resolutionMax'] = 1.2
        #Set search data equal to search residue parameters
        search.data = data
        search.save()
        self.assertEqual(
            # See that one protein is returned
            set(chainList),
            set(Search.parse_search(search)),
            "Resolution search failed on %f-%f   Expected results were %s   Returned results were %s"%(data['resolutionMin'],data['resolutionMax'],list(x.protein.resolution for x in chainList),list(y.protein.resolution for y in Search.parse_search(search)))
        )

        data['resolutionMin'] = 1.1
        data['resolutionMax'] = 1.1
        #Set search data to exactly the second protein's resolution
        search.data = data
        search.save()
        self.assertEqual(
            # See that one protein is returned
            set(chainList),
            set(Search.parse_search(search)),
            "Resolution search failed on %f-%f   Expected results were %s   Returned results were %s"%(data['resolutionMin'],data['resolutionMax'],list(x.protein.resolution for x in chainList),list(y.protein.resolution for y in Search.parse_search(search)))
        )

        data['resolutionMin'] = .3
        data['resolutionMax'] = .1
        #Set search data to exactly the second protein's resolution
        search.data = data
        search.save()
        self.assertEqual(
            # See that no proteins are returned
            set(),
            set(Search.parse_search(search)),
            "Resolution search failed on %f-%f   Expected results were %s   Returned results were %s"%(data['resolutionMin'],data['resolutionMax'],set(),list(y.protein.resolution for y in Search.parse_search(search)))
        )


    def testSearchThreshold(self):

        #locals = locals()
        chainList = {}
        for t in zip(range(1,6), (1,10,25,25,70)):
            i,thresh = t
            chainList[i] = self.create_bulk_residues( i, 1, threshold = thresh)

        #setattr(obj, name, value)

        # create Search
        search = Search(segmentLength=0)
        data = {}
        data['residues'] = 0
        data['threshold'] = 0
        search.data = data

        self.assertEqual(
            # See that the intended query is executed by parse_search
            set(),
            set(Search.parse_search(search)),
            "Threshold search failed on %i  Expected results were the %s  Returned results were %s"%(data['threshold'],set(), list(y.protein.threshold for y in Search.parse_search(search)))
        )

        data['threshold'] = 25
        search.data = data
        self.assertEqual(
            # See that the intended query is executed by parse_search
            set(chainList[1]+chainList[2]+chainList[3]+chainList[4]),
            set(Search.parse_search(search)),
            "Threshold search failed on %i  Expected results were the %s  Returned results were %s"%(data['threshold'],list(y.protein.threshold for y in chainList[1]+chainList[2]+chainList[3]+chainList[4]), list(y.protein.threshold for y in Search.parse_search(search)))
        )

        data['threshold'] = 75

        self.assertEqual(
            # See that the intended query is executed by parse_search
            set(chainList[1]+chainList[2]+chainList[3]+chainList[4]+chainList[5]),
            set(Search.parse_search(search)),
            "Threshold search failed on %i  Expected results were the %s  Returned results were %s"%(data['threshold'],list(y.protein.threshold for y in chainList[1]+chainList[2]+chainList[3]+chainList[4]+chainList[5]), list(y.protein.threshold for y in Search.parse_search(search)))
        )


        #for index in range(PRO_MIN,PRO_MAX):
            #search.threshold = index
            #search.save()

            #self.assertEqual(
                ## See that the intended query is executed by parse_search
                #set(Segment.objects.filter(protein__threshold=index).all()),
                #set(Search.parse_search(search).all()),
                #"Threshold search failed on %i"%(index)
            #)

    def testSearchCode(self):
        data = {}
        #create residues
        chainList = self.create_bulk_residues( 1, 3, code='gmez')
        chainList2 = self.create_bulk_residues( 2, 3, code='kats')
        # create Search
        search = Search(segmentLength=0)
        search.codes_include = True
        data['proteins_i'] = True
        data['residues'] = 0
        data['proteins'] = 'gmez'
        search.data = data

        self.assertEqual(
            set(chainList),
            set(Search.parse_search(search)),
            "Resolution search failed on %s   Expected results were %s   Returned results were %s"%(data['proteins'],list(x.protein.code for x in chainList),list(y.protein.code for y in Search.parse_search(search)))
        )

        #test inverse
        data['proteins_i'] = False
        self.assertEqual(
            set(chainList2),
            set(Search.parse_search(search)),
            "Resolution search failed on %s   Expected results were %s   Returned results were %s"%(data['proteins'],list(x.protein.code for x in chainList),list(y.protein.code for y in Search.parse_search(search)))
        )

        #testing multi-code and parsing robustness
        data['proteins'] = 'gmez,             kats'
        data['proteins_i'] = True
        self.assertEqual(
            set(chainList+chainList2),
            set(Search.parse_search(search)),
            "Resolution search failed on %s   Expected results were %s   Returned results were %s"%(data['proteins'],list(x.protein.code for x in chainList+chainList2),list(y.protein.code for y in Search.parse_search(search)))
        )

        #test inverse
        data['proteins_i'] = True
        data['proteins'] = 'none'
        self.assertEqual(
            set(),
            set(Search.parse_search(search)),
            "Resolution search failed on the empty set   Expected results were the empty set   Returned results were %s"%(list(y.protein.code for y in Search.parse_search(search)))
        )

        #test inverse
        data['proteins_i'] = False
        data['proteins'] = 'gmez,kats'
        self.assertEqual(
            set(),
            set(Search.parse_search(search)),
            "Resolution search failed on the empty set   Expected results were the empty set   Returned results were %s"%(list(y.protein.code for y in Search.parse_search(search)))
        )

    def testSearchAa(self):
        aa_range = ('a', 'r', 'n', 'd', 'c')
        chainList = {}
        for t in zip(range(1,6), aa_range):
            i,aa_choice = t
            chainList[i] = self.create_bulk_residues( i, 1, aa = aa_choice)

        # create Search
        search = Search(segmentLength=0)
        #search.save()
        data = {}
        data['residues'] = 1

        for aa_index,aa_choice in enumerate(aa_range):
            data['aa_i_0'] = 1
            data['aa_0'] = str(aa_choice)
            search.data = data

            self.assertEqual(
                # See that the intended query is executed by parse_search
                list(chainList[aa_index+1]),
                list(Search.parse_search(search)),
                "Specific AA search failed on '%s'  Expected results were %s  Returned results were %s with %s entries"%(aa_choice, (list(getattr(x, 'aa') for x in (chainList[aa_index+1]))), (list(str(getattr(y, 'aa')) for y in (Search.parse_search(search)))), len(Search.parse_search(search)))
            )

            #test the negation of a single index.
            data['aa_i_0'] = 0
            negatedList = []
            for k in range(1,6):
                if k != (aa_index+1):
                    negatedList += chainList[k]
            search.data = data
            self.assertEqual(
                list(negatedList),
                list(Search.parse_search(search)),
                "Negated AA search failed on '%s'  Expected results were %s  Returned results were %s with %s entries"%(aa_choice, (list(getattr(x, 'aa') for x in (negatedList))), (list(str(getattr(y, 'aa')) for y in (Search.parse_search(search)))), len(Search.parse_search(search)))
            )

        data['aa_i_0'] = 1
        aa_choice = ('a', 'r', 'n')
        data['aa_0'] = (aa_choice)
        search.data = data
        #data['aa_i'] = 0
        self.assertEqual(
            # See that the intended query is executed by parse_search
            list(chainList[1] + chainList[2] + chainList[3]),
            list(Search.parse_search(search)),
            "Specific AA search failed on '%s'  Expected results were %s  Returned results were %s with %s entries"%(aa_choice, (list(getattr(x, 'aa') for x in (chainList[1] + chainList[2] + chainList[3]))), (list(str(getattr(y, 'aa')) for y in (Search.parse_search(search)))), len(Search.parse_search(search)))
        )

        data['aa_i_0'] = 0
        aa_choice = ('a', 'r', 'n')
        data['aa_0'] = (aa_choice)
        search.data = data
        #data['aa_i'] = 0
        self.assertEqual(
            # See that the intended query is executed by parse_search
            list(chainList[4] + chainList[5]),
            list(Search.parse_search(search)),
            "Specific AA search failed on '%s'  Expected results were %s  Returned results were %s with %s entries"%(aa_choice, (list(getattr(x, 'aa') for x in (chainList[4] + chainList[5]))), (list(str(getattr(y, 'aa')) for y in (Search.parse_search(search)))), len(Search.parse_search(search)))
        )


    def testSearchSs(self):

        ss_range = ('H', 'G',  'E', 'T', 'S')
        chainList = {}
        for t in zip(range(1,6), ss_range):
            i,ss_choice = t
            chainList[i] = self.create_bulk_residues( i, 1, ss = ss_choice)

        # create Search
        search = Search(segmentLength=0)
        #search.save()
        data = {}
        data['residues'] = 1

        for ss_index,ss_choice in enumerate(ss_range):
            data['ss_i_0'] = 1
            data['ss_0'] = ss_choice
            search.data = data

            self.assertEqual(
                # See that the intended query is executed by parse_search
                list(chainList[ss_index+1]),
                list(Search.parse_search(search)),
                "Specific ss search failed on '%s'  Expected results were %s  Returned results were %s with %s entries"%(ss_choice, (list(getattr(x, 'ss') for x in (chainList[ss_index+1]))), (list(str(getattr(y, 'ss')) for y in (Search.parse_search(search)))), len(Search.parse_search(search)))
            )

            #test the negation of a single index.
            data['ss_i_0'] = 0
            negatedList = []
            for k in range(1,6):
                if k != (ss_index+1):
                    negatedList += chainList[k]
            search.data = data
            self.assertEqual(
                list(negatedList).sort(),
                list(Search.parse_search(search)).sort(),
                "Negated ss search failed on '%s'  Expected results were %s  Returned results were %s with %s entries"%(ss_choice, (list(getattr(x, 'ss') for x in (negatedList))), (list(str(getattr(y, 'ss')) for y in (Search.parse_search(search)))), len(Search.parse_search(search)))
            )

        data['ss_i_0'] = 1
        ss_choice = ('H', 'G',  'E')
        data['ss_0'] = (ss_choice)
        search.data = data
        #data['ss_i'] = 0
        self.assertEqual(
            # See that the intended query is executed by parse_search
            list(chainList[1] + chainList[2] + chainList[3]).sort(),
            list(Search.parse_search(search)).sort(),
            "Specific ss search failed on '%s'  Expected results were %s  Returned results were %s with %s entries"%(ss_choice, (list(getattr(x, 'ss') for x in (chainList[1] + chainList[2] + chainList[3]))), (list(str(getattr(y, 'ss')) for y in (Search.parse_search(search)))), len(Search.parse_search(search)))
        )

        data['ss_i_0'] = 0
        ss_choice = ('H', 'G',  'E')
        data['ss_0'] = (ss_choice)
        search.data = data
        #data['ss_i'] = 0
        self.assertEqual(
            # See that the intended query is executed by parse_search
            list(chainList[4] + chainList[5]).sort(),
            list(Search.parse_search(search)).sort(),
            "Specific ss search failed on '%s'  Expected results were %s  Returned results were %s with %s entries"%(ss_choice, (list(getattr(x, 'ss') for x in (chainList[4] + chainList[5]))), (list(str(getattr(y, 'ss')) for y in (Search.parse_search(search)))), len(Search.parse_search(search)))
        )

    def testSearchMultipleResidues(self):
        # create protein
        protein = self.create_protein(1)
        #create chain
        chain = self.create_chain(protein)
        #create residues
        chainList = [self.create_residue(1, z, protein, chain) for z in range(1,6)]

        # create Search
        search = Search(segmentLength=0)

        for i, field in enumerate(FIELDS):

            # create associated Search_residues
            data = {}
            data['residues'] = 1
            data['%s_i_-2'%field] = 1
            data['%s_-2'%field] = str(self.calculateResidueField(1, 1, field))
            data['%s_i_0'%field] = 1
            data['%s_0'%field] = str(self.calculateResidueField(1, 3, field))
            data['%s_i_2'%field] = 1
            data['%s_2'%field] = str(self.calculateResidueField(1, 5, field))
            search.data = data
            self.assertAlmostEqual(
                (getattr(chainList[2], '%s'%field)),
                (getattr((Search.parse_search(search).all()[0]), '%s'%field)),#this grabs the first object in the returned search, and attempts to grab the attribute in 'field' out of it
                4,
                "Multiple residue search failed on field %s  Expected result was %s  Returned result was %s"%(field, set((getattr(chainList[2], '%s'%field), )), set(getattr(x, '%s'%field) for x in Search.parse_search(search)))
            )

class SearchFieldValidationCase(unittest.TestCase):
    def setUp(self):
        pass

    def testFieldSyntaxParser(self):
        validFields = [
            '1',
            '1-1',
            '1,2,3',
            '1-2',

            '1-1',
            '1-2,5-6',
            '1,23,456',
            '1-23',
            '1-23,456-7890',
            '-1',
            '-23',
            '.5',
            '.123',
            '0.5',
            '0.5,0.6',
            '0.5-0.6',
            '0.5-0.6,0.5',
            '0.5-123.123',
            '-1-2',
            '1--2',
            '-1--2'
        ]

        invalidFields = []
        searchField = SearchSyntaxField()

        for value in validFields:
            self.assertNotEqual(searchField.syntaxPattern.match(value), None, "Valid Field Pattern Failed: '%s'" % value)

        for value in invalidFields:
            self.assertNotEqual(searchField.syntaxPattern.match(value), None)

#Selenium tests
class PersistingSearchOptions(unittest.TestCase):
    def setUp(self):

       # Create a new instance of the Firefox driver
        self.driver = webdriver.Firefox()

    def test_removed_options_persist(self):

        # Load search page
        self.driver.get("http://localhost:8000/search")

        # Select the box that indicates number of residues
        residues = self.driver.find_element_by_id("id_residues")
        for option in residues.find_elements_by_tag_name('option'):
            if option.text == "4":
                option.click()

        #Composition
        composition = self.driver.find_element_by_id("id_aa_choices_list_col_2")
        comp_option = composition.find_elements_by_tag_name('li')[0]
        comp_option.click()

        #Conformation
        conformation = self.driver.find_element_by_id("id_ss_choices_list_col_2")
        conf_option = conformation.find_elements_by_tag_name('li')[0]
        conf_option.click()

        #Mobility
        self.driver.find_element_by_id("mobility_header").click()
        mobility = self.driver.find_element_by_id("id_bm_2")
        mobility.clear()

        #Angles
        self.driver.find_element_by_id("angles_header").click()
        angles = self.driver.find_element_by_id("id_a1_2")
        angles.click()
        angles.send_keys("30")

        #Lengths
        self.driver.find_element_by_id("lengths_header").click()
        lengths = self.driver.find_element_by_id("id_L1_2")
        lengths.click()
        lengths.send_keys("25")

        #XAngles
        ele = self.driver.find_element_by_id("chi_header").click()
        xangle = self.driver.find_element_by_id("id_chi1_2")
        xangle.click()
        xangle.send_keys("20")

        #Hackish way to do it, but there doesn't seem to be any other
        #common ways to do it.
        residues = self.driver.find_element_by_id("id_residues")
        for option in residues.find_elements_by_tag_name('option'):
            if option.text == "3":
                option.click()

        for option in residues.find_elements_by_tag_name('option'):
            if option.text == "4":
                option.click()

        #Composition
        self.assertEquals(comp_option.get_attribute("class"), "selected")

        #Conformation
        self.assertEquals(conf_option.get_attribute("class"), "selected")

        #Mobility
        self.assertEquals(mobility.text, "<25")

        #Angles
        self.assertEquals(angles.text, "")

        #Lengths
        self.assertEquals(lengths.text, "")

        #XAngles
        self.assertEquals(xangles.text, "")
