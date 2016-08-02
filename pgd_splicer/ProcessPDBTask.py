#!/usr/bin/env python
from __future__ import division

if __name__ == '__main__':
    import sys
    import os

    #python magic to add the current directory to the pythonpath
    sys.path.append(os.getcwd())

    # ==========================================================
    # Setup django environment
    # ==========================================================
    if not 'DJANGO_SETTINGS_MODULE' in os.environ:
        os.environ['DJANGO_SETTINGS_MODULE'] = 'pgd.settings'
    # ==========================================================
    # Done setting up django environment
    # ==========================================================

from datetime import datetime
from django.utils import timezone
from gzip import GzipFile
from math import degrees, sqrt
from multiprocessing import Pool
from operator import itemgetter
import os
import sys
from tempfile import NamedTemporaryFile
import logging

import Bio.PDB
from Bio.PDB import calc_angle as pdb_calc_angle, PDBIO
from Bio.PDB import calc_dihedral as pdb_calc_dihedral
from Bio.PDB.PDBIO import Select
from django.db import transaction

from tools import localfile

from pgd_core.models import (Protein as ProteinModel, Chain as ChainModel,
                             Residue as ResidueModel, Sidechain_ARG,
                             Sidechain_ASN, Sidechain_ASP, Sidechain_CYS,
                             Sidechain_GLN, Sidechain_GLU, Sidechain_HIS,
                             Sidechain_ILE, Sidechain_LEU, Sidechain_LYS,
                             Sidechain_MET, Sidechain_PHE, Sidechain_PRO,
                             Sidechain_SER, Sidechain_THR, Sidechain_TRP,
                             Sidechain_TYR, Sidechain_VAL)

from pgd_splicer.chi import CHI_MAP, CHI_CORRECTIONS_TESTS, CHI_CORRECTIONS
from pgd_splicer.sidechain import bond_angles, bond_lengths


def NO_VALUE(field):
    """
    Helper function for determining the value to use when the field is an
    invalid value.
    """

    if field in ('bm', 'bg', 'bs'):
        return 0
    else:
        return None

AA3to1 =  {
    'ALA' : 'a',
    'ARG' : 'r',
    'ASN' : 'n',
    'ASP' : 'd',
    'CYS' : 'c',
    'GLU' : 'e',
    'GLN' : 'q',
    'GLY' : 'g',
    'HIS' : 'h',
    'ILE' : 'i',
    'LEU' : 'l',
    'LYS' : 'k',
    'MET' : 'm',
    'PHE' : 'f',
    'PRO' : 'p',
    'SER' : 's',
    'THR' : 't',
    'TRP' : 'w',
    'TYR' : 'y',
    'VAL' : 'v',
}


aa_class = {
    'r':(Sidechain_ARG,'sidechain_ARG'),
    'n':(Sidechain_ASN,'sidechain_ASN'),
    'd':(Sidechain_ASP,'sidechain_ASP'),
    'c':(Sidechain_CYS,'sidechain_CYS'),
    'q':(Sidechain_GLN,'sidechain_GLN'),
    'e':(Sidechain_GLU,'sidechain_GLU'),
    'h':(Sidechain_HIS,'sidechain_HIS'),
    'i':(Sidechain_ILE,'sidechain_ILE'),
    'l':(Sidechain_LEU,'sidechain_LEU'),
    'k':(Sidechain_LYS,'sidechain_LYS'),
    'm':(Sidechain_MET,'sidechain_MET'),
    'f':(Sidechain_PHE,'sidechain_PHE'),
    'p':(Sidechain_PRO,'sidechain_PRO'),
    's':(Sidechain_SER,'sidechain_SER'),
    't':(Sidechain_THR,'sidechain_THR'),
    'w':(Sidechain_TRP,'sidechain_TRP'),
    'y':(Sidechain_TYR,'sidechain_TYR'),
    'v':(Sidechain_VAL,'sidechain_VAL')
}


class InvalidResidueException(Exception):
    """
    Exception identifying something wrong while processing
    a protein residue.  this is used to jump to a common error
    handling routine.
    """
    pass


class ProcessPDBTask():
    """
    Task that takes a list of pdbs and processes the files extracting
    geometry data from the files.  The data is stored in Protein, Chain and
    Residue models and commited to the database.
    """
    total_proteins = 0
    finished_proteins = 0

    def progress(self):
        if not self.total_proteins:
            return 0
        return self.finished_proteins / self.total_proteins * 100


    def work(self, pdbs):
        """
        Work function - expects a list of pdb file prefixes.
        """
        # process a single protein dict, or a list of proteins

        if not isinstance(pdbs, list):
            pdbs = [pdbs]

        logging.info("processing {} proteins".format(len(pdbs)))

        self.total_proteins = len(pdbs)
        skipped = 0
        imported = 0

        pool = Pool(maxtasksperchild=50)
        started = datetime.now()

        results = pool.imap_unordered(workhorse, pdbs)
        for result in results:
            if result:
                imported += 1
            else:
                skipped += 1
            self.finished_proteins += 1

            percent = self.progress()
            now = datetime.now()
            elapsed = now - started
            if int(percent):
                remaining = elapsed * 100 // int(percent)
                remaining -= elapsed
            else:
                remaining = elapsed * 100
            logging.info("-" * 42)
            logging.info("processed protein {} of {}".format(self.finished_proteins, self.total_proteins))
            logging.info("{} imported, {} skipped ({:f}%)".format(imported, skipped, percent))
            logging.info("{} elapsed, {} remaining".format(elapsed, remaining))
            logging.info("-" * 42)

        logging.info('processing complete: {} elapsed, {} imported, {} skipped'.format(elapsed, imported, skipped))

        # return only the code of proteins inserted or updated
        # we no longer need to pass any data as it is contained in the database
        # for now assume everything was updated
        return pdbs


def workhorse(line):
    """
    Update a single protein entry.

    This is a viable entrypoint for parallelizing the process.
    """

    # Break out selection line data here.
    values = line.split(' ')

    data = {
        'code': values[0],
        'chains': [c for c in values[1]],
        'threshold': float(values[2]),
        'resolution': float(values[3]),
        'rfactor': float(values[4]),
        'rfree': float(values[5])
    }

    logger = logging.getLogger(values[0])
    pdb = {
        'data': data,
        'logger': logger
        }

    # only update pdbs if they are newer
    if pdb_file_is_newer(pdb):
        return process_pdb(pdb)
    else:
        logger.info('pdb file is not newer')
        return False


@transaction.commit_manually
def process_pdb(pdb):
    """
    Process an individual pdb file
    """

    # create a copy of the data.  This dict will have a large amount of data
    # added to it as the protein is processed.  This prevents memory leaks
    # due to the original dict having a reference held outside this method.
    # e.g. if it were looped over with a large list of PDBs
    # data = data.copy()
    data = pdb['data']
    logger = pdb['logger']

    try:
        residue_props = None
        code = data['code']
        chains_filter = data.get("chains")
        logger.debug('process_pdb')

        # update datastructure
        data['chains'] = {}

        # 1) parse with bioPython
        data = parseWithBioPython(code, data, chains_filter)

        # 2) Create/Get Protein and save values
        try:
            protein = ProteinModel.objects.get(code=code)
            logger.debug('protein exists')
        except ProteinModel.DoesNotExist:
            logger.debug('protein does not exist')
            protein = ProteinModel()
        protein.code       = code
        protein.threshold  = float(data['threshold'])
        protein.resolution = float(data['resolution'])
        protein.rfactor    = float(data['rfactor'])
        protein.rfree      = float(data['rfree'])
        protein.pdb_date   = data['pdb_date']
        protein.save()

        # 3) Get/Create Chains and save values
        chains = {}
        for chaincode, residues in data['chains'].items():
            chainId = '%s%s' % (protein.code, chaincode)
            try:
                chain = protein.chains.get(id=chainId)
                logger.debug('{}: chain exists'.format(chaincode))
            except ChainModel.DoesNotExist:
                logger.debug('{}: chain does not exist'.format(chaincode))
                chain = ChainModel()
                chain.id      = chainId
                chain.protein = protein
                chain.code    = chaincode
                chain.save()

                protein.chains.add(chain)
            #create dictionary of chains for quick access
            chains[chaincode] = chain

            # 4) iterate through residue data creating residues
            for residue_props in sorted(residues.values(), key=itemgetter("chainIndex")):

                # 4a) find the residue object so it can be updated or create a new one
                try:
                    residue = chain.residues.get(oldID=str(residue_props['oldID']))
                    logger.debug('{}: residue exists'.format(residue))
                except ResidueModel.DoesNotExist:
                    #not found, create new residue
                    residue = ResidueModel()
                    residue.protein = protein
                    residue.chain   = chain
                    residue.chainID = chain.id[4]
                    # logger.debug('{}: residue does not exist'.format(residue))

                # 4b) copy properties into a residue object
                #     property keys should match property name in object
                residue.__dict__.update(residue_props)

                # 4c) set previous
                if "prev" in residue_props:
                    residue.prev = old_residue

                # 4d) find and create sidechain if needed.  set the property
                #     in the residue for the correct sidechain type
                if 'sidechain' in residue_props:
                    klass, name = aa_class[residue_props['aa']]
                    try:
                        sidechain = getattr(residue, name)
                        if not sidechain:
                            sidechain = klass()
                    except:
                        sidechain = klass()
                    sidechain.__dict__.update(residue_props['sidechain'])
                    sidechain.save()
                    residue.__setattr__(name, sidechain)

                # 4e) save
                residue.save()
                chain.residues.add(residue)

                # 4f) Update old_residue.next
                if "prev" in residue_props:
                    old_residue.next = residue
                    old_residue.save()


                old_residue = residue
            logger.debug('{}: {} residues'.format(chainId, len(residues)))


    except Exception, e:
        import traceback
        exceptionType, exceptionValue, exceptionTraceback = sys.exc_info()
        from cStringIO import StringIO
        tb = StringIO()
        traceback.print_tb(exceptionTraceback, limit=10, file=tb)
        logger.error('{}: EXCEPTION in Residue: {} {}'.format(code, e.__class__, e))
        logger.error('{}: PROPS: {}'.format(code, residue_props))
        logger.error('{}: TRACEBACK: {}'.format(code, tb.getvalue()))
        tb.close()
        transaction.rollback()
        return False

    # 5) entire protein has been processed, commit transaction
    transaction.commit()
    return True


def pdb_file_is_newer(pdb):
    """
    Compares if the pdb file used as an input is newer than data already
    in the database.

    This is used to prevent processing proteins if they do not need to be
    processed.
    """

    data = pdb['data']
    logger = pdb['logger']
    code = data['code']
    lfile = localfile(code)
    if os.path.exists(lfile):
        pdb_date = timezone.make_aware(datetime.fromtimestamp(os.path.getmtime(lfile)), timezone.get_default_timezone())

    else:
        logger.error('File not found: {}'.format(lfile))
        return False
    try:
        protein = ProteinModel.objects.get(code=code)
    except ProteinModel.DoesNotExist:
        # Protein not in database, pdb is new
        data['pdb_date'] = pdb_date
        return True

    data['pdb_date'] = pdb_date

    return protein.pdb_date < pdb_date


def amino_present(s):
    """
    Whether a given line contains a valid amino acid.

    NB: "SEC" is a valid amino acid, regardless of its absence from
    AA3to1.

    See #17223 for details.
    """

    return any(amino in s for amino in AA3to1)

def hetatm_amino(s):
    """
    Whether a given line refers to a HETATM amino acid.

    These lines are undesirable and we should discard them.
    """

    return s.startswith("HETATM ") and amino_present(s)


def atom_noamino(s):
    """
    Whether a given line is an ATOM line with no valid amino acid.

    These lines are undesirable and we should discard them.
    """

    return s.startswith("ATOM ") and not amino_present(s)


# class ChainSelect(Select):
#     """
#     This class selects only the chains specified in the selection file.
#     """
# 
#     def __init__(self, chains_filter):
#         self.chains_filter = chains_filter
# 
#     def accept_chain(self, chain):
#         return (chain.get_id() in self.chains_filter)


class PGDSelect(Select):
    """
    This class filters incoming PDB files to include:
     - only the chains listed in the selection files
     - only the atoms with the greatest average occupancy per residue
    and to exclude:
     - all residues with residue names that are not amino acids
     - all residues with non-blank hetflags
    """

    def __init__(self, chains_filter):
        self.best_atoms = {}
        self.chains_filter = chains_filter

    def accept_chain(self, chain):
        """
        All chains which are not included in chains_filter are removed.
        """
        return (chain.get_id() in self.chains_filter)

    def accept_residue(self, residue):
        """
        This method is being overloaded to perform the calculations required
        to properly accept atoms in the model based on the highest average
        occupancy value per altloc.
        """

        # Only residues that are amino acids are allowed.
        resname = residue.get_resname()
        if resname not in AA3to1:
            return False

        # No non-blank hetflags are allowed.
        hetflag, resseq, icode = residue.get_id()
        if hetflag != ' ':
            return False

        # Calculate the highest average occupancy value per altloc.
        if residue.is_disordered():
            res_occ = {}
            for atom in residue:
                name = atom.name
                res_occ[name] = {}
                if atom.is_disordered():
                    altloc = atom.get_altloc()
                    occ = atom.get_occupancy()
                    try:
                        res_occ[name][altloc].append(occ)
                    except KeyError:
                        res_occ[name][altloc] = [occ]

            occ_list = {}
            for occ_vals in [res_occ[k] for k in res_occ if res_occ[k]]:
                for altloc in occ_vals:
                    for occ in occ_vals[altloc]:
                        try:
                            occ_list[altloc].append(occ)
                        except KeyError:
                            occ_list[altloc] = [occ]

            best_altlocs = sorted(occ_list, key=lambda k: sum(occ_list[k])/len(occ_list[k]), reverse=True)
            self.best_atoms[residue] = {atom: next(altloc for altloc in best_altlocs if altloc in res_occ[atom] or ' ') for atom in res_occ}

        return True

    def accept_atom(self, atom):
        """
        This method uses the information generated by accept_residue and
        stored in best_altlocs to accept only those atoms which are not
        disordered or with an altloc that matches the best altloc for that
        particular set of atoms.
        """

        if atom.is_disordered():
            return self.best_atoms[atom.get_parent()][atom.name] == atom.get_altloc()
        else:
            return True


def parseWithBioPython(code, props, chains_filter=None):
    """
    Parse values from file that can be parsed using BioPython library
    @return a dict containing the properties that were processed
    """

    logger = logging.getLogger(code)

    lfile = localfile(code)

    chains = props['chains']

    decompressed = NamedTemporaryFile()
    gunzipped = GzipFile(lfile)

    # Go through, one line at a time, and discard lines that have the bad
    # HETATM pattern. This is largely for 2VQ1, see #8319 for details.
    # Also remove any ATOM lines with invalid amino acids.
    # See #17223 for details.
    for line in gunzipped:
        if not hetatm_amino(line) and not atom_noamino(line):
            decompressed.write(line)

    # Be kind; rewind.
    decompressed.seek(0)

    # Remove unused chains from PDB file.
    # if chains_filter:
    #     logger.warning('yay chains_filter {}'.format(chains_filter))
    #     chains_structure = Bio.PDB.PDBParser().get_structure('pdbname', decompressed.name)
    #     chains_io = PDBIO()
    #     chains_io.set_structure(chains_structure)
    #     chains_io.save(decompressed.name, select=ChainSelect(chains_filter))

    # Open structure for pre-cleaned PDB file.
    select_structure = Bio.PDB.PDBParser().get_structure('pdbname',
                                                         decompressed.name)

    # write new PDB based on conformation changes
    io = PDBIO()
    io.set_structure(select_structure)
    io.save(decompressed.name, select=PGDSelect(chains_filter))

    # write PDB to current directory as well
    # import shutil
    # shutil.copy(decompressed.name, '{}-postocc.ent'.format(code))

    # Reopen structure from cleaned PDB file.
    structure = Bio.PDB.PDBParser().get_structure('pdbname',
                                                  decompressed.name)

    # dssp can't do multiple models. if we ever need to, we'll have to
    # iterate through them
    dssp = Bio.PDB.DSSP(model=structure[0], pdb_file=decompressed.name,
                        dssp='dsspcmbi')

    if not dssp.keys():
        raise Exception("No chains were parsed!")

    for chain in structure[0]:
        chain_id = chain.get_id()

        # only process selected chains
        # if chains_filter and not chain_id in chains_filter:
        #     logger.debug('skipping chain {}'.format(chain_id))
        #     continue

        # construct structure for saving chain
        if not chain_id in props['chains']:
            residues = {}
            props['chains'][chain_id] = residues
            logger.debug('processing chain {} ({} residues)'.format(chain, len(chain)))

        newID = 0

        #iterate residues
        res_old_id = None
        oldN       = None
        oldCA      = None
        oldC       = None
        prev       = None

        for res in chain:
            try:
                newID += 1
                terminal = False
                hetflag, res_id, icode = res.get_id()

                # We can't handle any hetflags. This is primarily to filter
                # out water, but there can be others as well.
                # if hetflag != ' ':
                #     raise InvalidResidueException("HetCode %r" % hetflag)

                resname = res.resname

                # We can't deal with residues that aren't of amino acids.
                # if resname not in AA3to1:
                #     raise InvalidResidueException("Bad amino acid %r" %
                #                                   resname)

                # XXX Get the dictionary of atoms in the Main conformation.
                # BioPython should do this automatically, but it does not
                # always choose the main conformation.  Leading to some
                # Interesting results
                # Occupancy selection code renders altloc choice unnecessary.
                atoms = {}
                for atom in res.get_unpacked_list():
                    atoms[atom.name] = atom

                # Exclude water residues
                # Exclude any Residues that are missing _ANY_ of the
                #     mainchain atoms.  Any atom could be missing
                all_mainchain = ('N' in atoms) and ('CA' in atoms) and ('C' in atoms) and ('O' in atoms)
                if not all_mainchain:
                    raise InvalidResidueException("Missing atom")

                # Create dictionary structure and initialize all values.  All
                # Values are required.  Values that are not filled in will retain
                # the NO_VALUE value.
                #
                # Store residue properties using OLD_ID as the key to ensure it is
                # unique.  We're including residues from all chains in the same
                # dictionary and chainindex may have duplicates.
                old_id = res_id if icode == ' ' else '%s%s' % (res_id, icode)
                if old_id in residues:
                    res_dict = residues[old_id]
                else:
                    # This residue doesn't exist yet.
                    res_dict = {}
                    residues[old_id] = res_dict
                    res_dict['oldID'] = old_id

                initialize_all_geometry(res_dict)

                # Get Properties from DSSP and other per residue properties
                chain = res.get_parent().get_id()
                key = chain, (hetflag, res_id, icode)
                if key in dssp:
                    (secondary_structure, accessibility,
                     relative_accessibility, phi, psi) = dssp[key][2:7]
                else:
                    raise InvalidResidueException("Key %r not in DSSP" %
                                                  (key,))

                res_dict['chain'] = chain
                res_dict['ss'] = secondary_structure
                res_dict['aa'] = AA3to1[resname]
                res_dict['h_bond_energy'] = 0.00

                # Get Vectors for mainchain atoms and calculate geometric angles,
                # dihedral angles, and lengths between them.
                N    = atoms['N'].get_vector()
                CA   = atoms['CA'].get_vector()
                C    = atoms['C'].get_vector()
                CB   = atoms['CB'].get_vector() if "CB" in atoms else None
                O    = atoms['O'].get_vector()

                if oldC:
                    # determine if there are missing residues by calculating
                    # the distance between the current residue and previous
                    # residue.  If the L1 distance is greater than 2.5 it
                    # cannot possibly be the correct order of residues.
                    L1 = calc_distance(oldC,N)

                    if L1 < 2.5:
                        # properties that span residues
                        residues[res_old_id]['a6'] = calc_angle(oldCA,oldC,N)
                        residues[res_old_id]['a7'] = calc_angle(oldO,oldC,N)
                        residues[res_old_id]['psi'] = calc_dihedral(oldN,oldCA,oldC,N)
                        residues[res_old_id]['ome'] = calc_dihedral(oldCA,oldC,N,CA)
                        residues[res_old_id]['next'] = newID
                        res_dict['prev'] = residues[res_old_id]['chainIndex']
                        res_dict['a1']     = calc_angle(oldC,N,CA)
                        res_dict['phi']    = calc_dihedral(oldC,N,CA,C)
                        res_dict['L1'] = L1

                        # proline has omega-p property
                        if prev and resname == 'PRO' and 'CD' in atoms:
                            CD = atoms['CD'].get_vector()
                            res_dict['omep'] = calc_dihedral(oldCA, oldC, N, CD)

                        terminal = False

                if terminal:
                    # break in the chain,
                    # 1) add terminal flags to both ends of the break so
                    #    the break can quickly be found.
                    # 2) skip a number in the new style index.  This allows
                    #    the break to be visible without checking the
                    #    terminal flag
                    newID += 1
                    res_dict['terminal_flag'] = True
                    residues[res_old_id]['terminal_flag'] = True

                # newID cannot be set until after we determine if it is terminal
                res_dict['chainIndex'] = newID

                res_dict['L2'] = calc_distance(N,CA)
                res_dict['L4'] = calc_distance(CA,C)
                res_dict['L5'] = calc_distance(C,O)
                res_dict['a3'] = calc_angle(N,CA,C)
                res_dict['a5'] = calc_angle(CA,C,O)

                if CB:
                    res_dict['a2'] = calc_angle(N,CA,CB)
                    res_dict['a4'] = calc_angle(CB,CA,C)
                    res_dict['L3'] = calc_distance(CA,CB)
                    res_dict['zeta'] = calc_dihedral(CA, N, C, CB)

                # Calculate Bg - bfactor of the 4th atom in Chi1.
                try:
                    atom_name = CHI_MAP[resname][0][3]
                    res_dict['bg'] = res[atom_name].get_bfactor()
                except KeyError:
                    # not all residues have chi
                    pass


                # Other B Averages
                #    Bm - Average of bfactors in main chain.
                #    Bm - Average of bfactors in side chain.
                main_chain = []
                side_chain = []
                for name in atoms:
                    if name in ('N', 'CA', 'C', 'O','OXT'):
                        main_chain.append(atoms[name].get_bfactor())
                    elif name in ('H'):
                        continue
                    else:
                        side_chain.append(atoms[name].get_bfactor())

                if main_chain != []:
                    res_dict['bm'] = sum(main_chain) // len(main_chain)

                if side_chain != []:
                    res_dict['bs'] = sum(side_chain) // len(side_chain)

                #    occscs - Min of side chain
                #    occm   - Min of main chain
                #    occ_m  - List containing main chain occupancy values
                #    occ_scs- List containing side chain occupancy values
                #    issue link - https://code.osuosl.org/issues/17565
                occ_m, occ_scs = [], []
                for name in atoms:
                    if name in ('N', 'CA', 'C', 'O','OXT', 'CB'):
                        occ_m.append(atoms[name].get_occupancy())
                    elif name in ('H'):
                        continue
                    else:
                        occ_scs.append(atoms[name].get_occupancy())

                if occ_m != []:
                    res_dict['occm'] = min(occ_m)

                if occ_scs != []:
                    res_dict['occscs'] = min(occ_scs)

                if res_dict['aa'] == 'g' or res_dict['aa'] == 'a' :
                    res_dict['occscs'] = 1.0

                # CHI corrections - Some atoms have symettrical values
                # in the sidechain that aren't guarunteed to be listed
                # in the correct order by within the PDB file.  The correct order can be
                # determined by checking for the larger of two Chi values.  If the 2nd value
                # listed in CHI_CORRECTIONS_TESTS is larger, then corrections are needed. There
                # may be multiple corrections per Residue which are listed in CHI_CORRECTIONS.
                # each pair of atoms will be swapped so that any future calculation will use
                # the correct values.
                #
                # Both angles are required to determine whether atoms are labeled correctly.
                # If only one angle is present then we check to see if it is less than 90
                # degrees.  if it is less than 90 degress then it is considered ChiN-1
                # See: Ticket #1545 for more details.
                if resname in CHI_CORRECTIONS_TESTS:
                    values = calc_chi(atoms, prev, CHI_CORRECTIONS_TESTS[resname])
                    correct = False
                    if not len(values):
                        for atom1, atom2 in CHI_CORRECTIONS[resname]:
                            if atom1 in atoms: del atoms[atom1]
                            if atom2 in atoms: del atoms[atom2]

                    elif len(values) == 1:
                        if 'chi1' in values and abs(values['chi1']) > 90:
                            correct = True
                        elif 'chi2' in values and abs(values['chi2']) < 90:
                            correct = True

                    elif abs(values['chi2']) < abs(values['chi1']):
                        correct = True

                    if correct:
                        for atom1, atom2 in CHI_CORRECTIONS[resname]:
                            old_atom1 = atoms.get(atom1)
                            old_atom2 = atoms.get(atom2)
                            if old_atom1:
                                atoms[atom2] = old_atom1
                            else:
                                # del atoms[atom2]
                                atoms.pop(atom2, None)
                            if old_atom2:
                                atoms[atom1] = old_atom2
                            else:
                                # del atoms[atom1]
                                atoms.pop(atom1, None)


                #Calculate CHI values.  The mappings for per peptide chi's are stored
                #in a separate file and a function is used to calculate the chi based
                #based on the peptide of this residue and the lists of atoms in the
                #chi mappings.
                if resname in CHI_MAP:
                    res_dict.update(calc_chi(atoms, prev, CHI_MAP[resname]))
                sidechain = {}
                if resname in bond_lengths:
                    calc_sidechain_lengths(atoms, sidechain, bond_lengths[resname])
                if resname in bond_angles:
                    calc_sidechain_angles(atoms, prev, sidechain, bond_angles[resname])
                if sidechain:
                    res_dict['sidechain'] = sidechain

                # Reset for next pass.  We save some relationships which span two atoms.
                res_old_id = old_id
                oldN       = N
                oldCA      = CA
                oldC       = C
                oldO       = O
                prev = res

            except InvalidResidueException, e:
                # something has gone wrong in the current residue
                # indicating that it should be excluded from processing
                # log a warning
                logger.warning('invalid residue! chain {} residue: {} exception: {}'.format(chain_id, res.get_id(), e))
                if oldC:
                    residues[res_old_id]['terminal_flag'] = True
                    newID += 1
                oldN       = None
                oldCA      = None
                oldC       = None
                prev       = None
                if res_id in residues:
                    del residues[res_id]
                # Disordered residues are stored under old_id, not res_id.
                # First check to see if old_id is defined.
                if 'old_id' in vars():
                    if old_id in residues:
                        del residues[old_id]

        logger.debug('processed {} residues'.format(len(residues)))

    return props


def initialize_geometry(residue, geometry_list):
    """
    Initialize the dictionary for geometry data
    """

    for item in geometry_list:
        if item not in residue:
            residue[item] = NO_VALUE(item)


def initialize_all_geometry(residue):
    """
    Set up a residue dictionary with all the appropriate keys for geometry
    data.
    """

    length_list = ['L1', 'L2', 'L3', 'L4', 'L5', 'L6', 'L7', 'bg', 'bs', 'bm']
    angles_list = ['a1', 'a2', 'a3', 'a4', 'a5', 'a6', 'a7']
    dihedral_list = ['psi', 'ome', 'phi', 'zeta', 'chi1', 'chi2', 'chi3', 'chi4']
    initialize_geometry(residue, length_list)
    initialize_geometry(residue, angles_list)
    initialize_geometry(residue, dihedral_list)


def calc_distance(atom1, atom2):
    """
    Calculates distance between atoms because this is not built into BIOPython

    scribed from http://www.scribd.com/doc/9816032/BioPython-for-Bioinfo
    """
    dx = atom1[0] - atom2[0]
    dy = atom1[1] - atom2[1]
    dz = atom1[2] - atom2[2]
    return sqrt(dx*dx + dy*dy + dz*dz)


def calc_angle(atom1, atom2, atom3):
    """
    overridding pdb version of function to return values converted to degress
    """
    return degrees(pdb_calc_angle(atom1, atom2, atom3))


def calc_dihedral(atom1, atom2, atom3, atom4):
    """
    overridding pdb version of function to return values converted to degress
    """
    return degrees(pdb_calc_dihedral(atom1, atom2, atom3, atom4))


def calc_chi(residue, residue_prev, mapping):
    """
    Calculates Values for CHI using the predefined list of CHI angles in
    the CHI_MAP.  CHI_MAP contains the list of all peptides and the atoms
    that make up their different chi values.  This function will process
    the values known to exist, it will also skip chi values if some of the
    atoms are missing

    @param mapping: list of mappings to use
    """
    residue_dict = {}
    try:
        for i in range(len(mapping)):
            chi_atom_names= mapping[i]
            try:
                chi_atoms = []
                for n in chi_atom_names:
                    if residue_prev and n[-2:]=='-1':
                        chi_atoms.append(residue_prev[n[:-2]].get_vector())
                    else:
                        chi_atoms.append(residue[n].get_vector())
                chi = calc_dihedral(*chi_atoms)
                residue_dict['chi%i'%(i+1)] = chi
            except KeyError:
                #missing an atom
                continue

    except KeyError:
        # this residue type does not have chi
        pass
    return residue_dict


def calc_sidechain_lengths(residue, residue_dict, mapping):
    """
    Calculates Values for sidechain bond lengths. Uses a predefined list
    from sidechain.py, specifically bond_lengths.
    """
    try:
        for i in range(len(mapping)):
            atom_names = mapping[i]
            try:
                sidechain_atoms = [residue[n].get_vector() for n in atom_names]
                sidechain_length = calc_distance(*sidechain_atoms)
                residue_dict['%s_%s' % (atom_names[0], atom_names[1])] = sidechain_length
            except KeyError:
                #missing an atom
                continue

    except KeyError:
        # this residue type does not have sidechain lengths
        pass


def calc_sidechain_angles(residue, residue_prev, residue_dict, mapping):
    """
    Calculates Values for sidechain bond angles. Uses a predefined list
    from sidechain.py, specifically bond_angles.
    """
    try:
        for i in range(len(mapping)):
            atom_names = mapping[i]
            try:
                sidechain_atoms = []
                for n in atom_names:
                    if residue_prev and n[-2:]=='-1':
                        sidechain_atoms.append(residue_prev[n[:-2]].get_vector())
                    else:
                        sidechain_atoms.append(residue[n].get_vector())
                sidechain_angle = calc_angle(*sidechain_atoms)
                angle_key = '%s_%s_%s' % (atom_names[0], atom_names[1], atom_names[2])
                angle_key = angle_key.replace('-','_')
                residue_dict[angle_key] = sidechain_angle
            except KeyError:
                #missing an atom
                continue

    except KeyError:
        # this residue type does not have sidechain angles
        pass


if __name__ == '__main__':
    """
    Run if file is executed from the command line
    """

    task = ProcessPDBTask()

    # Configure logging.
    logging.basicConfig(level=logging.DEBUG,
                        format='%(asctime)s %(name)-4s %(levelname)-8s %(message)s',
                        datefmt='%Y-%m-%d %H:%M:%S',
                        filename='ProcessPDB.log',
                        filemode='w')

    # Log all info messages to the console.
    console = logging.StreamHandler()
    console.setLevel(logging.INFO)
    formatter = logging.Formatter('%(name)-4s %(levelname)-8s %(message)s')
    console.setFormatter(formatter)
    logging.getLogger('').addHandler(console)

    pdbs = []

    argv = sys.argv
    if len(argv) == 1:
        print 'Usage:'
        print '   ProcessPDBTask code chains threshold resolution rfactor rfree [repeat]'
        print '       chains are a string of chain ids: ABCXYZ'
        print ''
        print '   <cmd> | ProcessPDBTask --pipein'
        print '   piped protein values must be separated by newlines'
        sys.exit(0)

    elif len(argv) == 2 and argv[1] == '--pipein':
        for line in sys.stdin:
            pdbs.append(line)

    else:
        for i in range(1,len(argv),6):
            try:
                line = " ".join(argv[i:i+6])
                pdbs.append(line)
            except IndexError, e:
                print e
                print 'Usage: ProcessPDBTask.py code chain threshold resolution rfactor rfree...'
                sys.exit(0)

    task.work(pdbs)
