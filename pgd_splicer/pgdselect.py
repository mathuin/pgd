import logging
from Bio.PDB.PDBIO import Select


class PGDSelect(Select):
    """
    This class performs three tasks:
     * it removes atoms and residues which do not work with this software, and
     * it selects the best atoms based on highest average occupancy.
    """

    def __init__(self, chains_filter=None, logger=None):
        self.best_atoms = {}
        self.chains_filter = chains_filter
        if logger:
            self.logger = logger
        else:
            self.logger = logging.getLogger('').addHandler(logging.NullHandler())
        self.not_in_AA3to1 = []
        self.has_hetflag = []
        self.missing_atom = []
        self.no_acceptable_altlocs = []
        self.not_disordered = []
        self.logger.debug("finished init")

    def __del__(self):
        if self.not_disordered:
            self.logger.info("residues not disordered: {}".format(len(self.not_disordered)))
        if self.has_hetflag:
            self.logger.info("residues with hetflags: {}".format(len(self.has_hetflag)))
        if self.missing_atom:
            self.logger.info("residues missing atoms: {}".format(len(self.missing_atom)))
        if self.no_acceptable_altlocs:
            self.logger.info("residues with no acceptable altlocs: {}".format(len(self.no_acceptable_altlocs)))
        if self.not_in_AA3to1:
            self.logger.info("residues not in AA3to1: {}".format(len(self.not_in_AA3to1)))

    def has_mainchain_atoms(self, residue):
        atoms = {atom.name: atom for atom in residue.get_unpacked_list()}
        return (('N' in atoms) and ('CA' in atoms) and ('C' in atoms) and ('O' in atoms))

    def accept_chain(self, chain):
        """
        Each disordered residue contains one or more disordered atom structures.
        These structures contain one or more atoms, each with its own altloc and
        occupancy values.  As each atom associated with a particular altloc has
        its own occupancy value, the altloc with the highest average occupancy
        value is identified.  Atoms with other altloc values are removed from
        the structure.

        Because the occupancy value is calculated across all residues with a
        given sequence number, the calculation must occur in the chain
        acceptance method.
        """

        best_atoms = {}
        best_altlocs = {}

        # only process selected chains
        # XXX: disabled
        if self.chains_filter and not chain.get_id() in self.chains_filter:
            self.logger.debug("chain {} not in chains_filter".format(chain.get_id()))
            return False

        for residue in chain.get_unpacked_list():

            # Only process residues without hetflags.
            hetflag, resseq, icode = residue.get_id()
            if hetflag != ' ':
                self.logger.debug("residue {} has hetflag".format(residue))
                self.has_hetflag.append(residue)
                continue

            # Only process residues in AA3to1.
            # resname = residue.resname
            # if resname not in AA3to1:
            #     self.logger.debug("residue {} not in AA3to1".format(residue))
            #     self.not_in_AA3to1.append(residue)
            #     continue

            # Only process residues with all atoms.
            if not self.has_mainchain_atoms(residue):
                self.logger.debug("residue {} missing atom".format(residue))
                self.missing_atom.append(residue)
                continue

            # Calculate occupancy for disordered residues.
            res_occ = {}
            atoms = []
            if residue.is_disordered():
                self.logger.debug("residue: {} is disordered".format(residue))
                for atom in residue.get_unpacked_list():
                    if atom.element != 'H':
                        name = atom.name
                        altloc = atom.get_altloc()
                        atoms.append('{}|{}'.format(name, altloc))
                        if atom.is_disordered():
                            try:
                                res_occ[name]
                            except KeyError:
                                res_occ[name] = {}
                            occ = atom.get_occupancy()
                            res_occ[name].update({altloc: occ})
                self.logger.debug("atom list: {}".format(res_occ))
                self.logger.debug("atoms: {}".format(atoms))

            # 4P3H A/154 has only one disordered atom: hydrogen.
            # Thus it is considered ordered for our purpose.
            if res_occ == {}:
                self.logger.debug("residue {} has no relevant disordered atoms".format(residue))
                self.best_atoms[residue] = {}
                self.not_disordered.append(residue)
                continue
            else:
                occ_list = {}
                for occ_vals in [res_occ[k] for k in res_occ if res_occ[k]]:
                    for altloc in occ_vals:
                        try:
                            occ_list[altloc].append(occ_vals[altloc])
                        except KeyError:
                            occ_list[altloc] = [occ_vals[altloc]]
                self.logger.debug("altloc list: {}".format(occ_list))

                # What altloc has the highest average occupancy in this residue?
                for best_altloc in sorted(occ_list, key=lambda k: sum(occ_list[k])/len(occ_list[k]), reverse=True):
                    # Does this altloc represent a complete amino acid?
                    self.logger.debug("best_altloc: {}".format(best_altloc))
                    occ_altloc = sum(occ_list[best_altloc])/len(occ_list[best_altloc])
                    self.logger.debug("occ_altloc: {}".format(occ_altloc))
                    best_atoms[residue] = {atom: best_altloc if best_altloc in res_occ[atom] else ' ' for atom in res_occ}
                    self.logger.debug("best_atoms: {}".format(best_atoms[residue]))
                    for atom, altloc in best_atoms[residue].iteritems():
                        key = '{}|{}'.format(atom, altloc)
                        if key not in atoms:
                            self.logger.debug("altloc {} missing atoms".format(best_altloc))
                            break
                    else:
                        # Store the best altloc for all residues sharing this sequence number.
                        try:
                            best_altlocs[resseq].update({residue: occ_altloc})
                        except KeyError:
                            best_altlocs[resseq] = {residue: occ_altloc}
                        self.logger.debug("best_altlocs[{}] = {}".format(resseq, best_altlocs[resseq]))
                        break
                else:
                    self.logger.debug("residue {} has no acceptable altlocs".format(residue))
                    self.no_acceptable_altlocs.append(residue)

        # The residue with the best altloc has the best atoms.
        for resseq, altlocs in best_altlocs.iteritems():
            altlocs = best_altlocs[resseq]
            best_residue = max(altlocs, key=altlocs.get)
            self.logger.debug("best residue for {} is {}".format(resseq, best_residue))
            self.best_atoms[best_residue] = best_atoms[best_residue]
            self.logger.debug("best atoms for {} are {}".format(resseq, best_atoms[best_residue]))

        return True

    def accept_residue(self, residue):
        """
        This method uses the information generated by accept_chain and stored
        in best_atoms to determine which disordered residue is selected to
        represent that residue in the chain.
        """

        return residue in self.best_atoms

    def accept_atom(self, atom):
        """
        This method uses the information generated by accept_chain and stored
        in best_atoms to determine which disordered atom is selected to
        represent that atom in the residue.
        """

        if atom.element != 'H':
            if atom.is_disordered():
                if self.best_atoms[atom.get_parent()][atom.name] == atom.get_altloc():
                    atom.set_altloc(' ')
                    return True
                else:
                    return False
            else:
                return True
        else:
            return False
