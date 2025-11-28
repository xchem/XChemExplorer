import coot


# return a list of molecule numbers (closed and open)
# The elements of the returned list need to be tested against
# is_valid_model_molecule_qm
#
def molecule_number_list():
    ret = []
    for mol_no in range(coot.graphics_n_molecules()):
        if valid_map_molecule_qm(mol_no) or valid_model_molecule_qm(mol_no):
            ret.append(mol_no)
    return ret


# The screen centre.
#
# return the rotation centre as a 3 membered list of numbers
# is python list [...] !!!
#


def rotation_centre():
    return [
        coot.rotation_centre_position(0),
        coot.rotation_centre_position(1),
        coot.rotation_centre_position(2),
    ]


# Return the molecule centre as a list of 3 numbers.
#
#  Note: mol_cen could contain values less than -9999.
#


def molecule_centre(imol):
    return [
        coot.molecule_centre_internal(imol, 0),
        coot.molecule_centre_internal(imol, 1),
        coot.molecule_centre_internal(imol, 2),
    ]


# Move the centre of molecule number imol to the current screen centre
#


def move_molecule_to_screen_centre(imol):
    if valid_model_molecule_qm(imol):
        rotate_centre = rotation_centre()
        coot.translate_molecule_by(
            imol,
            (rotate_centre[0] - molecule_centre(imol)[0]),
            (rotate_centre[1] - molecule_centre(imol)[1]),
            (rotate_centre[2] - molecule_centre(imol)[2]),
        )


# This is a short name for the above.
# deftexi move_molecule_here
move_molecule_here = move_molecule_to_screen_centre


# return a list of chain ids for given molecule number @var{imol}.
# return empty list on error
#


def chain_ids(imol):
    chain_id_is = []
    number_of_chains = coot.n_chains(imol)
    for chain_no in range(number_of_chains):
        chain_id_is.append(coot.chain_id_py(imol, chain_no))
    return chain_id_is


# python (schemeyish) interface to eponymous scripting interface function!?
# return True or False
#


def valid_model_molecule_qm(imol):
    if coot.is_valid_model_molecule(imol) == 1:
        return True
    else:
        return False


# python (schemeyish) interface to eponymous scripting interface function.
# return True or False
#


def valid_map_molecule_qm(imol):
    if coot.is_valid_map_molecule(imol) == 1:
        return True
    else:
        return False
