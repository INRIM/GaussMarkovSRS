import numpy as np
import numpy.lib.recfunctions as rfn

from gauss_markov_srs.cross_index import get_cross_index


# get a unique key for each measurement
def key(entry):
    s = "_".join((f'{entry["Id"]}', entry["Ref"], entry["Atom1"], entry["Atom2"], entry["Sup"]))
    return s.strip("_")  # get rid of last '_' is Sup is empty


def data_with_long_id(data):
    keys = np.array([key(d) for d in data])
    b = rfn.append_fields(x, "Long_id", keys, usemask=False)
    return b


def get_y(data, reference):
    ref_str_to_i, i_to_ref_str = get_cross_index(reference["Atom"])
    # dat_str_to_i, i_to_dat_str = get_cross_index(data['Id'])

    nu1 = reference[ref_str_to_i(data["Atom1"])]["nu0"]
    nu2 = reference[ref_str_to_i(data["Atom2"])]["nu0"]

    y = (data["Value"] / nu1 * nu2 - 1).astype(float)
    yu = (data["Unc"] / nu1 * nu2).astype(float)

    return y, yu


def get_model_matrix(data, reference):
    ref_str_to_i, i_to_ref_str = get_cross_index(reference["Atom"])

    nMeas = len(data)
    nAtoms = len(reference)

    M = np.zeros((nMeas, nAtoms))

    index1 = ref_str_to_i(data["Atom1"])
    index2 = ref_str_to_i(data["Atom2"])

    M[np.arange(nMeas), index1] = 1
    M[np.arange(nMeas), index2] = -1

    # now the first column is related to Cs, that is non and adjusted freq
    # so get rid of it
    M = np.delete(M, 0, axis=1)

    return M


def get_corr_matrix(data, cor_data):
    nMeas = len(data)
    dat_str_to_i, i_to_dat_str = get_cross_index(data["Id"])

    cyy = np.eye(nMeas)

    # mask correlation values not requested by the input data
    # in case the same correlation data is used on a subset of input data
    mask = np.isin(cor_data["Id1"], data["Id"]) & np.isin(cor_data["Id2"], data["Id"])

    index1 = dat_str_to_i(cor_data["Id1"][mask])
    index2 = (dat_str_to_i(cor_data["Id2"][mask]),)

    cyy[index1, index2] = cor_data["Corr"][mask]

    # symmetrize
    cyy = cyy + cyy.T - np.diag(cyy.diagonal())

    return cyy