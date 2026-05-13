# pylint: disable=missing-module-docstring
import requests


def fetch(pdb_id, path=None):
    """ Fetches the PDB file for the given PDB ID.
    :param pdb_id: The PDB ID to fetch.
    :param path: The path to save the file. Defaults to the PDB ID plus '.pdb'.
    """
    pdb_id = pdb_id.upper()
    pdb_filename = pdb_id + '.pdb'
    url = "https://files.rcsb.org/download/" + pdb_filename
    data = requests.get(url, timeout=10)
    if data.status_code != 200:
        raise FileNotFoundError(f"PDB ID {pdb_id} dont exist.")
    if path is None:
        path = pdb_filename
    with open(path, 'wb') as f:
        f.write(data.content)
