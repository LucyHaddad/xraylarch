import numpy as np
import re

class Structure3D():
    """""
    get 3d coordinates and paths associated with a given feff .inp
    args
    -----
    path: path to folder where feff has ran
    """""
    def __init__(self, path):
        self.path = path
        inp_path = f"{self.path}/feff.inp"
        with open(inp_path) as f:
            self.inp_text = f.read()
        self.inp_data, self.keys = self.parse_fefftext()
    

    def parse_fefftext(self):

        ftext = (self.inp_text.split("ATOMS")[1]).split("END")[0]
        ftext = ftext.split("\n")

        coords_idx, atoms_idx = ftext[1].index("ipot"), ftext[1].index("tag")
        ftext_coords, ftext_atoms = [l[:coords_idx] for l in ftext], [l[atoms_idx:atoms_idx+3] for l in ftext]
        coords_info = ftext_coords[2:-2]
        atoms_info = ftext_atoms[2:-2]

        atoms = re.findall("[A-Z][a-z]?[^_]", str(atoms_info))
        atoms = [(x.replace("'","")).replace(" ","") for x in atoms]
        atoms_list_enum = [a if not (s:=sum(j==a for j in atoms[:i])) else f"{a}{s+1}" for i,a in enumerate(atoms)]

        coords = np.array(re.findall("-?[0-9].\d{4}", str(coords_info)))
        coords = [float(i) for i in coords]
        split_coords = np.array_split(coords, len(atoms), axis=0)
        inp_data = dict(zip(atoms_list_enum, split_coords))

        keys = np.unique(atoms)

        return inp_data, keys


    def get_scattering_atoms(self, path_no):

        spath = f"{self.path}/{path_no}"

        with open(spath,"r") as f:
            sdata = f.read()
        sdata = (sdata.split('    k   real[2*phc]')[0]).split("edge")[1]

        coords = np.array(re.findall("-?[0-9].\d{4}", sdata))
        coords = [float(i) for i in coords]

        atoms_list = re.findall("[A-Z][a-z]?", sdata)
        atoms_list_enum = [a if not (s:=sum(j==a for j in atoms_list[:i])) else f"{a}{s+1}" for i,a in enumerate(atoms_list)]

        split_coords = np.array_split(coords, len(atoms_list), axis=0)
        scattering_coords = dict(zip(atoms_list_enum, split_coords))

        return scattering_coords
    
    def get_arrow_vectors(self, path_no):
        
        scattering_coords = self.get_scattering_atoms(path_no)
        atoms_list = list(scattering_coords.keys())
        coords_list = []
        for i in range(len(atoms_list)):
            x0, y0, z0 = scattering_coords[atoms_list[i]][0], scattering_coords[atoms_list[i]][1], scattering_coords[atoms_list[i]][2]
            x1, y1, z1 = 0, 0, 0
            if i < len(atoms_list)-1:
                x1, y1, z1 = scattering_coords[atoms_list[i+1]][0], scattering_coords[atoms_list[i+1]][1], scattering_coords[atoms_list[i+1]][2]
            if i == len(atoms_list)-1:
                x1, y1, z1 = scattering_coords[atoms_list[0]][0], scattering_coords[atoms_list[0]][1], scattering_coords[atoms_list[0]][2]
            dx, dy, dz = x1-x0,y1-y0,z1-z0
            coords_list = np.concatenate([coords_list, [x0,y0,z0,dx,dy,dz]])

        return coords_list