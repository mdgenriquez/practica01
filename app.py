import streamlit as st
import pandas as pd
from rdkit import Chem
from rdkit.Chem import Draw
from rdkit.Chem import AllChem
from rdkit.Chem import DataStructs
from rdkit.Chem.Fingerprints import FingerprintMols
import py3Dmol
from stmol import showmol
from io import BytesIO
from pubchempy import get_compounds, Compound
import sys
import pubchempy as pcp

# Configuración de la barra lateral
st.sidebar.image("plant.jpg") #, caption="Autor: Guadalupe Enriquez Asesor: -Dr.Jesus Alvarado"
st.sidebar.title("Predictor plant")
st.sidebar.markdown("Autor: Guadalupe Enriquez")


def load_data(file_path):
    
    try:
        return pd.read_csv(file_path)
    except Exception as e:
        st.error(f"Error al cargar el archivo CSV: {e}")
        return None

# Generar archivo SDF
def generate_sdf(mol):
    """Genera un archivo SDF 3D de la molécula dada."""
    mol = Chem.AddHs(mol)
    AllChem.EmbedMolecule(mol)
    AllChem.MMFFOptimizeMolecule(mol, maxIters=200)
    
    sdf_data = BytesIO()
    writer = Chem.SDWriter(sdf_data)
    writer.write(mol)
    writer.close()
    sdf_data.seek(0)
    return sdf_data


csv_path = "base_datos_plant_smile.csv"  # Nombre del archivo
data = load_data(csv_path)

if data is not None:
    # Validar columnas necesarias
    if "SMILES" not in data.columns or "Name" not in data.columns:
        st.error("El archivo CSV debe contener las columnas 'SMILES' y 'Name'.")
    else:
        st.title("Predictor plant")
        
        # Mostrar un slider para seleccionar una molécula
        molecule_index = st.slider("Selecciona una molécula", 0, len(data) - 1, 0)
        selected_row = data.iloc[molecule_index]
        
        # Extraer datos de la molécula seleccionada
        smiles = selected_row["SMILES"]
        name = selected_row["Name"]

        #st.markdown("### Nombre IUPAC")  
        #nombreiupac = pcp.get_compounds(smiles,'smiles')
        #st.text(nombreiupac[0].iupac_name)

        st.markdown("### Coeficiente de partición")
        coeficientedeparticion = get_compounds(smiles, 'smiles')
        st.text(coeficientedeparticion[0].xlogp)
        
        #st.subheader(f"ID: {name}")
        st.text(f"Código SMILES: {smiles}")
        
        # Generar la representación 2D de la molécula
        mol = Chem.MolFromSmiles(smiles)
        if mol:
            Draw.MolToFile(mol, "molecule.png")
            st.image("molecule.png", caption="Representación 2D de la molécula")
        else:
            st.error("No se pudo generar la representación molecular a partir del código SMILES.")
        
        # Visualización en 3D
        st.subheader("Visualización en 3D")
        def show_3d(smi):
            mol = Chem.MolFromSmiles(smi)
            mol = Chem.AddHs(mol)
            Chem.AllChem.EmbedMolecule(mol)
            Chem.AllChem.MMFFOptimizeMolecule(mol, maxIters=200)
            mblock = Chem.MolToMolBlock(mol)
            viewer = py3Dmol.view(width=500, height=400)
            viewer.addModel(mblock, "mol")
            viewer.setStyle({"stick": {}})
            viewer.zoomTo()
            showmol(viewer, height=400, width=500)

        show_3d(smiles)
        
# Calcular similitud de Tanimoto
def get_similar_compounds(query_smiles, df, threshold=0.7):
    query_mol = Chem.MolFromSmiles(query_smiles)
    query_fp = FingerprintMols.FingerprintMol(query_mol)
    
    similarities = []
    for index, row in df.iterrows():
        mol = Chem.MolFromSmiles(row["smiles"])
        if mol:
            fp = FingerprintMols.FingerprintMol(mol)
            similarity = DataStructs.FingerprintSimilarity(query_fp, fp)
            if similarity >= threshold:
                similarities.append((row["name"], row["smiles"], similarity))
    
    return sorted(similarities, key=lambda x: x[2], reverse=True)
    

else:
    st.warning("El archivo CSV no contiene datos válidos o no está cargado.")
