import pandas as pd
from multiprocessing import Pool
import os
from selenium import webdriver
import requests
import time 
chromedriver_path = "/usr/bin/chromedriver"

def alphafold(file):
    data = pd.read_csv(file)
    uniprot_id = data.loc[0].Uniprot_id.split(".")[0]
    
    if(uniprot_id + ".pdb" not in os.listdir("alphafold_pdb")):
        print("buscando "+ uniprot_id)
        op = webdriver.ChromeOptions()
        op.add_argument('headless')#Para no ver la ventana del navegador (recomendado)
        driver = webdriver.Chrome(chromedriver_path, options=op)
        url = "http://alphafold.ebi.ac.uk/entry/" + uniprot_id
        driver.get(url)
        time.sleep(1)
        link = driver.find_element_by_partial_link_text("PDB file").get_attribute('href')
        response = requests.get(link)
        f = open("alphafold_pdb/"+ uniprot_id + ".pdb", "w")
        f.write(response.text)
        f.close()
        print(uniprot_id + " pdb exportado")
        driver.close()
    else:
        print("ya se encuentra "+ uniprot_id)

def process(file):
    try:
        alphafold("./uniprot/"+file)
    except:
        f = open("Errores.txt", "a")
        f.write(file)

if __name__ == '__main__':
    uniprot_files = os.listdir("./uniprot/")
    #print(uniprot_files)
    with Pool(10) as p:
        p.map(process, uniprot_files)