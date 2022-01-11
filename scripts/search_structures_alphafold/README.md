# Buscador de secuencias utilizando Uniprot y Alphafold2.

## Requisitos:

 > pandas (pip install pandas)
 
 > biopython (pip install biopython)
 
 > numpy (pip install numpy)
 
 > selenium (pip install selenium)

## Ejecución:
Ejecutar en orden

 > python3 search_uniprot.py test.csv

 > python3 search_alphafold.py

 > python3 Build_Dataset.py test.csv

## Consideraciones:

Selenium requiere de drivers para el Browser. En el script search_alphafold.py se utiliza el driver de Chrome, ubicado en "/usr/bin/chromedriver". 

El archivo .csv con las secuencias tiene que tener las columnas [id,sequence].

Archivo de prueba: test.csv
