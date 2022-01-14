import pandas as pd
import sys
from sklearn import preprocessing
import joblib
from sklearn.decomposition import IncrementalPCA

#dataset a procesar, debe tener todos los componentes (antigeno, anticuerpo L y anticuerpo H) sin respuesta ni ID de interaccion o ID de secuencias
dataset = pd.read_csv(sys.argv[1])
path_output = sys.argv[2]

#apply scale
print("Scale dataset")
min_max_scaler = preprocessing.MinMaxScaler()
dataset_scale = min_max_scaler.fit_transform(dataset)

#Descomprimir cuando no exista el instance model de scaler
#print("Export scaler instance")
#name_export_scaler = "{}scaler_model.joblib".format(path_output)
#joblib.dump(min_max_scaler, name_export_scaler)

print("IncrementalPCA process")
pca = IncrementalPCA()
pca = pca.fit(dataset_scale)

dataset_pca = pca.fit_transform(dataset_scale)

print("Export data and IncrementalPCA-model")
header = ["P_{}".format(i+1) for i in range(len(dataset_pca[0]))]
df_export = pd.DataFrame(dataset_pca, columns=header)
df_export.to_csv("{}IncrementalPCA_results.csv".format(path_output), index=False)
joblib.dump(pca, "{}IncrementalPCA_model.joblib".format(path_output))