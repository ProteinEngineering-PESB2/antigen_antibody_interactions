import pandas as pd
from sklearn.feature_extraction.text import TfidfVectorizer
import sys
import clustering_algorithms
import evaluation_clustering

dataset = pd.read_csv(sys.argv[1])
path_export = sys.argv[2]

vec = TfidfVectorizer(stop_words="english")
vec.fit(dataset.description)
features = vec.transform(dataset.description)

matrix_summary = []
clustering_instance = clustering_algorithms.aplicateClustering(features)


print("Explore k-means")
for k in range(2,30):
    response = clustering_instance.aplicateKMeans(k)
    if response == 0:
        eval_clustering = evaluation_clustering.evaluationClustering(features.toarray(), clustering_instance.labels)
        labels = [str(value) for value in clustering_instance.labels]
        row_data = ["k-means", k, eval_clustering.calinski, eval_clustering.siluetas]
        matrix_summary.append(row_data)
print("Explore Birch")
for k in range(2,30):
    response = clustering_instance.aplicateBirch(k)
    if response == 0:
        eval_clustering = evaluation_clustering.evaluationClustering(features.toarray(), clustering_instance.labels)
        evaluation_cluster = {"algorithm": "Birch", "K": k, "calinski": eval_clustering.calinski, "siluetas": eval_clustering.siluetas}
        labels = [str(value) for value in clustering_instance.labels]
        evaluation_cluster.update({"labels":labels})
        row_data = ["Birch", k, eval_clustering.calinski, eval_clustering.siluetas]
        matrix_summary.append(row_data)

print("Apply DBScan")
response = clustering_instance.aplicateDBSCAN()
if response == 0:
    eval_clustering = evaluation_clustering.evaluationClustering(features.toarray(), clustering_instance.labels)
    evaluation_cluster = {"algorithm": "DBSCAN", "param": "-", "calinski": eval_clustering.calinski,
                          "siluetas": eval_clustering.siluetas}
    labels = [str(value) for value in clustering_instance.labels]
    evaluation_cluster.update({"labels": labels})
    row_data = ["DBScan", None, eval_clustering.calinski, eval_clustering.siluetas]
    matrix_summary.append(row_data)

print("Apply Affinity")
response = clustering_instance.aplicateAffinityPropagation()
if response == 0:
    eval_clustering = evaluation_clustering.evaluationClustering(features.toarray(), clustering_instance.labels)
    evaluation_cluster = {"algorithm": "Affinity", "param": "-", "calinski": eval_clustering.calinski,
                          "siluetas": eval_clustering.siluetas}
    labels = [str(value) for value in clustering_instance.labels]
    evaluation_cluster.update({"labels": labels})
    row_data = ["Affinity", None, eval_clustering.calinski, eval_clustering.siluetas]
    matrix_summary.append(row_data)

print("Apply Mean shift")
response = clustering_instance.aplicateMeanShift()
if response == 0:
    eval_clustering = evaluation_clustering.evaluationClustering(features.toarray(), clustering_instance.labels)
    evaluation_cluster = {"algorithm": "Mean Shift", "param": "-", "calinski": eval_clustering.calinski,
                          "siluetas": eval_clustering.siluetas}
    labels = [str(value) for value in clustering_instance.labels]
    evaluation_cluster.update({"labels": labels})
    row_data = ["Mean Shift", None, eval_clustering.calinski, eval_clustering.siluetas]
    matrix_summary.append(row_data)

print("Explore Agglomerative")
for linkage in ['ward', 'complete', 'average', 'single']:
    for affinity in ['euclidean', 'l1', 'l2', 'manhattan', 'cosine', 'precomputed']:
        for k in range(2, 30):
            response = clustering_instance.aplicateAlgomerativeClustering(linkage, affinity, k)
            if response == 0:
                eval_clustering = evaluation_clustering.evaluationClustering(features.toarray(),
                                                                             clustering_instance.labels)
                params = "{}-{}-{}".format(linkage, affinity, k)
                evaluation_cluster = {"algorithm": "Agglomerative", "param": params, "calinski": eval_clustering.calinski,
                                      "siluetas": eval_clustering.siluetas}
                labels = [str(value) for value in clustering_instance.labels]
                evaluation_cluster.update({"labels": labels})
                row_data = ["Agglomerative", params, eval_clustering.calinski, eval_clustering.siluetas]
                matrix_summary.append(row_data)

df_export = pd.DataFrame(matrix_summary, columns=['algorithm', 'params', 'calinski', 'siluetas'])
df_export.to_csv(path_export+"exploratory_clustering.csv", index=False)