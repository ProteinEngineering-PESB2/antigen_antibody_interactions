#Import all the dependencies
from gensim.models.doc2vec import Doc2Vec, TaggedDocument
from nltk.tokenize import word_tokenize
import pandas as pd
import sys
import clustering_algorithms
import evaluation_clustering

print("Reading data")
dataset = pd.read_csv(sys.argv[1])
path_export = sys.argv[2]

tagged_data = [TaggedDocument(words=word_tokenize(_d.lower()), tags=[str(i)]) for i, _d in enumerate(dataset.description)]

#define parameters
vec_size = [i for i in range(2, 30)]
alpha = 0.025

for vec in vec_size:
    print("Iteration: ", vec)
    print("Creating data")
    model = Doc2Vec(vector_size=vec, alpha=alpha, min_alpha=0.00025, min_count=1, dm=1)
    model.build_vocab(tagged_data)

    print("Training encoder")
    model.train(tagged_data, total_examples=model.corpus_count, epochs=500)


    print("Get embedding")
    embedding_matrix = []

    for i in range(len(tagged_data)):
        embedding = model.infer_vector(tagged_data[i][0])
        embedding_matrix.append(embedding)

    print("Export embedding")
    header = ["p_{}".format(i) for i in range(vec)]
    df_export = pd.DataFrame(embedding_matrix, columns=header)
    df_export.to_csv("{}{}_embedding_data.csv".format(path_export,vec), index=False)

    clustering_instance = clustering_algorithms.aplicateClustering(df_export)
    matrix_summary = []

    print("Explore k-means")
    for k in range(2,30):
        response = clustering_instance.aplicateKMeans(k)
        if response == 0:
            eval_clustering = evaluation_clustering.evaluationClustering(df_export, clustering_instance.labels)
            labels = [str(value) for value in clustering_instance.labels]
            row_data = ["k-means", k, eval_clustering.calinski, eval_clustering.siluetas]
            matrix_summary.append(row_data)

    print("Explore Birch")
    for k in range(2,30):
        response = clustering_instance.aplicateBirch(k)
        if response == 0:
            eval_clustering = evaluation_clustering.evaluationClustering(df_export, clustering_instance.labels)
            evaluation_cluster = {"algorithm": "Birch", "K": k, "calinski": eval_clustering.calinski, "siluetas": eval_clustering.siluetas}
            labels = [str(value) for value in clustering_instance.labels]
            evaluation_cluster.update({"labels":labels})
            row_data = ["Birch", k, eval_clustering.calinski, eval_clustering.siluetas]
            matrix_summary.append(row_data)

    print("Apply DBScan")
    response = clustering_instance.aplicateDBSCAN()
    if response == 0:
        eval_clustering = evaluation_clustering.evaluationClustering(df_export, clustering_instance.labels)
        evaluation_cluster = {"algorithm": "DBSCAN", "param": "-", "calinski": eval_clustering.calinski,
                              "siluetas": eval_clustering.siluetas}
        labels = [str(value) for value in clustering_instance.labels]
        evaluation_cluster.update({"labels": labels})
        row_data = ["DBScan", None, eval_clustering.calinski, eval_clustering.siluetas]
        matrix_summary.append(row_data)

    print("Apply Affinity")
    response = clustering_instance.aplicateAffinityPropagation()
    if response == 0:
        eval_clustering = evaluation_clustering.evaluationClustering(df_export, clustering_instance.labels)
        evaluation_cluster = {"algorithm": "Affinity", "param": "-", "calinski": eval_clustering.calinski,
                              "siluetas": eval_clustering.siluetas}
        labels = [str(value) for value in clustering_instance.labels]
        evaluation_cluster.update({"labels": labels})
        row_data = ["Affinity", None, eval_clustering.calinski, eval_clustering.siluetas]
        matrix_summary.append(row_data)

    print("Apply Mean shift")
    response = clustering_instance.aplicateMeanShift()
    if response == 0:
        eval_clustering = evaluation_clustering.evaluationClustering(df_export, clustering_instance.labels)
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
                    eval_clustering = evaluation_clustering.evaluationClustering(df_export,
                                                                                 clustering_instance.labels)
                    params = "{}-{}-{}".format(linkage, affinity, k)
                    evaluation_cluster = {"algorithm": "Agglomerative", "param": params, "calinski": eval_clustering.calinski,
                                          "siluetas": eval_clustering.siluetas}
                    labels = [str(value) for value in clustering_instance.labels]
                    evaluation_cluster.update({"labels": labels})
                    row_data = ["Agglomerative", params, eval_clustering.calinski, eval_clustering.siluetas]
                    matrix_summary.append(row_data)

    df_export = pd.DataFrame(matrix_summary, columns=['algorithm', 'params', 'calinski', 'siluetas'])
    df_export.to_csv("{}{}_exploratory_clustering.csv".format(path_export, vec), index=False)