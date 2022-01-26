
from sklearn import metrics

class evaluationClustering(object):

    def __init__(self, dataSet, labelsResponse):

        self.dataSet = dataSet
        self.labelsResponse = labelsResponse
        try:
            self.calinski = metrics.calinski_harabasz_score(self.dataSet, self.labelsResponse)
            self.siluetas = metrics.silhouette_score(self.dataSet, self.labelsResponse, metric='euclidean')
        except:
            self.calinski = "ERROR"
            self.siluetas = "ERROR"
            pass