'''docstring'''
from scipy import signal

class PeakIdentifier():
    """docstring for main"""
    def __init__(self, counts_table):
        #geno_table = dict(sorted(geno_table.items(), key=lambda x: x[1][0], reverse=True))
        counts_table = dict(sorted(counts_table.items(), key=lambda x: x[0], reverse=False))

        self.counts_array = self.get_counts_array(counts_table)
        #peaks = self.get_peaks()
        #self.get_alleles_using_peaks_from_genotable()

    def get_peaks(self):
        average = sum(self.counts_array)/len(self.counts_array)
        #print(average)
        #print(signal.find_peaks(self.counts_array, threshold=average, distance=1))
        ans = signal.find_peaks(self.counts_array, prominence=average)
        #print(self.counts_array)
        #print(ans[0])
        return(ans[0])

    def get_counts_array(self,counts_table):
        ''' get array from counts_table dictionary, the array index is the dictionary
            key (repeat counts), zero padding is done
        '''
        counts_array = []

        for key in counts_table.keys():
            while key > len(counts_array):
                counts_array.append(0)
            counts_array.append(counts_table[key])
        return(counts_array)      
