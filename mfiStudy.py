import pandas as pd
import numpy as np
import random
import math
from scipy import stats

from sklearn import preprocessing

class mfiStudy:
    """
        mfiStudy class introduces all the data parsing concepts presented in the paper
        
        Class parameters:
            self.X: Actual data
            self.Y: Labels (if any - or None)
            self.Z: Headers (if any or None) = Name for each dimension
            
        Supported methods (public):
            readTable: Get the columns of an excel/csv file
            prepareCrossValidation: Split a set into training and testing subsets.
            coResponsesValidation: Assess the used cut-off detection method by using a cross validation approach, 
                                   based on the antigen immune co-response.
    """

    #
    # Private
    #
    """
        Method: _distDivergence
    """ 
    def _distDivergence(self, _model1, _model2):
        """
            Compute the divergence between the two probability models (Variation of Jensen Shannon)
            
            _model1: The basic model.
            _model2: The test models.

            RETURN: The divergence score [0,1]
        """
        
        POINT_DENSITY = 10 # how many points to fit

        x = np.divide(list(range(0,(POINT_DENSITY + 1))), POINT_DENSITY-1)
        
        # basis fit
        p = stats.invgauss(mu=_model1[0], loc=_model1[1], scale=_model1[2]).pdf(x)
        p[p < .000000000000001] = .000000000000001
        p1 = np.array(np.divide(p, np.sum(p)))


        print(p1)
        avgP = p1

        # test fit
        p2 = np.zeros((_model2.shape[0],POINT_DENSITY+1))
        for i in range(0,_model2.shape[0]):
            p = stats.invgauss(mu=_model2[i,0], loc=_model2[i,1], scale=_model2[i,2]).pdf(x)
            p[p < .000000000000001] = .000000000000001
            p2[i,:] = np.array(np.divide(p, np.sum(p)))

            avgP = avgP + p2[i,:]

        # compute average
        avgP = avgP/(1+_model2.shape[0])
    
        print(avgP)
        gJS = np.sum(p1*np.log(np.divide(p1,avgP)))
        for i in range(0,_model2.shape[0]):
            gJS = gJS + np.sum(p2*np.log(np.divide(p2[i,:],avgP)))
        gJS = gJS/(1+_model2.shape[0])

        print(gJS)

        return (math.exp(-gJS))
    
    #
    # Public
    #

    """
        Method: readTable
    """ 
    def readTable(self, _path, _colsToRead=None, _linesToSkip=None, _headersRow=False, _labelsCol=None, _controlCol=None, _controlThreshold=None, _cleanCharacters=[]):
        """
            Read an excel/csv file to extract the table data.
            
            _path (str): FULL path of the excel file to read - Read permission is assumed
            _colsToRead (str): Columns to read, eg. D:AF meaning read cols from D to AF. None to read all.
                               IMPORTANT: if _headersRow is True then make sure the headers row is also include.
            _linesToSkip (int): How many lines to skip reading (ignored). Use None to include all lines.
                                IMPORTANT: if _linesToSkip is provided and _headersRow is True
                                           then make sure the ignored rows include the header row
            _headersRow (bool): TRUE = Use the first line as headers. Default is FALSE.
            _labelsCol (str): Which col to use for getting data labels. Default None = No labels.
            _controlCol (str): Define a column that will be used to threshold the input data
            _controlThreshold (float): The threshold that will be used when checking the _controlCol.
                                       If _controlCol[at position i] < _controlThreshold, then the row is rejected
            _cleanCharacters (str array): Remove these characters from header title.
        """    
        # Read the actual data
        if (_linesToSkip is None):
            if (_headersRow):
                #first row includes the headers of the dimensions
                self.X = pd.read_excel (_path, usecols=_colsToRead, skiprows=1).values                
            else:
                self.X = pd.read_excel (_path, usecols=_colsToRead).values
        else:
            self.X = pd.read_excel (_path, usecols=_colsToRead, skiprows=_linesToSkip).values

        # Normalize
        self.X = preprocessing.normalize(self.X, norm='l2')

        # Read the data headers (if any)
        if (_headersRow):
            self.Z = pd.read_excel (_path, usecols=_colsToRead, dtype=str, nrows=1, header=None).values
            self.Z = np.ndarray.flatten(self.Z[0])
            
            if (self.Z.shape[0] != self.X.shape[1]):
                raise NameError("Number of headers should match the dimension of data.")

            # Check for cleaning characters
            for s in _cleanCharacters:
                self.Z = np.array([z.replace(s, '') for z in self.Z])
        else:
            self.Z = None
            
        # Read the data labels (if any)
        if (_labelsCol is not None):
            self.Y = pd.read_excel (_path, usecols=_labelsCol, skiprows=_linesToSkip).values
            self.Y = np.ndarray.flatten(self.Y)
            
            if (self.Y.shape[0] != self.X.shape[0]):
                raise NameError("Number of labels should match the number of data.")
        else:
            self.Y = None
        
        # Thresholding
        if (_controlCol is not None and _controlThreshold is not None):
            cutCol = pd.read_excel (_path, usecols=_controlCol, skiprows=_linesToSkip).values
            ind = np.where(cutCol >= _controlThreshold)
            self.X = self.X[ind,:]
            if (self.Y is not None):
                self.Y = self.Y[ind,:]
    
    """
        Method: prepareCrossValidation
    """ 
    def prepareCrossValidation(self, _trainPct=.9, _allInd=None):
        """
            Suffle the training set and prepare two sets: Training & Testing
            
            _trainPct: (% in (0,1)) - percentage of data to use for the training set
                        Default is 10-fold cross validation
            _Ind: Return from these indices (If None, use all indices)

            RETURN: (trainInd, testInd): The indices of self.X to contain the training and testing observations respectively.
        """
        if _trainPct <=0 or _trainPct >=1:
            _trainPct = .9 # default is 10-fold cross validation
        
        if _allInd is None:
            _allInd = list(range(0,self.X.shape[0]))
        random.shuffle(_allInd)

        breakInd = math.floor(_trainPct*len(_allInd))
        trainInd = _allInd[0:breakInd]
        testInd  = _allInd[breakInd:]

        return(trainInd, testInd)
        
    """
        Method: coResponsesValidation
    """
    def coResponsesValidation(self, _ColBasis, _ColControl, _allInd=None):
        """
            Co-response cross validation of the detected cut-offs. Assess the quality of the
            detected cut-off according to the co-responses of the defined antigens.

            The method relies on the observation that there are certain rules regarding the
            response of groups of antigens (correlated or anti-correlated). Under that perspective,
            the studied cut-off detection method is applied on that list of grouped antigens and then,
            we compute if the defined cut-offs detect the anticipated co-expression.
            
            _ColBasis: The column of the matrix that will be used as basis to fit a model.
            _ColControl: An array of columns that will be used to test the model fitting (hypothesis testing).
            _allInd: The records to be used for testing. None (default) to use all records.

            Return: Similarity percentage
        """
        if _allInd is None:
            _allInd = list(range(0,self.X.shape[1]))

        # fit on observation data
        fitBasis = stats.invgauss.fit(self.X[_allInd, _ColBasis])

        # fit on test data
        fitTest = np.zeros((len(_ColControl),3))
        for i in range(0, len(_ColControl)):
            # fit on test data
            (mu, loc, scale) = stats.invgauss.fit(self.X[_allInd, _ColControl[i]])
            fitTest[i:,] = [mu, loc, scale]
        
        # hyphothesis testing
        d = self._distDivergence(fitBasis, fitTest)

        return (d)