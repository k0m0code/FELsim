def appendToList(self, xStd, yStd, xMean, yMean, x_axis, interval, matrixVariables):

        # Not used anymore
        # To be removed
        # Due to the use of pandas frame from the 6d distribution, we do no need to recalculate the standard devs
        '''
        Append updated values to five different arrays, used for plotBeamPositionTransform

        Parameters
        ----------
        xStd: list[float]
            standard deviation of particles' x position for each distance interval
        yStd: list[float]
            standard deviation of particles' y position for each distance interval
        xMean: list[float]
            average of particles' x position for each distance interval
        yMean: list[float]
            average of particles' y position for each distance interval
        x_axis: list[float]
            contains distance intervals to measure particle data over
        interval: float
            the amount between each distance interval
        matrixVariables: np.array(list[float][float])
            2D numPy list of particle elements to measure data from

        Returns
        -------
        xStd: list[float]
            updated standard deviation of x position list
        yStd: list[float]
            updated standard deviation of y position list
        xMean: list[float]
            updated average of x position list
        yMean: list[float]
            updated average of y positiion list
        x[axis]: list[float]
            updated distance interval list
        '''
        xStd.append(np.std(matrixVariables[:,0]))
        yStd.append(np.std(matrixVariables[:,2]))
        xMean.append(np.mean(matrixVariables[:,0]))
        yMean.append(np.mean(matrixVariables[:,2]))
        x_axis.append(round(x_axis[-1]+interval, self.DEFAULTINTERVALROUND))
        return xStd, yStd, xMean, yMean, x_axis
