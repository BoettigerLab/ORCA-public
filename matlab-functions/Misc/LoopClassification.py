def loopClassification(data, Mb, MATLAB=True):
    ## TODO: add sequential type (jet)
    ## TODO: make sure that cross and multicontact are mutually exclusive
    """
    loopClassification

    Parameters
    ----------
    data : a list or an array-like structure
        A list of coordinates of chromatin contacts (loop anchors) in a single cell
        These coordinates are all from lower triangle (ie, if coord = [x, y], x >= y)

    Mb : an integer
        The number of readout probes

    Return
    ------
    dict
        a dictionary contains all loop shapes in this single cell
        dict["nestedLoopDict"] is a dictionary whose keys are the loop base and values are all smaller loops
        dict["multiContactSet"] is a set containing the readout probes that are multicontact hubs
        dict["multiContactGroup"] is a list of lists. Each nested list contains coordinates that share a similar
            row or column
        dict["isolatedLoopSet"] is a set of isolated loops
        dict["crossLoopDict"] is a dictionary of cross loop whose keys are the reference loop and values are
            all the cross loops that are detected based on that reference loop.
        dict["numMultiContact"] is an integer represents how many unique readout probes that are multicontact hubs
    """
    import numpy as np

    if MATLAB:
        # convert to integer and zero-based indexing
        data = np.intc(np.around(data, 0) - 1)
        Mb = np.intc(Mb)
    else:
        data = np.intc(np.around(data, 0))

    # initialize variables
    nestedLoopDict = dict()
    crossLoopDict = dict()
    contactMatrix = np.zeros((Mb, Mb))

    # iterate through each loop anchor
    for loopCoord in data:
        x, y = loopCoord
        k = (x, y)
        contactMatrix[y, x] += 1
        # iterate through each loop anchor, called this loopCoord2
        for loopCoord2 in data:
            x2, y2 = loopCoord2
            # if loopCoord2 is inside loopCoord, loopCoord is said to be a nested loop where
            # loopCoord2 is a smaller loop inside it.
            if y < x2 < x and y < y2 < x:
                if k not in nestedLoopDict:
                    nestedLoopDict[k] = [(x2, y2)]
                else:
                    nestedLoopDict[k].append((x2, y2))

            # if loopCoord and loopCoord2 share some area under the dot, loopCoord and loopCoord2 are cross loops
            # the first set of logical operations make sure that loopCoord and loopCoord2 are not multicontact
            # the seconde set of logical operations make sure than loopCoord and loopCoord2 only partially share
            # the area under the dot using XOR operation (^).
            if (x2 != x and y2 != y and x2 != y and y2 != x) and ((y < x2 < x) ^ (y < y2 < x)):
                if k not in crossLoopDict:
                    crossLoopDict[k] = [(x2, y2)]
                else:
                    crossLoopDict[k].append((x2, y2))

    # find rows and columns that sum > 1
    rowSum = np.sum(contactMatrix, axis=1)
    colSum = np.sum(contactMatrix, axis=0)

    multiContactRow = np.arange(Mb)[rowSum > 1]
    multiContactCol = np.arange(Mb)[colSum > 1]

    multiContactSet = set()
    multiContactGroup = []

    for r in multiContactRow:
        coords = data[data[:, 1] == r]
        for coord in coords:
            multiContactSet.add(tuple(coord))
        multiContactGroup.append(coords)

    for c in multiContactCol:
        coords = data[data[:, 0] == c]
        for coord in coords:
            multiContactSet.add(tuple(coord))
        multiContactGroup.append(coords)

    # the remaining loop anchors are isolated
    isolatedLoopSet = set()
    for loopCoord in data:
        x, y = loopCoord
        k = (x, y)
        if k not in nestedLoopDict and k not in multiContactSet and k not in crossLoopDict:
            isolatedLoopSet.add(k)

    return {"nestedLoopDict": nestedLoopDict,
            "multiContactSet": multiContactSet,
            "multiContactGroup": multiContactGroup,
            "isolatedLoopSet": isolatedLoopSet,
            "crossLoopDict": crossLoopDict,
            "numMultiContact": len(set(multiContactRow) | set(multiContactCol))
            }


loopType = loopClassification(data, Mb)