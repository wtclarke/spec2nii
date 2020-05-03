import json
import os.path as op
import numpy as np

def writeJSON(filename,outdir,metaDict):
    # Form full path
    fullfilepath = op.join(outdir,filename+'.json')

    # Clean up dict
    for key in metaDict:
        if isinstance(metaDict[key],np.ndarray):
            metaDict[key] = metaDict[key].tolist()

        # print(f'{key} is type: {type(metaDict[key])}.')

    with open(fullfilepath, 'w') as fp:
        json.dump(metaDict, fp)
    

