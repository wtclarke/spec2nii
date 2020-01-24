import json
import os.path as op

def writeJSON(filename,outdir,metaDict):
    # Form full path
    fullfilepath = op.join(outdir,filename+'.json')

    with open(fullfilepath, 'w') as fp:
        json.dump(metaDict, fp)
    

