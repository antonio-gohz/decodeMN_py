# Generated with SMOP  0.41
from libsmop import *
# ..\decompositionproject_python\decomposition_lib\parsave.m

    
@function
def parsave(fname=None,struc=None,*args,**kwargs):
    varargin = parsave.varargin
    nargin = parsave.nargin

    cellfun(lambda x=None: assignin('caller',x,getattr(struc,(x))),fieldnames(struc))
    filenames=fieldnames(struc)
# ..\decompositionproject_python\decomposition_lib\parsave.m:3
    save(fname,filenames[arange()])
    
    return
    
if __name__ == '__main__':
    pass
    