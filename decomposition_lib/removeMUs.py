# Generated with SMOP  0.41
from libsmop import *
# ..\decompositionproject_python\decomposition_lib\removeMUs.m

    
@function
def removeMUs(signal=None,decompParameters=None,fs=None,ids_remove=None,manual_removal=None,*args,**kwargs):
    varargin = removeMUs.varargin
    nargin = removeMUs.nargin

    # removes mus in signal and decompoParamerters according to logical input
# ids_remove
    
    if logical_not(isempty(ids_remove)):
        structs_in=cellarray([signal,decompParameters])
# ..\decompositionproject_python\decomposition_lib\removeMUs.m:6
        structs_out=removeFromStruct(structs_in,ids_remove)
# ..\decompositionproject_python\decomposition_lib\removeMUs.m:7
        signal=structs_out[1]
# ..\decompositionproject_python\decomposition_lib\removeMUs.m:9
        decompParameters=structs_out[2]
# ..\decompositionproject_python\decomposition_lib\removeMUs.m:10
    
    if manual_removal:
        ids_remove=false(concat([1,min(size(signal.spikeTrains))]))
# ..\decompositionproject_python\decomposition_lib\removeMUs.m:14
        mnIndexTicks=cellarray([cellarray([]),cellarray(['',num2str(ravel(signal.COV),'%.2f')])])
# ..\decompositionproject_python\decomposition_lib\removeMUs.m:15
        labels=cellarray(['Time (s)','Motor neuron (#)','CoV','Torque (Nm)'])
# ..\decompositionproject_python\decomposition_lib\removeMUs.m:16
        figure
        ax=plotMUs(signal,fs,'flagPlotSpikeTrains',0,'flagPlotDR',1,'mnIndexTicks',mnIndexTicks,'labels',labels)
# ..\decompositionproject_python\decomposition_lib\removeMUs.m:18
        title('Clic on MU to remove')
        __,y=ginput
# ..\decompositionproject_python\decomposition_lib\removeMUs.m:21
        close_('gcf')
        removedMUs=round(unique(y.T))
# ..\decompositionproject_python\decomposition_lib\removeMUs.m:24
        ids_remove[removedMUs]=true
# ..\decompositionproject_python\decomposition_lib\removeMUs.m:25
        structs_in=cellarray([signal,decompParameters])
# ..\decompositionproject_python\decomposition_lib\removeMUs.m:26
        structs_out=removeFromStruct(structs_in,ids_remove)
# ..\decompositionproject_python\decomposition_lib\removeMUs.m:27
        signal=structs_out[1]
# ..\decompositionproject_python\decomposition_lib\removeMUs.m:28
        decompParameters=structs_out[2]
# ..\decompositionproject_python\decomposition_lib\removeMUs.m:29
        mnIndexTicks=cellarray([cellarray([]),cellarray(['',num2str(ravel(signal.COV),'%.2f')])])
# ..\decompositionproject_python\decomposition_lib\removeMUs.m:30
        ax=plotMUs(signal,fs,'flagPlotSpikeTrains',0,'flagPlotDR',1,'mnIndexTicks',mnIndexTicks,'labels',labels)
# ..\decompositionproject_python\decomposition_lib\removeMUs.m:32
    
    return signal,decompParameters
    
if __name__ == '__main__':
    pass
    
    ## remove from struct
    
@function
def removeFromStruct(structs_in=None,ids_remove=None,*args,**kwargs):
    varargin = removeFromStruct.varargin
    nargin = removeFromStruct.nargin

    # Number of input structs
    num_structs=length(structs_in)
# ..\decompositionproject_python\decomposition_lib\removeMUs.m:43
    # Initialize output structs
    structs_out=copy(structs_in)
# ..\decompositionproject_python\decomposition_lib\removeMUs.m:46
    # Loop through each struct
    for s in arange(1,num_structs).reshape(-1):
        # Get the current struct and its fields
        input_struct=structs_in[s]
# ..\decompositionproject_python\decomposition_lib\removeMUs.m:51
        fields=fieldnames(input_struct)
# ..\decompositionproject_python\decomposition_lib\removeMUs.m:52
        for i in arange(1,length(fields)).reshape(-1):
            field_data=getattr(input_struct,(fields[i]))
# ..\decompositionproject_python\decomposition_lib\removeMUs.m:56
            data_size=size(field_data)
# ..\decompositionproject_python\decomposition_lib\removeMUs.m:57
            if length(ids_remove) == data_size(1):
                # Remove rows
                setattr(structs_out[s],fields[i],field_data(logical_not(ids_remove),arange()))
# ..\decompositionproject_python\decomposition_lib\removeMUs.m:62
            else:
                if length(ids_remove) == data_size(2):
                    # Remove columns
                    setattr(structs_out[s],fields[i],field_data(arange(),logical_not(ids_remove)))
# ..\decompositionproject_python\decomposition_lib\removeMUs.m:65
                else:
                    if length(data_size) > 2:
                        if length(ids_remove) == data_size(3):
                            # Remove 3 dim e.g. eSIG
                            setattr(structs_out[s],fields[i],field_data(arange(),arange(),logical_not(ids_remove)))
# ..\decompositionproject_python\decomposition_lib\removeMUs.m:69
                        #         else
#             warning('Length of ids_remove does not match any dimension of the field "#s" in struct #d. Skipping.', fields{i}, s);
    
    return structs_out
    
if __name__ == '__main__':
    pass
    