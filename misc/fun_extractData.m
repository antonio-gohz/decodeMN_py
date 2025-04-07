function opData = fun_extractData(filename)
biomechproData = importdata(filename);

%% extract Marker data
% markers needed: 'RCAL'
mkrData  = biomechproData.Marker.MarkerData;
mkrLabel = biomechproData.Marker.MarkerDataLabel;

rcalInd  = find(contains(mkrLabel,'RCAL'));
rcalData = squeeze(mkrData(:,rcalInd,:));
opData.markers.data{1}  = rcalData;

lcalInd  = find(contains(mkrLabel,'LCAL'));
lcalData = squeeze(mkrData(:,lcalInd,:));
opData.markers.data{2}  = lcalData;

opData.markers.label    = {'RCAL','LCAL'};
%% extract IK data
% joints needed: 'ankle_angle_r'
ikData  = biomechproData.Marker.JointAngData;
ikLabel = biomechproData.Marker.JointAngDataLabel;

rankleAngleInd  = find(contains(ikLabel,'AJCR'));
rankleAngleData = squeeze(ikData(:,rankleAngleInd,:));
opData.ik.data{1}  = rankleAngleData;

lankleAngleInd  = find(contains(ikLabel,'AJCL'));
lankleAngleData = squeeze(ikData(:,lankleAngleInd,:));
opData.ik.data{2}  = lankleAngleData;

opData.ik.label    = {'ankle_angle_r','ankle_angle_l'};


%% extract ID data
% joints needed: 'ankle_angle_r_moment'
idData  = biomechproData.Marker.JointTrqData;
idLabel = biomechproData.Marker.JointTrqDataLabel;

rankleMomInd  = find(contains(idLabel,'AJCR'));
rankleMomData = squeeze(idData(:,rankleMomInd,:));
opData.id.data{1}  = rankleMomData;

lankleMomInd  = find(contains(idLabel,'AJCL'));
lankleMomData = squeeze(idData(:,lankleMomInd,:));
opData.id.data{2}  = lankleMomData;


opData.id.label    = {'ankle_angle_r_moment','ankle_angle_l_moment'};


%%
frameRate  = biomechproData.Marker.MarkerFrameRate;
dataLength = length(opData.markers.data{1});
dT = (1/frameRate);
opData.time = (dT:dT:(dataLength/frameRate))-dT;

end