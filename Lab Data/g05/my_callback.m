block = find_system('dataRecorder/Recorder','Name','To File');

if ~isempty(block)
    file_name = strcat('FlightDynLabPt1_', datestr(now,'yyyy-mm-dd_HH-MM-SS'),'.mat')
    set_param(block{1}, 'Filename', file_name);
end