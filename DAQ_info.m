%% DAQ info                                                     02/04/2018
%  - This script obtains the information from a connected DAQ and stores
%  the following values in a .mat file:
%    - Vendor
%    - Model
%    - Subsystems

function DAQ_info(name)
    info = daq.getDevices();
    
    device.vendor = info.Vendor;
    device.model = info.Model;
    device.subsystems = info.Subsystems;
    device.subsystems(1,1) = 
    
    save(name, 'device');
end
