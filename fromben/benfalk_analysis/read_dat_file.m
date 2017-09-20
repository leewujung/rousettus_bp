[filename,directory_name] = uigetfile('*.dat', 'Select a file');
fullname = fullfile(directory_name, filename);
data = load(fullname);