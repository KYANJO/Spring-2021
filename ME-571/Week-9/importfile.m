function [x, y, u, u_exact] = importfile(filename, dataLines)
%IMPORTFILE Import data from a text file
%  [X, Y, U, U_EXACT] = IMPORTFILE(FILENAME) reads data from text file
%  FILENAME for the default selection.  Returns the data as column
%  vectors.
%
%  [X, Y, U, U_EXACT] = IMPORTFILE(FILE, DATALINES) reads data for the
%  specified row interval(s) of text file FILENAME. Specify DATALINES as
%  a positive scalar integer or a N-by-2 array of positive scalar
%  integers for dis-contiguous row intervals.
%
%  Example:
%  [x, y, u, u_exact] = importfile("/Volumes/GoogleDrive/My Drive/teaching/ME471:571/ME471_571 Spring 2021/Labs/Week 9/results_poisson_serial.dat", [2, Inf]);
%
%  See also READTABLE.
%
% Auto-generated by MATLAB on 09-Mar-2021 22:56:43

%% Input handling

% If dataLines is not specified, define defaults
if nargin < 2
    dataLines = [2, Inf];
end

%% Setup the Import Options and import the data
opts = delimitedTextImportOptions("NumVariables", 4);

% Specify range and delimiter
opts.DataLines = dataLines;
opts.Delimiter = ",";

% Specify column names and types
opts.VariableNames = ["x", "y", "u", "u_exact"];
opts.VariableTypes = ["double", "double", "double", "double"];

% Specify file level properties
opts.ExtraColumnsRule = "ignore";
opts.EmptyLineRule = "read";

% Import the data
tbl = readtable(filename, opts);

%% Convert to output type
x = tbl.x;
y = tbl.y;
u = tbl.u;
u_exact = tbl.u_exact;
end