function outstr = pointToDashed(parameter)
parastring = num2str(parameter);
if contains(parastring,'.')
    strarr = split(parastring,'.');
    outstr = join(strarr,'-');
    outstr = char(outstr);
else
    outstr = char(parastring);
end
end