function save2file(array, filename)

s = size(array);

if(or(s(1) == 1,s(2) == 1))
    n = length(array);
    f = fopen(filename, 'w+');
    for I = 1:n
        st = strtrim(sprintf('%22.16g', array(I)));
%        fprintf(f, '%22.16e\n', array(I));
        fprintf(f, '%s\n', st);
    end
    fclose(f);
    return
end

n1 = s(1);
n2 = s(2);

f = fopen(filename, 'w+');
for I = 1:n1
    for J = 1:n2
        st = strtrim(sprintf('%22.16g', array(I,J)));
%        fprintf(f, '%22.16e\t', array(I,J));
        fprintf(f, '%s\t', st);
    end
    fprintf(f, '\n');
end

fclose(f);

return