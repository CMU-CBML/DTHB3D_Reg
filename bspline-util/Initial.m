function T1= Initial(knv1,knv2)
T1 = zeros(length(knv2)-1,length(knv1)-1);
for i = 1: length(knv2)-1,
    for j = 1:length(knv1)-1,
        if (knv2(i) >= knv1(j) && knv2(i) < knv1(j+1)),
            T1(i,j) = 1;
        else
            T1(i,j) = 0;
        end
    end
end
end