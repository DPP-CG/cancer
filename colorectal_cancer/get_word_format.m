function [word_format, inci_val] = get_word_format(polypOnsetRateAgeGroup, age_dist, AgeArrayLower, AgeArrayUpper)

% word_format is just for the onset rate (see "y" in the diagram for
% colorectal)
word_format = zeros(length(AgeArrayLower),1);

for i = 1:length(AgeArrayLower)
    word_format(i) = polypOnsetRateAgeGroup(i,1)*sum(age_dist(AgeArrayLower(i):AgeArrayUpper(i)))*1000;
end

inci_val = sum(sum(word_format(7:end)));

end