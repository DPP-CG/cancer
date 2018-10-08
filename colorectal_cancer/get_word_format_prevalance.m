function [word_format_prevalence, prev_val] = get_word_format_prevalance(prevalence,PAR,AgeArrayLower,AgeArrayUpper)

% word_format is just for the prevalence 
word_format_prevalence = zeros(size(prevalence));
age_dist = PAR{1};

for i = 1:length(AgeArrayLower)
    word_format_prevalence(i,:) = prevalence(i,:)*sum(age_dist(AgeArrayLower(i):AgeArrayUpper(i)));
end

prev_val = sum(sum(word_format_prevalence(7:end,:)));

end