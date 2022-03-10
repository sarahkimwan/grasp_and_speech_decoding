function [labels_text] = label_number_to_grasp_name(labels_number)
%LABEL_NUMBER_TO_GRASP_NAMES transforms a vector of numbers into a cell of grasp names
label_size = length(labels_number);
labels_text = cell(label_size,1);

for i = 1:label_size
    labels_text{i} = util.image2class_simple(labels_number(i)); 
end 

end

