clc
clear

input = [1 1 3 1 3 3 2 2];
in_perm = perms(input);
sizeof_in_perm = size(in_perm);
perm_row = sizeof_in_perm(1);
row = 0;
col = sizeof_in_perm(2);

for i=1:perm_row
    same_flag = 0;
    for j=i+1:perm_row
       if in_perm(i,:) == in_perm(j,:)
           same_flag = 1;
           break;
       end
    end
    if same_flag == 0
        row = row+1;
        in_perm_final(row,:) = in_perm(i,:);
    end
end
% in_perm_final


for i=1:row
    input_time(i) = 0;
    for j=2:col
        if in_perm_final(i,j-1) == 1 & in_perm_final(i,j) == 1
            input_time(i) = input_time(i) + 79;
        elseif in_perm_final(i,j-1) == 1 & in_perm_final(i,j) == 2
            input_time(i) = input_time(i) + 159;
        elseif in_perm_final(i,j-1) == 1 & in_perm_final(i,j) == 3
            input_time(i) = input_time(i) + 185;
        elseif in_perm_final(i,j-1) == 2 & in_perm_final(i,j) == 1
            input_time(i) = input_time(i) + 79;
        elseif in_perm_final(i,j-1) == 2 & in_perm_final(i,j) == 2
            input_time(i) = input_time(i) + 106;
        elseif in_perm_final(i,j-1) == 2 & in_perm_final(i,j) == 3
            input_time(i) = input_time(i) + 132;
        elseif in_perm_final(i,j-1) == 3 & in_perm_final(i,j) == 1
            input_time(i) = input_time(i) + 79;
        elseif in_perm_final(i,j-1) == 3 & in_perm_final(i,j) == 2
            input_time(i) = input_time(i) + 79;
        elseif in_perm_final(i,j-1) == 3 & in_perm_final(i,j) == 3
            input_time(i) = input_time(i) + 79;
        end
    end
    
    if i==1
        min = input_time(i);
        min_row = i;
    else
        if input_time(i) < min
            min = input_time(i);
            min_row = i;
        end
    end
end

for i=1:row
   if  input_time(i) == min
       in_perm(i,:)
   end
end

% for i=1:perm_row
%     input_time(i) = 0;
%     for j=2:col
%         if in_perm(i,j-1) == 1 & in_perm(i,j) == 1
%             input_time(i) = input_time(i) + 79;
%         elseif in_perm(i,j-1) == 1 & in_perm(i,j) == 2
%             input_time(i) = input_time(i) + 159;
%         elseif in_perm(i,j-1) == 1 & in_perm(i,j) == 3
%             input_time(i) = input_time(i) + 185;
%         elseif in_perm(i,j-1) == 2 & in_perm(i,j) == 1
%             input_time(i) = input_time(i) + 79;
%         elseif in_perm(i,j-1) == 2 & in_perm(i,j) == 2
%             input_time(i) = input_time(i) + 106;
%         elseif in_perm(i,j-1) == 2 & in_perm(i,j) == 3
%             input_time(i) = input_time(i) + 132;
%         elseif in_perm(i,j-1) == 3 & in_perm(i,j) == 1
%             input_time(i) = input_time(i) + 79;
%         elseif in_perm(i,j-1) == 3 & in_perm(i,j) == 2
%             input_time(i) = input_time(i) + 79;
%         elseif in_perm(i,j-1) == 3 & in_perm(i,j) == 3
%             input_time(i) = input_time(i) + 79;
%         end
%     end
%     
%     if i==1
%         min = input_time(i);
%         min_row = i;
%     else
%         if input_time(i) < min
%             min = input_time(i);
%             min_row = i;
%         end
%     end
% end

% for i=1:perm_row
%    if  input_time(i) == min
%        in_perm(i,:)
%    end
% end

% in_perm(min_row,:)
min