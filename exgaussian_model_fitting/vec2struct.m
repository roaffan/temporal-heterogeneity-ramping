%Converts a vector that was passed to the particle swarming function into a
%structure so that labels of each element of the vector become available.
%It uses lb strucre to infer the labels
function [y] = vec2struct(lb,x)

fields=fieldnames(lb);
for i=1:length(x)
    y.(fields{i})=x(i);
end

% % old version with vectors of np, nt and crossterms within the structure 
%function [y] = vec2struct(lb,np,nt,x,with_T)
% %tun vector x inta a structure so we can compute the terms more easily
% fields=fieldnames(lb);
% field_count_vec=1; %index for the lb and ub vectors
% field_count_str=1; %index for the lb and ub structures
% %constant
% y(1).(fields{1})=x(1);
% %postion
% for i=1:6
%     field_count_str=field_count_str+1;
%     for j=1:np
%         field_count_vec=field_count_vec+1;
%         y(j).(fields{field_count_str})=x(field_count_vec);
%     end
% end
% %displacement
% field_count_str=field_count_str+1;
% field_count_vec=field_count_vec+1;
% y(1).(fields{field_count_str})=x(field_count_vec);
% %crosterms
% field_count_str=field_count_str+1;
% if with_T==0
%     field_count_vec=field_count_vec+1;
%     y(1).(fields{field_count_str})=x(field_count_vec);
% else
%     for i=1:4
%         field_count_vec=field_count_vec+1;
%         y(i).(fields{field_count_str})=x(field_count_vec);
%     end
% end
% %time
% if with_T==1
%     for i=1:3
%         field_count_str=field_count_str+1;
%         for j=1:nt
%             field_count_vec=field_count_vec+1;
%             y(j).(fields{field_count_str})=x(field_count_vec);
%         end
%     end
% end
