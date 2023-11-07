function y = Mod(x,q)
    y = mod(x,q); 
    y = y - double(y>q/2)*q;
end


% if ~ismatrix(x)
%     for i=1:length(y)
%         if double(y(i)>q/2)
%             y(i)=y(i)-q;
%         end
%     end
% else
%     for i=1:size(y,1)
%         for j=1:size(y,2)
%             if double(y(i,j)>q/2)
%                 y(i,j)=y(i,j)-q;
%             end
%         end
%     end
% end