function [foundInRange, newStart, newEnd] =  withinRange(vec, startVal, endVal)

vec = sortrows(vec, 1);
flag = false;
newStart = startVal*ones(1,1);
newEnd = endVal*ones(1,1);

for j = 1:length(vec(:,1))
    for h = 1:length(newStart)
    %if the new segement is contained within prev. sequence. 
    if (vec(j,1) <= newStart(h)) && (vec(j,2) >= newEnd(h))
        flag = true;
        newStart(1) = vec(j,2) + 1;
        newEnd(1) = vec(j,2);
    %if the new segment contains old segment
    elseif (vec(j,1) > newStart(h)) && (vec(j,2) < newEnd(h))
        flag = true;
        oldEnd = newEnd(h);
        newEnd(h) = vec(j,1);
        vec(end+1,1) = newStart(h);
        vec(end,2) = vec(j,1);
        [~, addStart, addEnd] = withinRange(vec, vec(j,2), oldEnd);
        for f = 1:length(addStart)
            newStart(h+f,1) = addStart(f);
            newEnd(h+f,1) = addEnd(f);
        end
    %if it is between prexisting addresses:
    elseif (vec(j,1) <= newStart(h)) && (vec(j,2) >= newStart(h)) 
       flag = true;
       if vec(j,2) > newStart
           newStart(1) = vec(j,2);
       end
    elseif (vec(j,1) <= newEnd(h)) && (vec(j,2) >= newEnd(h))
        flag = true;
        if vec(j,1) < newEnd
            newEnd(1) = vec(j,1);
        end
    end
    end
end
foundInRange = flag;

end