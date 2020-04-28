function xyz = MapToXYZ(map)
% convert a nxn distance map to a nx3 xyz table
% input "map" - nxn matrix
% output "xyz" - nx3 matrix
% 
% sample code:
% poly = RandomPolymer(10);
% map = squareform(pdist(poly));
% figure(2); clf; imagesc(map);

n = length(map);
xyz = zeros(n,3);

xyz(1,:) = [0,0,0]; % this explicit notation is easier to read
xyz(2,:) = [map(2,1),0,0]; % first step arbitrarily decided is along x. 
% position 1 and 2 are just declared.
% solve 3, a two-dimensional constraint problem
try
    b = map(2,1)+map(3,2); % dealing with linear interpolation errors
    if map(3,1) > b
        xyz(3,:) = [b,0,0]; 
        map(3,1) = b;
    else
        syms x y;
        eq1 = map(3,1) == sqrt( (x-xyz(1,1)).^2 + (y-xyz(1,2)).^2 );
        eq2 = map(3,2) == sqrt( (x-xyz(2,1)).^2 + (y-xyz(2,2)).^2 );
        a = solve([eq2,eq1],[x,y],'Real',true);
        xyz(3,:) = [eval(a.x(1)),eval(a.y(1)),0]; % 
    end
    % solve m, 3 constraining variables  (could have just started here)
    for m=4:n
        b = map(m,m-1)+ map(m-1,m-2)+ map(m-2,m-3);
        if map(m,1) > b
            xyz(m,:)= [b,0,0];
            map(m,1) = b;
        else            
            syms x y z;
            eq1 = map(m,1) == sqrt( (x-xyz(1,1))^2 + (y-xyz(1,2))^2 + (z-xyz(1,3))^2);
            eq2 = map(m,2) == sqrt( (x-xyz(2,1))^2 + (y-xyz(2,2))^2 + (z-xyz(2,3))^2);
            eq3 = map(m,3) == sqrt( (x-xyz(3,1))^2 + (y-xyz(3,2))^2 + (z-xyz(3,3))^2);
            a = solve([eq1,eq2,eq3],[x,y,z],'Real',true);
            xyz(m,:) = [eval(a.x(1)),eval(a.y(1)),eval(a.z(1))];
        end
    end
catch er
    disp(er.message);
    disp('stop here');
end
% map2 = squareform(pdist(xyz));
% figure(2); clf; 
% subplot(1,3,1); imagesc(map); colorbar;
% subplot(1,3,2); imagesc(map2); colorbar;
% subplot(1,3,3); imagesc(map-map2); colorbar;
% 
% figure(3); clf; 
% subplot(1,2,1); PlotPolymer(poly);
% subplot(1,2,2); PlotPolymer(xyz);

