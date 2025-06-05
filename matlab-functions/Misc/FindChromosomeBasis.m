function [ChrBasis, gof] = FindChromosomeBasis(poly, NaNoption, angleThreshold)
    % poly is an (nx3) matrix. 
    % where row i corresponds to an xyz coordinate of Barcode i
    
    % if the user does not specify NaNoption 
    % set it to ignore NaN
    if nargin == 1
        NaNoption = 0;
        angleThreshold = 0;
    elseif nargin == 2
        angleThreshold = 0;
    end

    % if NaNoption = 0, this function ignores all NaN 
    % only work on barcodes that have coordinates 
    if NaNoption == 0
        nanHybs = any(isnan(poly), 2);
        poly = poly(~nanHybs, :);
    % else, this function uses linear interpolation to 
    % infer coordinates of the missing barcodes.
    else 
        poly = fillmissing(poly, 'linear', 1);
    end 
    % disp('poly')
    % disp(poly)

    % now we are finding the vector that connects between
    % barcode i and barcode i+1 
    dirVec = diff(poly);
    dirVec = dirVec ./ vecnorm(dirVec, 2, 2);
    
    % disp('dirVec')
    % disp(dirVec)

    % now guess that each vector is best describing the 
    % chromosome axis.
    % the vector that best describing the chromosome axis 
    % will maximize abs(the sum of sgn(v_i dot v_j)) for
    % all i and j
    numSameDirMat = zeros(height(dirVec), 1); 
    for i = 1:height(dirVec)
        vecI = dirVec(i, :);
        direction = vecI * dirVec';
        % disp('direction')
        % disp(direction)
        numSameDir = sum(direction >= angleThreshold); 
        numSameDirMat(i) = numSameDir; 
    end 
    
    [maxSameDir, bestBasisIndex] = max(numSameDirMat);
    % goodness of fit is the number of maxSameDir / total number of vectors
    gof = numSameDirMat(bestBasisIndex)/height(dirVec);

    bestBasis = dirVec(bestBasisIndex, :)';
    % then normalize to make it a unit vector
    bestBasis = bestBasis./norm(bestBasis);
    
    % disp('bestBasis')
    % disp(bestBasis)
    % disp('gof')
    % disp(gof)

    % now we are going to perturb this by a bit until 
    % it converges to some data
    % perturb by the range = 1 degree 
    stepSize = 45; 
    
    % need to fix perturbation 
    % of the vectors 

    % set maxIteration in case the function does not 
    % converge 
    maxIteration = 20; 
    
    % if poly is in 2D space then 
    % rotation matrix is R(2x2) 
    Rotation2D = @(x) ([cosd(x) -sind(x);
                        sind(x) cosd(x)]);
    
    if width(poly) == 2 
        % iterate for maxIteration times 
        for iter = 1:maxIteration 
            % generate 10 random angle
            randomAngle = -stepSize + 2*stepSize.*rand(10, 1);
            numSameDirMat = zeros(height(randomAngle), 1);

            for iAngle = 1:size(randomAngle)
                currAngle = randomAngle(iAngle);
                
                perturbVec = Rotation2D(currAngle)*bestBasis;
                
                % disp('perturbVec')
                % disp(perturbVec)
                % disp('dirVec')
                % disp(dirVec)

                direction = dirVec * perturbVec;
                numSameDir = sum(direction > angleThreshold);

                numSameDirMat(iAngle) = numSameDir; 
            end 

            betterVec = numSameDirMat > maxSameDir; 
            % disp('betterVec')
            % disp(betterVec)
            if sum(betterVec) > 0
                [maxSameDir, bestBasisIndex] = max(numSameDirMat);
                bestBasis = Rotation2D(randomAngle(bestBasisIndex))*bestBasis;
                gof = numSameDirMat(bestBasisIndex)/height(dirVec);
            else 
                break 
            end 
        end 
        
        % then the second orthonormal basis, by convention 
        % is achieved by rotating the first basis by 90 deg CCW 
        secondBasis = Rotation2D(90) * bestBasis;
        % disp('secondBasis')
        % disp(secondBasis)
        ChrBasis = [bestBasis secondBasis]; 
    % 3D case
    elseif width(poly) == 3
        % iterate for maxIteration times 

        rotx = @(x) ([1 0 0; ...
                      0 cosd(x) -sind(x); ...
                      0 sind(x) cosd(x)]);

        roty = @(x) ([cosd(x) 0 sind(x); ...
                      0 1 0; ...
                      -sind(x) 0 cosd(x)]);

        rotz = @(x) ([cosd(x) -sind(x) 0; ...
                      sind(x) cosd(x) 0; ...
                      0 0 1]);

        for iter = 1:maxIteration 
            % generate 10 random angle
            randomAngle = -stepSize + 2*stepSize.*rand(100, 3);
            numSameDirMat = zeros(height(randomAngle), 1);
           
            for iAngle = 1:size(randomAngle)
                currAngle = randomAngle(iAngle, :);
                
                Rx = rotx(currAngle(1));
                Ry = roty(currAngle(2));
                Rz = rotz(currAngle(3));
                
                perturbVec = Rz*Ry*Rx*bestBasis;
                
                % disp('perturbVec')
                % disp(perturbVec)
                % disp('dirVec')
                % disp(dirVec)

                direction = dirVec * perturbVec;
                numSameDir = sum(direction > angleThreshold);

                numSameDirMat(iAngle) = numSameDir; 
            end 

            betterVec = numSameDirMat > maxSameDir; 
            % disp('betterVec')
            % disp(betterVec)
            if sum(betterVec) > 0
                % disp('improving')
                % disp(iter)
                [maxSameDir, bestBasisIndex] = max(numSameDirMat);
                bestBasis = Rz*Ry*Rx*bestBasis;
                gof = numSameDirMat(bestBasisIndex)/height(dirVec);
            end 
        end 

        % guess the second basis 
        secondBasis = rotz(90)*bestBasis;
        
        % guess the third basis 
        thirdBasis = roty(90)*secondBasis;

        % then perform Gram-Schmidt diagonalization
        % to find orthonormal basis
        ChrBasis = GramSchmidt([bestBasis secondBasis thirdBasis]);

        % make sure that this is a proper orthogonal matrix
        % det(Q) = 1, right hand basis
        % if det(Q) = -1 switch 2nd basis and 3rd basis
        if abs(det(ChrBasis) + 1) < 0.0001
            ChrBasis = [ChrBasis(:, 1) ChrBasis(:, 3) ChrBasis(:, 2)];
        end 
    end
    
    

    


     
    