%% Reaction Order Calculator given Endothermic Decomposition
clear; clc; close all;

x = (-2:.05:2)';
Poly = (x+1).*(x+2).*(x-1).*(x-1);

data = horzcat(x,Poly);

% Load DSC data
%data = load('Me_When_I_finally_Get_to_tuch_the_DSC'); % DSC OFF LIMITS :(
temperature = data(:, 1);
heatFlow = data(:, 2);

% Smooth data to reduce noise
smoothedHeatFlow = smooth(heatFlow, 5);

% Use the findpeaks function to detect peaks
[pks, locs, widths, proms] = findpeaks(-smoothedHeatFlow, temperature,'MinPeakProminence', 0.1);

pks = -pks;     % pks - Height of detected peakrs // locs - Index of peaks

% Second derivative for finding points of inflection
secondDerivative = gradient(gradient(smoothedHeatFlow));

% Empty Data Sets
inflectionPoints = zeros(length(pks), 2);
S = zeros(length(pks),1);
n = zeros(length(pks),1);

% find points of inflection on left and right of peak
for i = 1:length(pks)
    % Get the location of the current peak
    peakIndex = find(temperature == locs(i));

    % Find the left side inflection point
    j = 1;
    while true
        leftInflectionIndex = (peakIndex - j);
        if secondDerivative(leftInflectionIndex) <= 0
            break
        elseif peakIndex - j == 1
            break
        else
            j = j+1;
        end
    end

    % Find the right side inflection point
    j = 1;
    while true
        rightInflectionIndex = (peakIndex + j);
        if secondDerivative(rightInflectionIndex) <= 0
            break
        elseif height(x) - peakIndex == j
            break
        else
            j = j+1;
        end
    end

    % Store the left and right points of inflection in the matrix
    inflectionPoints(i, 1) = temperature(leftInflectionIndex);
    inflectionPoints(i, 2) = temperature(rightInflectionIndex);

    % Left inflection point is the max left we can go
    % Right inflection point is the max right we can go

    % Recursion where we scan entire right side then use left + 1 and repeat
    % Create a empty matrix to store S1 & S2 data 
    S1Matrix = zeros(peakIndex-(leftInflectionIndex+1),(rightInflectionIndex-1)-peakIndex);
    S2Matrix = S1Matrix;

    for k = leftInflectionIndex+1:peakIndex
        % For each index from left inflection point to peakIndex, calculate tangent line
        % Calculate tangent line for point k
        index = k;
        tangentlineEqn = tangentLineFunc(temperature, smoothedHeatFlow, index);
        leftTangentLineEqn = tangentlineEqn;
    
        for l = peakIndex:rightInflectionIndex-1
            % Calculate the shape ratio S1 between second derivatives
            % Calculate the shape ratio S1 between second derivatives
            rowIndex = k - (leftInflectionIndex + 1) + 1;  % Adjusted index
            colIndex = l - peakIndex + 1;  % Adjusted index
            S1Matrix(rowIndex, colIndex) = secondDerivative(k) / secondDerivative(l);

            % Calculate tangent line for point l
            index = l;
            tangentlineEqn = tangentLineFunc(temperature, smoothedHeatFlow, index); 
            rightTangentLineEqn = tangentlineEqn;

            % Calculate the intersection point
            [intercept] = lineIntersection(leftTangentLineEqn, rightTangentLineEqn);

            % Difference in Temperature
            B = intercept(1) - temperature(k);
            A = temperature(l) - intercept(1);

            % Shape ratio given distances
            S2Matrix(rowIndex, colIndex) = A./B; % Store the computed S2 value = A./B;
        end
    end
    
    % Find the most accurate S value 
    minDistance = 0.5;  % Set a minimum distance (or temperature difference)

    DifferenceS = abs(S2Matrix - S1Matrix);
    minDifference = min(DifferenceS(:));
    [minRow, minCol] = find(DifferenceS == minDifference);

    % Ensure that the points are at least 'minDistance' apart
    while abs(temperature(minCol + peakIndex) - temperature(minRow + leftInflectionIndex)) < minDistance
        % Remove this pair and find the next closest match
        DifferenceS(minRow, minCol) = inf;
        minDifference = min(DifferenceS(:));
        [minRow, minCol] = find(DifferenceS == minDifference);
    end
    nR = length(minRow);


    if mod(nR, 2) == 1 % If the number of elements is odd, take the middle value
        medianRow = minRow((nR+1)/2);
        medianCol = minCol((nR+1)/2);
    else
        medianCol = minCol(nR/2 + 1); % If even, take the upper middle value
    end

    S(i) = (S1Matrix(medianRow, medianCol) + S2Matrix(medianRow, medianCol))./2;
    n(i) = 1.26*S(i).^1.26;

    % Track back to find the actual temperature indices
    actualRowIndex = medianRow + (leftInflectionIndex + 1); % Convert minRow to the actual index in temperature
    actualColIndex = medianCol + peakIndex; % Convert minCol to the actual index in temperature
    % Convert to symbolic expressions
    [tangentLineL, strL] = tangentLineFunc(temperature, smoothedHeatFlow, actualRowIndex);
    [tangentLineR, strR] = tangentLineFunc(temperature, smoothedHeatFlow, actualColIndex);
    fL = str2func(['@(x)',strL]);
    fR = str2func(['@(x)',strR]);
    % Calculate the intersection point
    [intercept] = lineIntersection(tangentLineL, tangentLineR);

    % Plot the DSC data to visualize the peak
    figure;
    plot(temperature, smoothedHeatFlow);
    xlabel('Temperature (°C)');
    ylabel('Heat Flow (W/g)');
    title('DSC Curve');
    hold on;
    % Highlight tangent lines
    fplot(fL,[temperature(1),temperature(end)])
    fplot(fR,[temperature(1),temperature(end)])
    % Highlight the tangent points
    plot(temperature(actualRowIndex), smoothedHeatFlow(actualRowIndex), 'ro', 'MarkerFaceColor', 'r', 'DisplayName', 'Point on Left');
    plot(temperature(actualColIndex), smoothedHeatFlow(actualColIndex), 'go', 'MarkerFaceColor', 'g', 'DisplayName', 'Point on Right');
    %plot(intercept,'bo')

end

