% Function to calculate the slope at a given index
function slope = calculateSlope(temperature, smoothedHeatFlow, index)
    % Handle edge cases (first or last points)
    if index == 1
        % Forward difference for the first point
        slope = (smoothedHeatFlow(index+1) - smoothedHeatFlow(index)) / (temperature(index+1) - temperature(index));
    elseif index == length(temperature)
        % Backward difference for the last point
        slope = (smoothedHeatFlow(index) - smoothedHeatFlow(index-1)) / (temperature(index) - temperature(index-1));
    else
        % Central difference for points in the middle
        slope = (smoothedHeatFlow(index+1) - smoothedHeatFlow(index-1)) / (temperature(index+1) - temperature(index-1));
    end
end