% Function to calculate the tangent line equation
function [tangentLineEqn, str] = tangentLineFunc(temperature, smoothedHeatFlow, index)
    % Calculate slope at the given index
    slope = calculateSlope(temperature, smoothedHeatFlow, index);
    
    % Get the point of interest (x1, y1)
    x1 = temperature(index);
    y1 = smoothedHeatFlow(index);
    
    % Equation of the tangent line: y = m(x - x1) + y1
    % We'll return the equation in the form of a string
    tangentLineEqn = sprintf('y = %.4f*(x - (%.4f)) + %.8f', slope, x1, y1);

    % Function for plotting
    str = sprintf('%.4f*(x - (%.4f)) + %.8f', slope, x1, y1);
    
    % Display the equation for reference
    %fprintf('Tangent line at index %d (Temperature = %.2f °C, Heat Flow = %.4f W/g):\n', index, x1, y1);
    %fprintf('  %s\n', tangentLineEqn);
end
