% simpleParforTest.m
function simpleParforTest()

    % Define the size of the array
    N = 1e3; % Adjust as needed for your test
    
    % Generate an array of random numbers
    data = rand(N, 1);
    
    % Preallocate the result array
    result = zeros(size(data));
    
    % Start a timer
    tic;
    
    % Calculate the square of each element in parallel
    parfor idx = 1:N
        result(idx) = data(idx)^2;
    end
    
    % Stop the timer and display elapsed time
    elapsedTime = toc;
    fprintf('Elapsed time: %.3f seconds\n', elapsedTime);
    
    % Optionally, display some of the results for verification
    disp('First 10 results:');
    disp(result(1:10));
end
