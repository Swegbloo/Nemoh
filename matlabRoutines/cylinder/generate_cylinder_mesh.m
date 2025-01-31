function generate_cylinder_mesh(draft, radius_in, radius_out, M)
    % Create or update the outer cylinder files


    for i = 1:M
        angle = (i-0.5) * (360 / M);
        %new_filename = sprintf('%dD', angle);

        % Read the input file
        filename = 'Cylinder_8x.dat';
        data = readmatrix(filename, 'FileType', 'text', 'Delimiter', '\t');
        
        % Get the number of rows and columns
        [rows, cols] = size(data);
        
        % Define the parameters (you need to set these values)
        r_0 = radius_in; % Example value, replace with actual value
        x_0 = radius_out*cos(angle); % Example value, replace with actual value
        y_0 = radius_out*sin(angle); % Example value, replace with actual value
        d_0 = draft; % Example value, replace with actual value
        
        % Process the data
        for j = 2:rows
            % Find the first zero in the row
            zero_index = find(data(j,:) == 0, 1);
            if isempty(zero_index)
                zero_index = cols + 1;
            end
            
            % Perform operations on non-zero elements
            if zero_index > 2
                data(j, 2) = r_0 * data(j, 2) + x_0;
            end
            if zero_index > 3
                data(j, 3) = r_0 * data(j, 3) + y_0;
                data(j, 4) = d_0 * data(j, 4);
            end
        end
        
        % Save the modified data to a new file
        new_filename = sprintf('%dD', angle); % Replace with desired output filename
        writematrix(data, new_filename, 'Delimiter', 'tab');
        % Read the data from the file
        data = readmatrix('filename.txt', 'Delimiter', '\t');
        
        % Replace NaN values in 3rd and 4th columns of 1st row with empty strings
        data(1, 3:4) = {''};
        
        % Open the file for writing
        fid = fopen('filename.txt', 'w');
        
        % Write the data back to the file
        for i = 1:size(data, 1)
            for j = 1:size(data, 2)
                if isempty(data{i,j})
                    fprintf(fid, '\t');
                else
                    fprintf(fid, '%.8f\t', data{i,j});
                end
            end
            fprintf(fid, '\n');
        end
        
        % Close the file
        fclose(fid);
        
    %     fid = fopen(filename, 'w');
    %     fprintf(fid, '2 0\n');
    %     fprintf(fid, '1\n');
    %     fprintf(fid, '1 %.8f %.8f 0\n', draft, radius_in);
    % 
    %     x = radius_out * cosd(angle);
    %     y = radius_out * sind(angle);
    %     fprintf(fid, '%.8f %.8f 0\n', x, y);
    % 
    %     % Generate cylinder mesh points
    %     n_points = 96;
    %     for j = 1:n_points
    %         theta = (j-1) * (360 / n_points);
    %         for k = 1:7
    %             z = -draft * (k-1) / 6;
    %             r = radius_in + x;
    %             px = r * cosd(theta);
    %             py = r * sind(theta) + y;
    %             fprintf(fid, '%.8f %.8f %.8f\n', px, py, z);
    %         end
    %     end
    %     fclose(fid);
    % end

    % Remove any excess cylinder files
    end
