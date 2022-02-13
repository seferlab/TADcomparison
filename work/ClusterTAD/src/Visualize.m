function [ output_args ] = Visualize( Data )
% This function allows for the visualization of data in Heat Map
colormap('hot');
imagesc(Data);
colorbar;
xlabel('Genomic bin (resolution: 40kb)');
ylabel('Genomic bin (resolution: 40kb)');
end

