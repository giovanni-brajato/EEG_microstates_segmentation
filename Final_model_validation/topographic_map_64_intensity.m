function [ rendered_image ] = topographic_map_64_intensity( electrodes_signals,max_value,scale_factor, intensity )
% FUNCTION: produce a matrix as large as the default image of the brain,
% with an overlay coloured layer which shows the intensity of the signal
 % INPUT:
 %      -electrodes_signals: 64 values of the electrodes
 %      -max_value: choose in which measure each electrode affect the image
 %      composition
%       - scale_factor: specify the spread of the signal in the map
%       - intensity: intensity of the coloured layer
%  OUTPUT:
%       - rendered_image: image of the brain with signals on top
 %
 % load base image and normalize
base_image = double(imread('Brain_top_view.jpg'));
base_image = base_image./255;
% positions of the electrodes on the image
% the first is x position (horizontal), the second is y (vertical),
% starting from top
electrodes_positions = [105,164; %1
                        130,166; %2
                        154,168; %3
                        181,171; %4
                        208,168; %5
                        233,166; %6
                        258,164; %7
                        101,204; %8
                        128,203; %9
                        154,203; %10
                        182,204; %11
                        209,204; %12
                        236,203; %13
                        262,203; %14
                        106,143; %15
                        130,242; %16
                        155,239; %17
                        182,237; %18
                        209,240; %19
                        234,242; %20
                        258,243; %21
                        151,69; %22
                        183,66; %23
                        215,69; %24
                        115,92; %25
                        145,100; %26
                        182,99; %27
                        219,99; %28
                        247,93; %29
                        96,120; %30
                        115,127; %31
                        138,132; %32
                        159,132; %33
                        185,135; %34
                        203,131; %35
                        227,132; %36
                        250,128; %37
                        268,121; %38
                        81,159; %39
                        283,160; %40
                        75,203; %41
                        289,204; %42
                        50,203; %43
                        316,204; %44
                        80,248; %45
                        284,248; %46
                        95,285; %47
                        115,279; %48
                        139,276; %49
                        160,275; %50
                        181,272; %51
                        203,274; %52
                        227,274; %53
                        249,280; %54
                        268,285; %55
                        118,316; %56
                        145,308; %57
                        183,310; %58
                        219,308; %59
                        247,314; %60
                        152,337; %61
                        183,344; %62
                        214,336; %63
                        183,383; %64
                        ];

% how should we characterize the level? we may thing as a gaussian
% interpolation. The signal level will determine how big is the gaussian(larger and taller)

x_dimension = 1:size(base_image,2);
y_dimension = 1:size(base_image,1);
[x,y] = meshgrid(x_dimension,y_dimension);
% initialize overlay image
overlay_image = zeros(size(base_image,1),size(base_image,2));
% we are interested only in the magnitude
electrodes_signals = abs(electrodes_signals);
for s = 1:length(electrodes_signals)
mean_x = electrodes_positions(s,1);
mean_y =  electrodes_positions(s,2);
sigma = abs((electrodes_signals(s)))*scale_factor;
levels = 1/(sqrt(2*pi*sigma^2))*max_value*exp(-(1/2*(x - mean_x).^2/sigma^2 + 1/2*(y -mean_y).^2/sigma^2));
overlay_image = overlay_image + levels; 
end
%  wanna black as the lowest color
personal_jet = jet(255);
personal_jet(1:32,:) = personal_jet(1:32,:).*repmat((logspace(-5,0,32))',1,3);
% image(overlay_image)
jet_overlay = ind2rgb(gray2ind(overlay_image,255),personal_jet);

rendered_image = base_image + intensity*jet_overlay;


end

