% this is a basic demo for our spatial vision model
% In this demo we run the example image through the model and test when our
% model thinks humans should be able to detect a grating against the image
% as a background.

% First, we need to load the image and convert it to luminance:
image = imread('DSC_0339.JPG');
image_bw = convert_lum_spectra(image);
% The conversion here uses measured spectra from a specific monitor used in
% the lab a while ago. If you have conversion methods for your screen you
% can replace the convert function here. The output should be in cd/m^2,
% but rescaling all images with a constant factor will not change model
% predictions.


% Now we can run the early vision model in two settings, with or without
% foveation:
% With foveation is the version we used for all detection experiments. This
% smoothly cuts out a 2x2 dva window around fixation and processes this cut
% out further. You could also use this version for foveal detection
% experiments:
out = early_vision_model(image_bw);
% additional inputs are:
% degSize: size of the image in dva
% timecourse: presentation time in ms, will be used to select parameters
% fixation: position of the fixation in dva relative to the upper left
%           corner.
% pars: parameters as returned by getPars
% sizePx: size in pixels for processing, defaults to 256x256
% Gradidx: which gradients to compute, if any
% refLum: reference luminance if you do not want to use the mean luminance
%         of the input.
% V1Mode: Will switch between the model versions described in the paper.
% csfSelector: allows selection of other CSFs as defined by getCsf in
%              private
out = early_vision_model(image_bw, [2,2], 100);

% the out is a 4D tensor with channel activations of size
% image_x x image_y x n_orientations x n_frequencies
% this will for example show the activations of the second orientation and
% third frequency:
imagesc(out(:,:,2,3))


% We also provide a model version that skips the cut out and processes the
% whole image with constant processing depth. When you want to process
% images as a whole or do not know where subjects are looking this is
% usually more sensible. We can for example handle the image as a 10x10 dva
% image with this function:
out = early_vision_model_NoFovea(image_bw, [10,10], 100);
% the output has the same format as above.

% Both model variants can provide three additional outputs:
% imageNoise: the noise variance for each location & channel
% outGrad: the gradient of the output with respect to the parameters
% noiseGrad: the gradient of the noise with respect to the parameters
%
% imageNoise is necessary for comparisons, the gradients are mostly for
% optimising the parameters.


%% Detection experiment simulation
% to make predictions how well humans con differentiate images, the most
% direct way is to pass 2 images through the model and input the results
% into the compare_images function.

% As an example let's take the image we used above and the image + a
% grating.
grating = image_hanningGrating(size(image_bw),[2,2],7,1,[1,1],0,1,[1,1],0,2)-1;

image1 = image_bw;
image2 = image_bw + 0.1 * mean(image_bw(:)) * grating;

% Run both through the model:
[out1, out1n] = early_vision_model(image1, [2, 2], 100);
[out2, out2n] = early_vision_model(image2, [2, 2], 100);

% Compare:
[diff, noise, p] = compare_images(out1, out1n, out2, out2n);
% this will return the expected difference, its noise variance and the
% probability of detection (without lapses)
% For constant noise, which we usually assume, diff=noise and diff is
% d'^2, because d' = diff/sqrt(noise) = sqrt(diff)
% We did test variants of the model with noise that scales with the signal,
% for which this relationship is no longer true.

% The same comparison can also be made using detection experiment, which
% takes the background image and the signal to be detected as input:
[diff, noise, p] = detection_experiment(grating, image1, [2,2], 100, 0.1);
% This yields a slightly different result here, because detection
% experiment guarantees that the luminances of the combined image stay
% positive. 

% A similar function grating_experiment yields results for comparisons
% between gratings based directly on the grating parameters

% And calc_threshold and calc_threshold_image, use the grating_experiment
% function or the detection_experiment function respectively to test
% different contrasts and arrive at a threshold contrast for a signal
% against a fixed background.
