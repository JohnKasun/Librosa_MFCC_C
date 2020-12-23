/*
  ==============================================================================

    onsetDetection.h
    Created: 14 Oct 2020 2:38:45pm
    Author:  Walter Kopacz

  ==============================================================================
*/

#include "classifyAudio.h"
#include <cmath>
#include <math.h>

#pragma once
class onsetDetection
{

public:

    static constexpr auto fftOrder = 10; // The order of the fft; nfft = 2^order
    static constexpr auto fftSize = 1 << fftOrder; // Size of fft in binary
    static constexpr auto hopLength = 768;

    onsetDetection();
    ~onsetDetection();
    
    std::vector<float> makeNoveltyFunction(juce::AudioBuffer<float>buffer, int audioNumOfSamples, double sampleRate);
    std::vector<int> pickPeaks(std::vector<float>noveltyFunction);
    std::vector<float> convertIndiciesToTime(std::vector<int>peaksInIndicies);
    std::vector<int> convertTimeToSamples(std::vector<float>peaksInTime);
    void testSegmentation(std::vector<float>noveltyFunction, std::vector<int> peaks, std::vector<float> peaksInSeconds);
    
    
    /*std::vector<std::vector<float>>fftData;
    * std::vector<float> noveltyFunction
    std::vector<int>peaks;
    std::vector<float>peaksInSeconds;
    std::vector<float>peaksInSamples;*/
    

private:
    
    double mSampleRate{ 0.0 };
   
    juce::dsp::FFT forwardFFT; // FFT object to perform forward fft on
    juce::dsp::WindowingFunction<float> hannWindow; 
    classifyAudio classification;
        
    
};
