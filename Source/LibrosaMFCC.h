/*
  ==============================================================================

    LibrosaMFCC.h
    Created: 8 Dec 2020 5:08:20pm
    Author:  JohnK

  ==============================================================================
*/
#pragma once
#include <JuceHeader.h>

class mfcc
{
    static constexpr auto fftSize = 2048;
    static constexpr auto melFilterNum = 128;
    static constexpr float pi = 3.14159265358f;

public:
    mfcc();
    ~mfcc();

    //master function
    std::vector<std::vector<float>> doMfcc(std::vector<float> y, int sampleRate = 22050, int n_mfcc = 20, int hopLength = 512);

    //helper functions
    double freqToMel(double freq);
    double melToFreq(double mel);
    std::vector<double> linspace(double start_in, double end_in, int num_in);
    std::vector<double> arange(double start_in, double end_in, double spacing);
    std::vector<std::vector<float>> normalize(std::vector<std::vector<float>> weights, std::vector<double> mel_f);
    std::vector<std::vector<float>> dotProduct(std::vector<std::vector<float>> matrix1, std::vector<std::vector<float>> matrix2);

    //process functions
    std::vector<float> padAudio(std::vector<float> y);
    std::vector<std::vector<float>> getMelFilterBank(double sampleRate);
    std::vector<std::vector<float>> doFFT(std::vector<float> audio, int hopLength);
    std::vector<std::vector<float>> signalPower(std::vector<std::vector<float>> fftData);
    std::vector<std::vector<float>> doFilter(std::vector<std::vector<float>> signal_power, std::vector<std::vector<float>> mel_basis);
    std::vector<std::vector<float>> doDCT(std::vector<std::vector<float>> signal_filtered, int n_mfcc);

private:

    juce::dsp::FFT forwardFFT; // FFT object to perform forward fft on
};
