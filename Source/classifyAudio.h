/*
  ==============================================================================

    classifyAudio.h
    Created: 11 Nov 2020 8:04:20pm
    Author:  JohnK

  ==============================================================================
*/

#include <JuceHeader.h>
#include <cmath>
#include <math.h>
#include "LibrosaMFCC.h"

#pragma once

class classifyAudio
{

    static constexpr auto fftOrder = 10; // The order of the fft; nfft = 2^order
    static constexpr auto fftSize = 1 << fftOrder; // Size of fft in binary
    static constexpr auto hopLength = 768;
    static constexpr auto melFilterNum = 128;
    static constexpr auto dctFilterNum = 18;


public:
    classifyAudio();
    ~classifyAudio();

    std::vector<double> normalizeFeatures(std::vector<double> featureVec);
    std::vector<int> splitAudio(juce::AudioBuffer<float>buffer, std::vector<int>peaks, double sampleRate);
    std::vector<std::vector<float>> doFFT(std::vector<float> audio);
    double freqToMel(double freq);
    double melToFreq(double mel);
    std::vector<std::vector<float>> getMelFilterBank(double sampleRate);
    std::vector<double> linspace(double start_in, double end_in, int num_in);
    std::vector<double> arange(double start_in, double end_in, double spacing);
    std::vector<std::vector<float>> doFilter(std::vector<std::vector<float>> signal_power, std::vector<std::vector<float>> mel_basis);
    std::vector<std::vector<float>> signalPower(std::vector<std::vector<float>> fftData);
    std::vector<std::vector<float>> constructDCT(std::vector<std::vector<float>> signal_filtered);
    std::vector<std::vector<float>> normalize(std::vector<std::vector<float>> weights, std::vector<double> mel_f);
    std::vector<std::vector<float>> dotProduct(std::vector<std::vector<float>> matrix1, std::vector<std::vector<float>> matrix2);
    std::vector<double> meanMfcc(std::vector<std::vector<float>> matrix);
    void tester(juce::AudioBuffer<float> buffer, double sampleRate);
    void testAccuracy2D(std::vector<std::vector<float>> cepCoeff);
    void testAccuracy1D(std::vector<float> section);

private:
    juce::dsp::FFT forwardFFT; // FFT object to perform forward fft on

    double mSampleRate{ 0.0 };
    double pi = 3.1415926535897932385;
    std::vector<double> minVals44100 = {-648.47772217, -41.23974609, -152.00854492, -23.61817169, -61.95689011, -38.56236649, -37.30766678, -67.04342651, -40.2208519, -41.70504379, -31.27592278, -24.66895866, -28.49080849, -37.02992249, -32.1645546, -24.17269897, -21.82624054, -20.05064392};
    std::vector<double> maxVals44100 = {-23.18715858, 237.5743866, 44.8717041, 132.29042053, 69.93689728, 76.80852509, 70.79082489, 48.25970459, 41.93119049, 37.43262863, 36.47388077, 39.49528122, 21.26286507, 28.64596558, 33.27782822, 31.86940575, 21.46536827, 20.09148598};
    std::vector<double> minVals48000 = {-656.27661133, -23.91659355, -156.65887451, -26.94819069, -56.11130905, -55.74367142, -53.92170715, -53.11208725, -53.85018539, -34.81391525, -35.0021019, -21.57038116, -26.25063133, -37.53669357, -30.06274223, -19.2313633, -30.5000782, -21.15400314};
    std::vector<double> maxVals48000 = {-36.69165039, 252.31211853, 44.37088013, 133.45358276, 70.93629456, 79.21143341, 74.84476471, 49.84890747, 40.43525314, 46.55800247, 32.27348709, 34.62980652, 18.92416573, 29.59581566, 28.95998383, 25.31793213, 25.91612244, 23.99448586};
    juce::AudioFormatManager mFormatManager;

    struct svm_model* model44100;
    struct svm_model* model48000;

    Lib_Mfcc Mfcc;
};
